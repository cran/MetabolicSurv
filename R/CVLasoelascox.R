#'   Cross Validations for Lasso Elastic Net Survival predictive models and Classification
#'
#' The function does cross validation for Lasso, Elastic net and Ridge regressions models before the survial analysis and classification. The survival analysis is based on the selected metabolites in the presence or absene of prognostic factors.
#'
#' The function performs the cross validations for Lasso, Elastic net and Ridge regressions models for Cox proportional hazard model. Metabolites are selected at each iteration and then use for the classifier. Which implies that predictive metabolites signature is varied from one cross validation to the other depending on selection. The underline idea is to investigate the Hazard Ratio for the train and test data based on the optimal lambda selected for the non-zero shrinkage coefficients, the nonzero selected metabolites will thus be used in the survival analysis and in calculation of the risk scores for each sets of data.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Censor A vector of censoring indicator
#' @param Mdata A large or small metabolic profile matrix. A matrix with metabolic profiles where the number of rows should be equal to the number of metabolites and number of columns should be equal to number of patients.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Quantile The cut off value for the classifier, default is the median cutoff
#' @param Metlist A list of metabolites to be considered in the model usually smaller than the metabolites in the Mdata . Default is to use all metabolites available and it is advisable to be greater than 17.
#' @param Standardize A Logical flag for the standardization of the metabolite matrix, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is standardize=TRUE.
#' @param Reduce A boolean parameter indicating if the metabolic profile matrix should be reduced, default is TRUE and larger metabolic profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of metabolites (default is 15) to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE
#' @param Alpha The mixing parameter for glmnet (see \code{\link[glmnet]{glmnet}}). The range is 0<= Alpha <= 1. The Default is 1
#' @param Fold number of folds to be used for the cross validation. Its value ranges between 3 and the numbe rof subjects in the dataset
#' @param Ncv Number of validations to be carried out. The default is 25.
#' @param nlambda The number of lambda values - default is 100 as in glmnet.
#' @return A object of class \code{\link[MetabolicSurv]{cvle}} is returned with the following values
#'   \item{Coef.mat}{A matrix of coefficients with rows equals to number of cross validations and columns equals to number of metabolites.}
#'   \item{Runtime}{A vector of runtime for each iteration measured in seconds.}
#'   \item{lambda}{A vector of estimated optimum lambda for each iterations.}
#'   \item{n}{A vector of the number of selected metabolites}
#'   \item{HRTrain}{A matrix of survival information for the training dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.}
#'   \item{HRTest}{A matrix of survival information for the test dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.}
#'   \item{pld}{A vector of partial likelihood deviance at each cross validations.}
#'   \item{Met.mat}{A matrix with 0 and 1. Number of rows equals to number of iterations and number of columns equals to number of metabolites. 1 indicates that the particular metabolite was selected or had nonzero coefficient and otherwise it is zero.}
#'  \item{Mdata}{The Metabolite data matrix that was used for the analysis either same as Mdata or a reduced version.}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}},
#' \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[glmnet]{glmnet}}, \code{\link[MetabolicSurv]{Lasoelacox}}
#' @examples
#' \donttest{
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#' Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## USING THE FUNCTION
#' Results = CVLasoelacox(Survival = Data$Survival,Censor = Data$Censor,
#' Mdata = t(Data$Mdata),Prognostic = Data$Prognostic, Quantile = 0.5,
#' Metlist = NULL,Standardize = TRUE, Reduce=FALSE, Select=15,
#' Alpha = 1,Fold = 4,Ncv = 10,nlambda = 100)
#'
#' ## NUMBER OF SELECTED METABOLITES PER CV
#' Results@n
#'
#' ## GET THE MATRIX OF COEFFICIENTS
#' Results@Coef.mat
#'
#' ## SURVIVAL INFORMATION OF THE TRAIN DATASET
#' Results@HRTrain
#'
#' ## SURVIVAL INFORMATION OF THE TEST DATASET
#' Results@HRTest
#' }

#' @export CVLasoelacox
#' @import superpc
#' @import stats
#' @import methods
#' @import grDevices
#' @import graphics

CVLasoelacox <- function (Survival,
                        Censor,
                        Mdata,
                        Prognostic,
                        Quantile = 0.5,
                        Metlist = NULL,
                        Standardize = TRUE,
                        Reduce=TRUE,
                        Select=15,
                        Alpha = 1,
                        Fold = 4,
                        Ncv = 10,
                        nlambda = 100)
{

  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Mdata)) stop("Argument 'Mdata' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")

  if(is.null(Metlist)){
    Mdata <- Mdata
  } else {
    Mdata <- Mdata[is.element(rownames(Mdata),Metlist),]
  }


  if (Reduce) {
    Dat <- Mdata
    DataForReduction<-list(x=Dat,y=Survival, censoring.status=Censor, metabolitenames=rownames(Dat))
    TentativeList<-names(sort(abs(superpc::superpc.train(DataForReduction, type="survival")$feature.scores),decreasing =TRUE))[1:Select]
    Mdata<-Dat[TentativeList,]
  } else {
    Mdata<-Mdata
  }



  metnames <- rownames(Mdata)

  n.mets<-nrow(Mdata)
  n.patients<-ncol(Mdata)
  Run.Time<-rep(NA, Ncv)
  n.train<-(n.patients-floor(n.patients/Fold))
  n.test<-floor(n.patients/Fold)
  cv.train <-matrix(0,Ncv,n.train)
  cv.test  <-matrix(0,Ncv,n.test)
  n.g<-rep(0, Ncv)

  #optimum lambda
  lambda<-rep(NA, Ncv)
  pld<-rep(NA, Ncv)

  #HR--------
  HRTrain<-matrix(NA,nrow=Ncv,ncol=3)
  HRTest<-matrix(NA,nrow=Ncv,ncol=3)



  if(is.null(Prognostic)){
    Data <- Mdata
    Data.Full <- t(Mdata)
    Penalty <- rep(1,row(Data.Full))
  } else{
    Data <- t(Mdata)
    Data.Full <- cbind(Data,Prognostic)
    Penalty <- c(rep(1,ncol(Data)),rep(0,ncol(Prognostic)))
  }

  # Survival times must be larger than 0
  Survival[Survival <= 0] <- stats::quantile(Survival, probs = 0.01)
  perPrognostic<-NULL
  coef.mat<-met.mat<-matrix(0,nrow=Ncv,ncol=nrow(Mdata))

  pIndex <- c(1:n.patients)

  for (i in 1:Ncv){
    message('Cross validation loop ',i)

    cv.train[i,] <-sort(sample(pIndex,n.train,replace=F) )
    cv.test[i,] <-c(1:n.patients)[-c(intersect(cv.train[i,] ,c(1:n.patients)))]

    Stime=Survival[cv.train[i,]]
    sen= Censor[cv.train[i,]]
    Data.Full2 <-Data.Full[cv.train[i,],]

    Run.Time[i] <-system.time(Lasso.Cox.CV <- glmnet::cv.glmnet(x = as.matrix(Data.Full2),y = survival::Surv(as.vector(Stime),as.vector(sen) == 1),
                            family = 'cox',
                            alpha = Alpha,
                            nfolds= Fold,
                            nlambda = nlambda,
                            penalty.factor = Penalty,
                            standardize = Standardize))[3]

  # Results of the cv.glmnet procedure

  Lambda <- Lasso.Cox.CV$lambda.min
  Alllamda <- Lasso.Cox.CV$lambda
  Coefficients <- stats::coef(Lasso.Cox.CV, s =Lambda)
  Coefficients.NonZero <- Coefficients[Coefficients[, 1] != 0, ]


  pld[i]<-Lasso.Cox.CV$cvm[Lasso.Cox.CV$lambda==Lasso.Cox.CV$lambda.min]  # partial likelihood deviance

  if (!is.null(dim(Coefficients.NonZero)))
  {
    Selected.mets <- setdiff(colnames(Data.Full2),colnames(Prognostic))
    Coefficients.NonZero <- Coefficients[c(Selected.mets,colnames(Prognostic)),]
    Select <- "F"
  }else{
    Select <- "T"
  }

  if (is.null(Prognostic))
  {
    Selected.mets <- names(Coefficients.NonZero)
  }else{
    Selected.mets <- setdiff(names(Coefficients.NonZero), colnames(Prognostic))
  }

  n.g[i] <- length(Selected.mets)

  # What to do if no metabolite is selected?
  # Decrease lambda until a metabolite is selected
  # Going through the list of lambda values and repeat the above procedure

  Lambda <- Lasso.Cox.CV$lambda.min
  Lambda.Sequence <- Lasso.Cox.CV$lambda
  Lambda.Index <- which(Lambda.Sequence == Lasso.Cox.CV$lambda.min)

  Lambda.Index.Add = 0

  while (n.g[i] <= 0) {
    Lambda.Index.Add <- Lambda.Index.Add + 1
    Coefficients <- stats::coef(Lasso.Cox.CV, s = 0.12)
    Coefficients.NonZero <- Coefficients[Coefficients[, 1] != 0,]

    if (!is.null(dim(Coefficients.NonZero)))
    {
      Selected.mets <- setdiff(colnames(Data.Full),colnames(Prognostic))
      Coefficients.NonZero <- Coefficients[c(Selected.mets,colnames(Prognostic)),]
      Select <- "F"
    }else{
      Select <- "T"
    }

    if (is.null(Prognostic))
    {
      Selected.mets <- names(Coefficients.NonZero)
    }
    else
    {
      Selected.mets <- setdiff(names(Coefficients.NonZero), colnames(Prognostic))
    }
    Lambda <- Lambda.Sequence[Lambda.Index + Lambda.Index.Add]
    n.g[i] <- length(Selected.mets)
  }

  lambda[i]<-Lambda
  met.mat[i,is.element(rownames(Mdata),Selected.mets)]<-1
  coef.mat[i,is.element(rownames(Mdata),Selected.mets)]<-Coefficients.NonZero[Selected.mets]

  scores.train <- as.vector(Coefficients.NonZero[Selected.mets] %*% t(Data[cv.train[i,],Selected.mets]))
  scores.test <- as.vector(Coefficients.NonZero[Selected.mets] %*% t(Data[cv.test[i,],Selected.mets]))


  #######################
  ## train set ###########
  Sdata<-data.frame(Survival=Survival[cv.train[i,]],Censor=Censor[cv.train[i,]])

  if (!is.null(Prognostic)) perPrognostic<-as.data.frame(Prognostic[cv.train[i,],])
colnames(perPrognostic) <- colnames(Prognostic)

  Results1<-EstimateHR(Risk.Scores=scores.train, Data.Survival =Sdata, Prognostic = perPrognostic, Plots = FALSE, Quantile = Quantile)

  if (!is.null(Prognostic)) HRTrain[i,]<-(summary(Results1$SurvResult)[[8]])[1,][-2]
  if ( is.null(Prognostic)) HRTrain[i,]<- summary(Results1$SurvResult)[[8]][-2]

  #######################
  ## test set ###########
  Sdata<-data.frame(Survival=Survival[cv.test[i,]],Censor=Censor[cv.test[i,]])

  if (!is.null(Prognostic)) perPrognostic<-as.data.frame(Prognostic[cv.test[i,],])
  colnames(perPrognostic) <- colnames(Prognostic)

  Results2<-EstimateHR(Risk.Scores =scores.test, Data.Survival = Sdata, Prognostic = perPrognostic, Plots = FALSE, Quantile = Quantile)

  if (!is.null(Prognostic)) HRTest[i,]<-(summary(Results2$SurvResult)[[8]])[1,][-2]
  if ( is.null(Prognostic)) HRTest[i,]<- summary(Results2$SurvResult)[[8]][-2]


}

  return(methods::new("cvle",Coef.mat=coef.mat,Runtime=Run.Time,lambda=lambda,n=n.g,Met.mat=met.mat,HRTrain=HRTrain,HRTest=HRTest,pld=pld,Mdata=Mdata))
}
