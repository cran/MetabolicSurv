#'Cross validation for the Metabolite specific analysis
#'
#' The function performs cross validation for each metabolite depending the number of fold which guides the division into the train and testing dataset. The classifier is then obtained on the training dataset to be validated on the test dataset
#'
#' This function performs the cross validation for metabolite by metabolite analysis. The data will firstly be divided into data train dataset and test datset. Furthermore, a metabolite-specific model is fitted on train data and a classifier is built. In addition, the classifier is then evaluated on test dataset for each particular metabolite. The Process is repeated for all the full or reduced metabolites to obtaind the HR statistics of the low risk group. The following steps depends on the number of cross validation specified.
#' @param Fold Number of times in which the dataset is divided. Default is 3 which implies dataset will be divided into three groups and 2/3 of the dataset will be the train datset and 1/3 will be to train the results.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Mdata A large or small metabolic profile matrix. A matrix with metabolic profiles where the number of rows should be equal to the number of metabolites and number of columns should be equal to number of patients.
#' @param Censor A vector of censoring indicator
#' @param Reduce A boolean parameter indicating if the metabolic profile matrix should be reduced, default is TRUE and larger metabolic profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of metabolites (default is 15) to be selected from supervised PCA. This is valid only if th argument Reduce=TRUE
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Quantile The cut off value for the classifier, default is the median cutoff
#' @param Ncv The Number of cross validation loop. Default is 50 but it is recommended to have at least 100.
#' @return A object of class \code{\link[MetabolicSurv]{cvmm}} is returned with the following values
#'   \item{HRTrain}{The Train dataset HR statistics for each metabolite by the number of CV}
#'   \item{HRTest}{The Test dataset HR statistics for each metabolite by the number of CV}
#'   \item{train}{The selected subjects for each CV in the train dataset}
#'   \item{train}{The selected subjects for each CV in the test dataset}
#' \item{n.mets}{The number of metabolite used in the analysis}
#'  \item{Ncv}{The number of cross validation performed}
#'  \item{Rdata}{The Metabolite data matrix that was used for the analysis either same as Mdata or a reduced version.}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}},
#' \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[MetabolicSurv]{MSpecificCoxPh}},
#' @examples
#' \donttest{
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#' Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## USING THE FUNCTION
#' Result = CVMetSpecificCoxPh(Fold=3,Survival=Data$Survival,
#' Mdata=t(Data$Mdata),Censor= Data$Censor,Reduce=TRUE,
#' Select=150,Prognostic=Data$Prognostic,Quantile = 0.5,Ncv=3)
#'
#' ## GET THE CLASS OF THE OBJECT
#' class(Result)     # An "cvmm" Class
#'
#' ##  METHOD THAT CAN BE USED FOR THE RESULT
#' show(Result)
#' summary(Result)
#' plot(Result)
#' }

#' @export CVMetSpecificCoxPh

CVMetSpecificCoxPh<-function(Fold=3,
                             Survival,
                             Mdata,
                             Censor,
                             Reduce=TRUE,
                             Select=150,
                             Prognostic=NULL,
                             Quantile = 0.5,
                             Ncv=3)
{
  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Mdata)) stop("Argument 'Mdata' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")

  if (Reduce) {
    DataForReduction<-list(x=Mdata,y=Survival, censoring.status=Censor, metabolitenames=rownames(Mdata))
    TentativeList<-names(sort(abs(superpc.train(DataForReduction, type="survival")$feature.scores),decreasing =TRUE))[1:Select]

    ReduMdata<-Mdata[TentativeList,]
  } else {
    ReduMdata<-Mdata
  }

  n.mets<-nrow(ReduMdata)
  n.patients<-ncol(ReduMdata)
  gNames<-rownames(ReduMdata)

  n.train<-(n.patients-floor(n.patients/Fold))
  n.test<-floor(n.patients/Fold)
  train <-matrix(0,Ncv,n.train)
  test  <-matrix(0,Ncv,n.test)

  #     THE HAZARD RATIO
  HRTrain<-array(NA,dim=c(n.mets,4,Ncv))
  HRTest<-array(NA,dim=c(n.mets,4,Ncv))

  set.seed(123)
  pIndex <- c(1:n.patients)
  options(warn=-1)


  Onemetcox<-function(imet,Prognostic,Survival,Censor,index){

    if (is.null(Prognostic)) {

      scdata <- data.frame(Survival=Survival[index],Censor=Censor[index],imet=imet[index])
      model1 <- survival::coxph(Surv(Survival, Censor==1) ~ imet,data=scdata)
    }
    if (!is.null(Prognostic)) {
      if (is.data.frame(Prognostic)) {
        nprog<-ncol(Prognostic)
        scdata <- data.frame(Survival=Survival[index],Censor=Censor[index],imet=imet[index],Prognostic[index,])
        prognames<-colnames(Prognostic)
        eval(parse(text=paste( "model1 <-survival::coxph(Surv(Survival, Censor==1) ~ imet",paste("+",prognames[1:nprog],sep="",collapse =""),",data=scdata)" ,sep="")))
      } else {
        stop(" Argument 'Prognostic' is NOT a data frame ")
      }

    }

    return(list(model1=model1,imet=imet,scdata=scdata))
  }


  for (j in 1:Ncv){
    message('Cross validation loop ',j)
    train[j,] <-sort(sample(pIndex,n.train,replace=F))
    test[j,] <-c(1:n.patients)[-c(intersect(train[j,] ,c(1:n.patients)))]

    for (i in 1:n.mets){  #---------------------------  STRAT FOR LOOP ------------------

      #training set ---------------------------------------------------------
      TrainTemp<-Onemetcox(ReduMdata[i,train[j,]],Prognostic,Survival,Censor,train[j,])

      m1 <- TrainTemp$model1
      Trtandmetabolite<-summary(m1)[[7]][c("imet"),1]
      p1.train     <- Trtandmetabolite[1]*ReduMdata[i,train[j,]]
      p1.test<- Trtandmetabolite[1]*ReduMdata[i,test[j,]]

      Tempimet <-EstimateHR(p1.train, Data.Survival = data.frame(Survival=Survival[train[j,]],Censor=Censor[train[j,]]), Prognostic = Prognostic[train[j,],],Plots = FALSE, Quantile = Quantile)

      HRTrain[i,,j]<-summary(Tempimet$SurvResult)[[8]][1,]

      #testing set ---------------------------------------------------------

      TempimetTE <-EstimateHR(p1.test,Data.Survival=data.frame(Survival=Survival[test[j,]],Censor=Censor[test[j,]]), Prognostic=Prognostic[test[j,],],Plots = FALSE, Quantile = Quantile)
      HRTest[i,,j]<-summary(TempimetTE$SurvResult)[[8]][1,]

    }#---------------------------  END OF  FOR LOOP for metabolite ----------------------

  }#--------------------------END OF loop over CV ---------------------------------------

  return(new("cvmm",HRTrain=HRTrain,HRTest=HRTest,train=train,test=test,n.mets=n.mets,Ncv=Ncv,Rdata = ReduMdata))
}
