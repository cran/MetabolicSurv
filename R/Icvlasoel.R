#'   Inner and Outer Cross Validations for Lasso Elastic Net Survival predictive models and Classification
#'
#' The function does cross validation for Lasso, Elastic net and Ridge regressions models based on fixed or top selected metabolites from \code{\link[MetabolicSurv]{CVLasoelacox}} with classifier validated on a independent sample for the survial analysis and classification. The survival analysis is based on the selected metabolites in the presence or absene of prognostic factors.
#'
#' The function does cross validation for Lasso, Elastic net and Ridge regressions models based on fixed or top selected metabolites from \code{\link[MetabolicSurv]{CVLasoelacox}} with classifier validated on a independent sample for the survial analysis and classification. The survival analysis is based on the selected metabolites in the presence or absene of prognostic factors. The classifier is built on the weights obtain from the inner cross validations results and it is tested on out-of-bag data. These weights can be fixed or can be updated at each outer iteration. If weights are not fixed then patients are classified using majority votes. Otherwise, weights obtained from the inner cross validations are summarized by mean weights and used in the classifier. Inner cross validations are performed by calling to function \code{\link[MetabolicSurv]{CVLasoelacox}}. Hazard ratio for low risk group is estimated using out-of-bag data.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Censor A vector of censoring indicator
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Mdata A large or small metabolic profile matrix. A matrix with metabolic profiles where the number of rows should be equal to the number of metabolites and number of columns should be equal to number of patients.
#' @param Fold number of folds to be used for the cross validation. Its value ranges between 3 and the numbe rof subjects in the dataset
#' @param Ncv Number of validations to be carried out. The default is 25.
#' @param Nicv Number of validations to be carried out for the inner loop. The default is 5.
#' @param Alpha The mixing parameter for glmnet (see \code{\link[glmnet]{glmnet}}). The range is 0<= Alpha <= 1. The Default is 1
#' @param TopK Top list of metabolites. Usually this can be mostly selected metabolites by function \code{\link[MetabolicSurv]{CVLasoelacox}}.
#' @param Weights A logical flag indicating if a fixed or non-fixed weights should be used during the classifier evaluations. Default is FALSE.
#' @return A object of class \code{\link[MetabolicSurv]{fcv}} is returned with the following values
#'   \item{Runtime}{A vector of runtime for each iteration measured in seconds.}
#'   \item{Fold}{Number of folds used.}
#'   \item{Ncv}{Number of outer cross validations used.}
#'   \item{Nicv}{Number of inner cross validations used.}
#'   \item{TopK}{The Top metabolites used}
#'   \item{HRInner}{A 3-way array in which first, second, and third dimensions correspond to Nicv, 1, and Ncv respectively. This contains estimated HR for low risk group on the out of bag data.}
#'   \item{HRTest}{A matrix of survival information for the test dataset based on the out of bag data. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.}
#'  \item{Weight}{A matrix with columns equals number of TopK metabolites and rows Ncv. Note that Weights are estimated as colMeans of coefficients matrix return from the inner cross validations.}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{CVLasoelacox}},
#' \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[glmnet]{glmnet}}, \code{\link[MetabolicSurv]{Lasoelacox}}
#' @examples
#' \donttest{
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#' Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## USING THE FUNCTION
#' Results = Icvlasoel(Data$Survival, Data$Censor, Data$Prognostic,
#' t(Data$Mdata), Fold = 3,Ncv = 5, Nicv = 7, Alpha = 1,
#' TopK = colnames(Data$Mdata[,80:100]), Weights = FALSE)
#'
#' ## NUMBER OF Outer CV
#' Results@Ncv
#' ## NUMBER OF Inner CV
#' Results@Nicv
#'
#' ## HR of low risk group for the Inner CV
#' Results@HRInner
#'
#' ## HR of low risk group for the out of bag dataset
#' Results@HRTest
#'
#' ## The weight for the analysis
#' Results@Weight
#' }

#' @export Icvlasoel
#' @import tidyr
#' @import rms
#' @import matrixStats
#' @import Rdpack

Icvlasoel<-function(Survival, Censor,Prognostic=NULL,Mdata,Fold=3,Ncv=50,Nicv=100,Alpha=0.1
                    ,TopK,Weights=FALSE
                    )
{

  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Mdata)) stop("Argument 'Mdata' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")

  TopK<-TopK
  if (length(TopK)!=sum(is.element(rownames(Mdata),TopK))) stop("TopK should be a subset of Mdata")

  #TopK<-rownames(ReduMdata)

  options(warn=-1)
  t1 <- proc.time()
  n.patients<-length(Survival)
  n.train<-(n.patients-floor(n.patients/Fold))
  n.test<-floor(n.patients/Fold)
  ind.train <-matrix(0,Ncv,n.train)
  ind.test  <-matrix(0,Ncv,n.test)

  HRInner<-array(NA,c(Nicv,1,Ncv))
  WeightMat<-matrix(NA,ncol=length(TopK),nrow=Ncv)

  #HR--------

  HRTE<-matrix(NA,nrow=Ncv,ncol=3)

  #----------record majority vote results
  MajorV<-matrix(NA,nrow=n.test,ncol=Ncv)



  #------------------------  Begin FOR LOOP :1 ---------------------------
  pIndex <- c(1:n.patients)
  Survival[Survival<=0]<-quantile(Survival,probs =0.02)

  for (j in 1:Ncv){
    ind.train[j,] <-sort(sample(pIndex,n.train,replace=F) )
    ind.test[j,] <-c(1:n.patients)[-c(intersect(ind.train[j,] ,c(1:n.patients)))]
  }



  for (j in 1:Ncv){
    message(paste(" ", j, " Outer Iteration ",sep=""))
    perProg<-NULL
    if (!is.null(Prognostic))     perProg<-Prognostic[ind.train[j,],]
    #inner cross valdiation
    ResT<-CVLasoelacox(Survival=Survival[ind.train[j,]],Censor=Censor[ind.train[j,]],
                       Mdata[is.element(rownames(Mdata),TopK),ind.train[j,]],Prognostic=perProg,Quantile=0.5,
                       Metlist = NULL, Standardize = TRUE, Reduce = FALSE, Select = 15,
                       Alpha = 1, Fold = 4, Ncv = Nicv, nlambda = 100)

    HRInner[,1,j]<-ResT@HRTest[,1]

    MeanWeightMet<-sapply(1:length(TopK),function(i) mean(ResT@Coef.mat[,i]))
    names(MeanWeightMet)<-TopK
    WeightMat[j,]<-MeanWeightMet
    #end of inner cross validation

    #weights are NOT fixed

    if (!Weights) {
      MajorityVotes<-matrix(NA,ncol=length(ind.test[j,]),nrow=Nicv)
      for (ii in 1:Nicv){
        #message('Inner Cross validation loop ',ii)
        WeightsInner<-ResT@Coef.mat[ii,]
        names(WeightsInner)<-TopK
        # Risk score for Test set
        ScoreInner<- WeightsInner%*%Mdata[is.element(rownames(Mdata),TopK),ind.test[j,]]

        MedianScoreInner<-median(ScoreInner)
        MajorityVotes[ii,]<-ifelse(ScoreInner<=MedianScoreInner,1,0)  # 1 => Responsive, 0 => non responsive
      }

      ResNres<-   sapply(1:length(ind.test[j,]),function(kk) sum(MajorityVotes[,kk]))<(Nicv/2) # patients should get more than 50% votes to be classified as Responsive
      MajorV[,j]<-ResNres

      #---------------------- Create indicator variable ----------------
      ##HR low
      gid<- NULL
      n.patients<-n.test
      gid[c(1:n.patients)[ResNres==FALSE]]<- "low"
      gid[c(1:n.patients)[ResNres==TRUE]]<- "high"



      #extract out of bag sample-----------------------------------------------------------------------------

      if (is.null(Prognostic)) {

        cdata <- data.frame(Survival=Survival[ind.test[j,]],Censor=Censor[ind.test[j,]],gid=as.factor(gid))
        m1 <- coxph(Surv(Survival, Censor==1) ~ gid,data=cdata)
      }

      if (!is.null(Prognostic)) {
        if (is.data.frame(Prognostic)) {
          nPrgFac<-ncol(Prognostic)
          cdata <- data.frame(Survival=Survival[ind.test[j,]],Censor=Censor[ind.test[j,]],gid=as.factor(gid),Prognostic[ind.test[j,],])
          NameProg<-colnames(Prognostic)
          eval(parse(text=paste( "m1 <-coxph(Surv(Survival, Censor==1) ~ gid",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
        } else {

          stop(" Argument 'Prognostic' is NOT a data frame ")
        }

      }

      HRTE[j,]<-(summary(m1)[[8]][1,])[-2]
    }  else {


      #-------------------------  now weights are fixed ----------------------------------------
      # Risk score for Test set

      RiskScore<-MeanWeightMet%*%Mdata[is.element(rownames(Mdata),TopK),ind.test[j,]]
      gid<- NULL
      n.patients<-n.patients<-n.test
      gid[c(1:n.patients)[RiskScore < median(RiskScore)]]<- "low"
      gid[c(1:n.patients)[RiskScore >= median(RiskScore)]]<- "high"

      if (is.null(Prognostic)) {

        cdata <- data.frame(Survival=Survival[ind.test[j,]],Censor=Censor[ind.test[j,]],gid=as.factor(gid))
        m1 <- coxph(Surv(Survival, Censor==1) ~ gid,data=cdata)
      }

      if (!is.null(Prognostic)) {
        if (is.data.frame(Prognostic)) {
          nPrgFac<-ncol(Prognostic)
          cdata <- data.frame(Survival=Survival[ind.test[j,]],Censor=Censor[ind.test[j,]],gid=as.factor(gid),Prognostic[ind.test[j,],])
          NameProg<-colnames(Prognostic)
          eval(parse(text=paste( "m1 <-coxph(Surv(Survival, Censor==1) ~ gid",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
        } else {

          stop(" Argument 'Prognostic' is NOT a data frame ")
        }

      }
      HRTE[j,]<-(summary(m1)[[8]][1,])[-2]

    }
  }

  t2 <- proc.time()
  RunTime<-as.vector(t2-t1)[3]
  #----------- results -------------
  return(new("fcv",Runtime=RunTime, Fold=Fold,  Ncv=Ncv,  Nicv=Nicv,TopK=TopK,
             HRInner=HRInner,HRTest=HRTE,Weight=WeightMat))
}
#----------------------------------------------------------- END
