#'Cross validation for sequentially increases metabolites
#'
#' This function does cross validation for the metabolite by metabolite analysis while sequentially increasing the number of metabolites as specified.
#'
#' The function is a cross validation version of the function \code{\link[MetabolicSurv]{SIMet}}. This function firstly processes the cross validation for the metabolite by metabolite  analysis results, and then sequentially considers top k metabolites. The function recompute first PCA or PLS on train data and estimate risk scores on both test and train data only on the metabolite matrix with top k metabolites. Patients are then classified as having low or high risk based on the test data where the cutoff used is median of the risk score. The process is repeated for each top K metabolite sets.
#' @param Object An object of class \code{\link[MetabolicSurv]{cvmm}}
#' @param Top The Top k number of metabolites to be used
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Censor A vector of censoring indicator
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @return A object of class \code{\link[MetabolicSurv]{cvsim}} is returned with the following values
#'    \item{HRpca}{A 3-way array in which first, second, and third dimensions correspond to number of metabolites, Hazard ratio infromation(Estimated HR, LowerCI and UpperCI), and number of cross validation respectively. This contains the estimated HR on test data and dimension reduction method is PCA.}
#'    \item{HRpls}{A 3-way array in which first, second, and third dimensions correspond to number of metabolites, Hazard ratio infromation(Estimated HR, LowerCI and UpperCI), and number of cross validation respectively. This contains the estimated HR on test data and dimension reduction method is PLS.}
#'    \item{Nmets}{The number of metabolites in the reduced matrix}
#'   \item{Ncv}{The number of cross validation done}
#'   \item{Top}{A sequence of top k metabolites considered. Default is Top=seq(5,100,by=5)}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{MSpecificCoxPh}}, \code{\link[MetabolicSurv]{SIMet}}
#' @examples
#' \donttest{
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#' Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## GETTING THE cvmm OBJECT
#' Result = CVMetSpecificCoxPh(Fold=3,Survival=Data$Survival,
#' Mdata=t(Data$Mdata),Censor= Data$Censor,Reduce=TRUE,Select=150,
#' Prognostic=Data$Prognostic,Quantile = 0.5,Ncv=3)
#'
#' ## USING THE FUNCTION
#'  Result2 = CVSim(Result, Top = seq(5, 100, by = 5), Data$Survival,
#'  Data$Censor,Prognostic = Data$Prognostic)
#'
#' ## GET THE CLASS OF THE OBJECT
#' class(Result2)     # An "cvsim" Class
#' }
#' @export CVSim

CVSim<-function(Object,Top=seq(5,100,by=5),Survival,Censor, Prognostic=NULL){

   Decrease=FALSE
  if (class(Object)!="cvmm") stop("Invalid class object.")
  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")


  met.mat<-matrix(0,Object@Ncv,Object@n.mets)

  Mdata<-Object@Rdata
  met.names <- rownames(Mdata)

  n.mets<-Object@n.mets
  n.patients<-ncol(Mdata)


  if (Object@n.mets<max(Top)) stop("The max(Top) should be less than or equal to total number of metabolites")


  HRPC<-array(NA,dim=c(Object@Ncv,3,length(Top)))
  HRPL<-array(NA,dim=c(Object@Ncv,3,length(Top)))

  #i=1
  options( warn = -1)
  for (j in 1: length(Top)){
    message('Analysis for ',paste0("Top",Top[j]))
    for (i in 1:Object@Ncv){
      message('Cross validation loop ',j)
      hrsetTE<-Object@HRTest[,c(1,3,4),i]
      hrsetmT<-Object@HRTrain[,c(1,3,4),i]

      Names.KMetsT<-rownames(Mdata)[order(hrsetTE[,1],decreasing=Decrease)]
      index.Top.KMetsT <- order(hrsetTE[,1],decreasing=Decrease)
      index.Top.KMets<-index.Top.KMetsT[1:Top[j]]

      meti <-Mdata[intersect(rownames(Mdata),Names.KMetsT[1:Top[j]]),]

      #PCA ----------------

      TrainTemp<-IntermediatePCA(meti,Prognostic,Survival,Censor,Object@train[i,])
      TestTemp<- IntermediatePCA(meti,Prognostic,Survival,Censor,Object@test[i,])
      m1 <- TrainTemp$m0
      TrtandMet1<-summary(m1)[[7]][c("pc1"),1]
      rm(m1)
      p1.train     <- TrtandMet1[1]*TrainTemp$pc1
      p1.test     <-  TrtandMet1[1]*TestTemp$pc1

      TempmetiTE <-EstimateHR(p1.test,Data.Survival=data.frame(Survival=Survival[Object@test[i,]],Censor=Censor[Object@test[i,]]),Prognostic=Prognostic[Object@test[i,],],Plots = FALSE, Quantile = 0.5)

      HRPC[i,,j]<- (summary(TempmetiTE$SurvResult)[[8]][1,])[-2]

      #PLS -------------------
      TrainTemp1<-IntermediatePLS(meti,Prognostic,Survival,Censor,Object@train[i,])
      TestTemp1 <-IntermediatePLS(meti,Prognostic,Survival,Censor,Object@test[i,])
      m1 <- TrainTemp1$m0
      TrtandMet2<-summary(m1)[[7]][c("pc1"),1]
      rm(m1)
      p1.train     <- TrtandMet2[1]*TrainTemp1$pc1
      p1.test     <- TrtandMet2[1]*TestTemp1$pc1
      TempmetiTE <-EstimateHR(p1.test,Data.Survival=data.frame(Survival=Survival[Object@test[i,]],Censor=Censor[Object@test[i,]]),Prognostic=Prognostic[Object@test[i,],],Plots = FALSE, Quantile = 0.5)

    HRPL[i,,j]<- (summary( TempmetiTE$SurvResult)[[8]][1,])[-2]

    }

  }

  Ncv<-Object@Ncv
  n.mets<-Object@n.mets
  return(new("cvsim",HRpca=HRPC,HRpls=HRPL,Nmets=n.mets,Ncv=Ncv,Top=Top))
  }
