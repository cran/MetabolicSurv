#' Generate Artificial Metabolic Survival Data
#'
#' The Function generates metabolic profile of any number of patients and also their survival information.
#'
#' The function generates the metabolic profile where small set of metabolites (30) are informative and rest of them are set as noisy metabolites. Next to that Survival time and Censoring information are generated based on first right singular vectors of \code{\link[base]{svd}} of the metabolic profile matrix. It also generates other prognostic factors such as Age, Stage and Gender which are slightly correlated with survival time.
#' @param nPatients The number of patients
#' @param nMet The number of metabolites
#' @param Prop The proportion of patients having low risk
#' @return An object of class list is returned with the following items .
#'   \item{Censor}{The censoring/event indicator}
#'   \item{Survival}{The Survival time}
#'   \item{Met.names}{The vector of metabolites}
#'   \item{Mdata}{The metabolic profile matrix}
#'   \item{Prognostic}{A data frame with prognostic factors.}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}}, \code{\link[MetabolicSurv]{DataHR}}
#' @examples
#' #GENERATE SOME METABOLIC SURVIVAL DATA WITH PROGNOSTIC FACTORS
#'
#' Data<-MSData(nPatients=100,nMet=150,Prop=0.5)
#'
#' SurvTime<-Data$Survival
#' Censor<-Data$Censor
#' ProgFact<-Data$Prognostic
#' MetData<-Data$Mdata
#' Metnames<-Data$Met.names
#' @export MSData

MSData<-function(nPatients=100,nMet=150,Prop=0.5){

  set.seed(425)
  Mdata<-matrix(rnorm(nMet*nPatients),ncol=nPatients)
  Mdata[1:10,]<-matrix(rnorm(10*nPatients,sd=1,mean=3),ncol=nPatients)
  Mdata[11:30,]<-matrix(rnorm(20*nPatients,sd=1,mean=2),ncol=nPatients)

  Survival<-c(abs(10+svd(Mdata[1:floor(nPatients*Prop),])$v[1:floor(nPatients*Prop),1]+
                    2*rnorm(floor(nPatients*Prop))),
              abs(8 +svd(Mdata[1:floor(nPatients*Prop),])$v[(floor(nPatients*Prop)+1):nPatients,1]+
                    3*rnorm(nPatients-floor(nPatients*Prop))))

  CenTime<-abs(10+svd(Mdata[1:floor(nPatients),])$v[1:floor(nPatients),1]+
                 1*rnorm(floor(nPatients)))

  Censor<- ifelse(Survival>CenTime,0,1)

  metabolitenames <- paste("Met",as.character(1:nMet),sep="")
  rownames(Mdata)<-metabolitenames

  ProgFact<-data.frame(Age=floor(Survival*0.68+rnorm(nPatients,30,10)),
                       Stage=sample(1:4,nPatients,replace=T),Gender=rbinom(nPatients, 1, 0.5))

  Mdata = t(Mdata)
  return(list(Censor=Censor,Survival=Survival,Met.names=metabolitenames,Mdata=Mdata,Prognostic=ProgFact))
}
