#'Metabolite by metabolite Cox proportional analysis
#'
#' The Function fits cox proportional hazard model and does classification for each metabolite
#'
#' This function fits  metabolite by metabolite Cox proportional hazard model and perform the classification based on a speciied quantile risk score which has been estimated using a single  metabolite. Function is useful for majority vote classification method and  metabolite by  metabolite analysis and also for top K  metabolites.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Mdata A large or small metabolic profile matrix. A matrix with metabolic profiles where the number of rows should be equal to the number of metabolites and number of columns should be equal to number of patients.
#' @param Censor A vector of censoring indicator
#' @param Reduce A boolean parameter indicating if the metabolic profile matrix should be reduced, default is TRUE and larger metabolic profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of metabolites (default is 15) to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Quantile The cut off value for the classifier, default is the median cutoff
#' @return A object of class \code{\link[MetabolicSurv]{ms}} is returned with the following values
#'   \item{Result}{The cox proportional regression result for each metabolite}
#'   \item{HRRG}{The hazard ratio statistics (Hazard-ratio, Lower confidence interval and upper confidence interval) of the riskgroup based on the riskscore and the cut off value for each metabolite}
#'   \item{Group}{The classification of the subjects based on each metabolite analysis}
#'   \item{Metnames}{The names of the metabolites for the analysis}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}},  \code{\link[MetabolicSurv]{EstimateHR}}
#' @examples
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#'Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## USING THE FUNCTION
#'Example1 = MSpecificCoxPh(Survival = Data$Survival,
#'Mdata = t(Data$Mdata), Censor = Data$Censor, Reduce = FALSE,
#'Select = 15,Prognostic = Data$Prognostic, Quantile = 0.5)
#'
#' ## KNOWLING THE CLASS OF THE OUTPUT
#'class(Example1)
#'
#' ## EXTRACTING THE COMPONENT OF THE FUNCTION
#' ### HAZARD RATIO INFORMATION FOR EACH METABOLITES
#'Example1@HRRG
#'
#'### COX MODEL RESULT FOR EACH METABOLITES
#'Example1@Result
#'
#'### CLASSIFICATION FOR EACH METABOLITES
#'Example1@Group
#' @export MSpecificCoxPh

MSpecificCoxPh<-function(Survival,
                         Mdata,
                         Censor,
                         Reduce=FALSE,
                         Select=15,
                         Prognostic=NULL,
                         Quantile = 0.5){

  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Mdata)) stop("Argument 'Mdata' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")

  if (Reduce) {
    DataForReduction<-list(x=Mdata,y=Survival, censoring.status=Censor, metabolitenames=rownames(Mdata))
    TentativeList<-names(sort(abs(superpc.train(DataForReduction, type="survival")$feature.scores),decreasing =TRUE))[1:Select]
    TentativeList

    ReduMdata<-Mdata[TentativeList,]
  } else {
    ReduMdata<-Mdata
  }


  n.metabolite<-nrow(ReduMdata)
  n.patients<-ncol(ReduMdata)
  gNames<-rownames(ReduMdata)

  HRp <- matrix(0,n.metabolite,3)
  gr <- matrix(0, n.metabolite,n.patients)
  res <-  vector("list", n.metabolite)

  for (i in 1:n.metabolite){  #---------------------------  STRAT FOR LOOP ------------------------

    meti <- ReduMdata[i,]

    if (is.null(Prognostic)) {

      cdata <- data.frame(Survival,Censor,meti)
      m0 <- survival::coxph(survival::Surv(Survival, Censor==1) ~ meti,data=cdata)
    }
    if (!is.null(Prognostic)) {
      if (is.data.frame(Prognostic)) {
        nPrgFac<-ncol(Prognostic)
        cdata <- data.frame(Survival,Censor,meti,Prognostic)
        NameProg<-colnames(Prognostic)
        eval(parse(text=paste( "m0 <-survival::coxph(survival::Surv(Survival, Censor==1) ~ meti",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
      } else {

        stop(" Argument 'Prognostic' is NOT a data frame ")
      }

    }
    #risk Score
    TrtandPC1<-summary(m0)[[7]][c("meti"),1]
    p1 <- TrtandPC1*meti

    Tempmeti <-EstimateHR(Risk.Scores =p1,Data.Survival =cdata, Prognostic=Prognostic,Plots = FALSE, Quantile = Quantile )
    res[[i]]<-Tempmeti$SurvResult
    HRp[i,]<-summary(Tempmeti$SurvResult)[[8]][1,c(1,3,4)]
    gr[i,]<-Tempmeti$Riskgroup

  }#---------------------------  END OF  FOR LOOP ------------------------

  colnames(HRp) = c("HR","LowerCI","UpperCI")
  Metnames = rownames(as.matrix(ReduMdata))
  rownames(gr) = Metnames
  colnames(gr) = paste0("Subject", 1: ncol(ReduMdata))

  return(new("ms",Result=res,HRRG=HRp,Group=gr,Metnames = Metnames))
}
