#'Classifiction for Majority Votes
#'
#' The Function fits cox proportional hazard model and does classification based on the majority votes.
#'
#' The Function fits cox proportional hazard model and does classification based on the majority votes while estimating the Hazard ratio of the low risk group. The function firstly count the number of Low risk classification for each subject based on the metabolite specific analysis which determ,ines the majority votes. In addition, It visualizes the metabolic specific calssification for the subjects. 25 ssubjects is taken for visualization purpose.
#' @param Result An object obtained from the metabolite specific analysis (\code{\link[MetabolicSurv]{MSpecificCoxPh}}) which is of class "ms"
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Censor A vector of censoring indicator
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param J The jth set of patients required for the visualization. The default is J=1 which is the first set of patients. For visualization, J should be less than the number of patients divided by 25
#' @return A list is returned with the following values
#'   \item{Model.result}{The cox proportional regression result based on the majority vote classification}
#'   \item{N}{The majority vote for each subject}
#'   \item{Classif}{The majority vote classification for each subjects}
#'   \item{Group}{The classification of the subjects based on each metabolite analysis}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{MSpecificCoxPh}}, \code{\link[survival]{coxph}},  \code{\link[MetabolicSurv]{EstimateHR}}
#' @references
#' \insertRef{tib}{MetabolicSurv}
#' @examples
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#'Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## RUNNING THE METABOLITE SPECIFIC FUNCTION
#' Example1 = MSpecificCoxPh(Survival = Data$Survival,
#'Mdata = t(Data$Mdata), Censor = Data$Censor, Reduce = FALSE,
#'Select = 15,Prognostic = Data$Prognostic, Quantile = 0.5)
#'
#' ## USING THE FUNCTION
#' Result2 = Majorityvotes(Example1,Data$Prognostic, Data$Survival,Data$Censor,J=2)
#'
#' ## THE SURVIVAL ANALYSIS FOR MAJORITY VOTE RESULT
#'  Result2$Model.result
#'
#'### THE MAJORITY VOTE FOR EACH SUBJECT
#' Result2$N
#'
#'### THE MAJORITY VOTE CLASSIFICATION FOR EACH SUBJECT
#' Result2$Classif
#'
#'### THE GROUP FOR EACH SUBJECT BASED ON THE METABOLITE SPPECIFIC ANALYSIS
#' Result2$Group
#' @export Majorityvotes


Majorityvotes<-function(Result,Prognostic, Survival,Censor,J=1){


  if (class(Result)!="ms") stop("Invalid class object.")
  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")


  ggr <- per.R<-per.NR <- NULL


  Group = Result@Group
  Group = data.frame(Group)

   for (i in 1:length(Survival)){
    per.R[i]<-sum(Group[,i]=="Low risk")
    ggr[i]<-ifelse((nrow(Group)-per.R[i])>per.R[i],"High risk","Low risk")
  }

  ggr<-as.factor(ggr)

  if (is.null(Prognostic)) {

    cdata <- data.frame(Survival,Censor,ggr)
    m0 <- survival::coxph(Surv(Survival, Censor==1) ~ ggr,data=cdata)
  }

  if (!is.null(Prognostic)) {
    if (is.data.frame(Prognostic)) {
      nProg<-ncol(Prognostic)
      cdata <- data.frame(Survival,Censor,ggr,Prognostic)
      NameProg<-colnames(Prognostic)
      eval(parse(text=paste( "m0 <-survival::coxph(Surv(Survival, Censor==1) ~ ggr",paste("+",NameProg[1:nProg],sep="",collapse =""),",data=cdata)" ,sep="")))
    } else {

      stop(" Argument 'Prognostic' is NOT a data frame ")
    }

  }

  VoteMat<-Group
  ng<-nrow(VoteMat)
  np<-ncol(VoteMat)
  Jmax<-floor(np/25)
  if (J<Jmax) {
    slist<-(1+(J-1)*25):(J*25)
    par(mar=c(5.1, 4.1, 4.1, 8.1), xpd=TRUE)
    plot(0,0, xlim=c(-2,ng),ylim=c(1,25),xlab="Metabolites",ylab="Patient Index",axes=FALSE,type="n", main="Metabolite-Specific Classification\n of 25  selected patients",cex.main=0.9)
    for(i in 1:ng){
      points(rep(i, 25)[VoteMat[i,slist]=="Low risk"],which(VoteMat[i,slist]=="Low risk"),col=4,pch=15)
      points(rep(i, 25)[VoteMat[i,slist]=="High risk"],which(VoteMat[i,slist]=="High risk"),col=5,pch=15)
    }
    axis(1);axis(2,at=1:25,slist);box()
    legend("topright",inset=c(-0.15,0),c("Low Risk","High Risk"),pch=c(15,15),col=c(4,5))
  } else {stop("J should be less than (no.of patients / 25)")
  }

  return(list(Model.result=m0, N= per.R, Classif=ggr,Group=VoteMat))

}
