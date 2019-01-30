#'Survival PCA and Classification for metabolic data
#'
#' The function performs principal component analysis (PCA) on Metabolomics matrix and fit Cox proportional hazard model with covariates using also the first PCA as covariates.
#'
#' This function can also be used to perform the grid analysis where the grid will be several quantile values and default is 0.5 which is the median cut-off. This function can handle single and multiple metabolites. For larger Metabolomics matrix, this function will reduce largerMetabolomics matrix to smaller version using supervised pca approach and this is by default done and can be control by using the argument Reduce. Other prognostic factors can be included to the model.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Mdata A large or small metabolic profile matrix. A matrix with metabolic profiles where the number of rows should be equal to the number of metabolites and number of columns should be equal to number of patients.
#' @param Censor A vector of censoring indicator
#' @param Reduce A boolean parameter indicating if the metabolic profile matrix should be reduced, default is TRUE and larger metabolic profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of metabolites (default is 15) to be selected from supervised PCA. This is valid only if th argument Reduce=TRUE
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Plots A boolean parameter indicating if the plots should be shown. Default is FALSE
#' @param Quantile The cut off value for the classifier, default is the median cutoff
#' @return A object of class SurvPca is returned with the following values
#'   \item{Survfit}{The cox proportional regression result using the first PCA}
#'   \item{Riskscores}{A vector of risk scores which is equal to the number of patents.}
#'   \item{Riskgroup}{The classification of the subjects based on the PCA into low or high risk group}
#'   \item{pc1}{The First PCA scores based on either the reduced Metabolite matrix or the full matrix}
#' \item{KMplot}{The Kaplan-Meier survival plot of the riskgroup}
#'  \item{SurvBPlot}{The distribution of the survival in the riskgroup}
#'  \item{Riskpca}{The plot of  Risk scores vs first PCA}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}},
#' \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[stats]{princomp}},
#'  \code{\link[MetabolicSurv]{SurvPlsClass}}
#' @references
#' \insertRef{ye1}{MetabolicSurv}
#' @examples
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#'Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## USING THE FUNCTION
#'Result = SurvPcaClass(Survival=Data$Survival, Mdata=t(Data$Mdata),
#'Censor=Data$Censor, Reduce = FALSE, Select = 150,
#'Prognostic = Data$Prognostic, Plots = FALSE, Quantile = 0.5)
#'
#' ## GETTING THE SURVIVAL REGRESSION OUTPUT
#' Result$SurvFit
#'
#' ## GETTING THE RISKSCORES
#'Result$Riskscores
#'
#' ### GETTING THE RISKGROUP
#'Result$Riskgroup
#'
#'### OBTAINING THE FIRST PRINCIPAL COMPONENT SCORES
#'Result$pc1
#' @export SurvPcaClass


SurvPcaClass<-function(
  Survival,
  Mdata,
  Censor,
  Reduce=TRUE,
  Select=150,
  Prognostic=NULL,
  Plots = FALSE,
  Quantile = 0.5

){

  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Mdata)) stop("Argument 'Mdata' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")


  if (Reduce) {
    DataForReduction<-list(x=Mdata,y=Survival, censoring.status=Censor, metnames=rownames(Mdata))
    TentativeList<-names(sort(abs(superpc::superpc.train(DataForReduction, type="survival")$feature.scores),decreasing =TRUE))[1:Select]
    TentativeList

    ReduMdata<-Mdata[TentativeList,]
  } else {
    ReduMdata<-Mdata
  }


  n.mets<-nrow(ReduMdata)
  n.patients<-ncol(ReduMdata)

  if (is.matrix(ReduMdata)) {
    pc1 <- f.pca(as.matrix(ReduMdata))[[6]][,1]
  } else {
    pc1<-ReduMdata
  }

  if (is.null(Prognostic)) {

    cdata <- data.frame(Survival,Censor,pc1)
    m0 <- survival::coxph(Surv(Survival, Censor==1) ~ pc1,data=cdata)
  }
  if (!is.null(Prognostic)) {
    if (is.data.frame(Prognostic)) {
      nPrgFac<-ncol(Prognostic)
      cdata <- data.frame(Survival,Censor,pc1,Prognostic)
      NameProg<-colnames(Prognostic)
      eval(parse(text=paste( "m0 <-survival::coxph(survival::Surv(Survival, Censor==1) ~ pc1",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
    } else {

      stop(" Argument 'Prognostic' is NOT a data frame ")
    }

  }
  Riskscores <- Riskgroup <- NULL

  #risk Score
  TrtandPC1<-summary(m0)[[7]][c("pc1"),1]
  p1 <- TrtandPC1*pc1

  TempRes<- EstimateHR(Risk.Scores = p1, Data.Survival = cdata, Prognostic = Prognostic, Plots = TRUE, Quantile = Quantile)
  gg <- data.frame(Riskscores = p1,Riskgroup = TempRes$Riskgroup,pc1 = pc1)
  ab <- ggplot2::ggplot(gg, ggplot2::aes(x=Riskscores, y=pc1, shape=Riskgroup, color=Riskgroup)) + ggplot2::geom_point()

  tempp<-list(SurvFit=TempRes$SurvResult,Riskscores = p1, Riskgroup=TempRes$Riskgroup,pc1=pc1)
  class(tempp)<-"SurvPca"

  temp<-list(SurvFit=TempRes$SurvResult,Riskscores = p1, Riskgroup=TempRes$Riskgroup,pc1=pc1, KMplot = TempRes$KMplot, SurvBPlot = TempRes$SurvBPlot, Riskpca = ab)
  class(temp)<-"SurvPca"

  if (Plots){
    return(temp)
  } else{
    return(tempp)
  }
}
