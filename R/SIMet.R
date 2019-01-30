#'Sequential Increase in Metabolites for the PCA or PLS classifier
#'
#' The Function fits cox proportional hazard model and does classification by sequentially increasing the metabolites  using either PCA or PLS based on the topK metabolites specified.
#'
#' This function sequentially increase the number of top K metabolites to be used in the PCA or PLS methods in order to obtain the risk score. This function internally calls \code{\link[MetabolicSurv]{MSpecificCoxPh}} to rank the metabolites based on HR for each metabolite. Therefore metabolites can be ordered based on increasing order of the HR for low risk group. Thereafter, the function takes few top K (15 is the default) to be used in the sequential analysis.
#' @param TopK 	Top K metabolites (15 by default) to be used in the sequential analysis.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Mdata A large or small metabolic profile matrix. A matrix with metabolic profiles where the number of rows should be equal to the number of metabolites and number of columns should be equal to number of patients.
#' @param Censor A vector of censoring indicator
#' @param Reduce A boolean parameter indicating if the metabolic profile matrix should be reduced, default is TRUE and larger metabolic profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of metabolites to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Plot A boolean parameter indicating if Plot should be shown. Default is FALSE
#' @param DimMethod	 Dimension reduction method which can either be PLS or PCA.
#' @param ...	 Additinal arguments for plotting and only valid  if Plot=TRUE
#' @return A list containing a data frame with estimated HR along with 95\% CI at each TopK value for the sequential analysis.
#'   \item{Result}{The hazard ratio statistics (HR, Lower confidence interval and upper confidence interval) of the lower riskgroup based for each sequential metabolite analysis}
#'   \item{TopKplot}{A graphical representation of the Result containing the hazard ratio statistics}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}},  \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[MetabolicSurv]{MSpecificCoxPh}}, \code{\link[MetabolicSurv]{SurvPcaClass}}, \code{\link[MetabolicSurv]{SurvPcaClass}}
#' @references
#' \insertRef{ye2}{MetabolicSurv}
#'
#' \insertRef{ye1}{MetabolicSurv}
#' @examples
#' \donttest{
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#'Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## USING THE FUNCTION
#' Example1 = SIMet(TopK = 10, Survival=Data$Survival,
#' Mdata=t(Data$Mdata), Censor=Data$Censor, Reduce = TRUE,
#' Select = 50,Prognostic = Data$Prognostic, Plot = TRUE, DimMethod ="PLS")
#'
#' ## FOR THE HR STATISTICS
#' Example1$Result
#'
#' ## FOR THE GRAPHICAL OUTPUT
#' Example1$TopKplot
#' }
#' @export SIMet


SIMet<-function(TopK=15,
                       Survival,
                       Mdata,
                       Censor,
                       Reduce=TRUE,
                       Select=50,
                       Prognostic=NULL,
                       Plot = FALSE,
                       DimMethod=c("PLS","PCA"),...)
{
  Decrease=FALSE
  DimMethod <- match.arg(DimMethod)

  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Mdata)) stop("Argument 'Mdata' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")


  if (Reduce) {
    DataForReduction<-list(x=Mdata,y=Survival, censoring.status=Censor, metabolitenames=rownames(Mdata))
    TentativeList<-names(sort(abs(superpc::superpc.train(DataForReduction, type="survival")$feature.scores),decreasing =TRUE))[1:Select]
    TentativeList

    ReduMdata<-Mdata[TentativeList,]
  } else {
    ReduMdata<-Mdata
  }


  n.metabolite<-nrow(ReduMdata)
  n.patients<-ncol(ReduMdata)

  object<- MSpecificCoxPh(Survival,ReduMdata,Censor, Reduce = FALSE, Select = 150,Prognostic, Quantile = 0.5)

  Names.Kmetabolites<-object@Metnames
  index.Top.Kmetabolites <- order(object@HRRG[,1],decreasing =Decrease)
  index.Top.Kmetabolites<-index.Top.Kmetabolites[1:TopK]

  TopSet<-Names.Kmetabolites[index.Top.Kmetabolites]

  Result<-matrix(NA,TopK,4)

  for (i in 1:length(TopSet[1:TopK]) ){

    mlist<-1:i

    if (DimMethod=="PLS") {
      Temp<- SurvPlsClass(Survival, ReduMdata[intersect(rownames(ReduMdata),TopSet[mlist]),], Censor, Reduce = FALSE, Select = 15,Prognostic = NULL, Plots = FALSE, Quantile = 0.5)
    } else {
      Temp<-  SurvPcaClass(Survival, ReduMdata[intersect(rownames(ReduMdata),TopSet[mlist]),], Censor, Reduce = FALSE, Select = 150,Prognostic = NULL, Plots = FALSE, Quantile = 0.5)
    }

    if (is.null(Prognostic)) Result[i,]<-c(i,(summary(Temp$SurvFit)[[8]][1,])[-2] )
    if (!is.null(Prognostic)) Result[i,]<-c(i,(summary(Temp$SurvFit)[[8]][1,])[-2] )
  }
  HR <- LowerCI <- UpperCI <- NULL
  colnames(Result)<-c("Topk","HR","LowerCI","UpperCI")
  Result <- data.frame(Result)

  if (Plot) {
    TopKplot <- ggplot2::ggplot(data=Result, ggplot2::aes(x=1:nrow(Result)),HR) +
      ggplot2::geom_errorbar(ggplot2::aes(x =1:nrow(Result), ymax = UpperCI, ymin = LowerCI)) +
      ggplot2::ylab(paste("HR Interval Range Based on ",DimMethod,sep="")) +
      ggplot2::xlab("Top K")+ ggplot2::theme_classic() + ggplot2::geom_point(ggplot2::aes(y=Result[,2],x=1:nrow(Result)),colour="red")

    return(list(Result=Result, TopKplot=TopKplot))
  }

  else
    return(Result=Result)

}
