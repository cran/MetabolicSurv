
#'Quantile sensitivity analysis
#'
#' The function performs sensitivity of the cut off quantile for obtaining the risk group obtained under \code{\link[MetabolicSurv]{SurvPlsClass}}, \code{\link[MetabolicSurv]{SurvPcaClass}} or \code{\link[MetabolicSurv]{Lasoelacox}} requires for the survival analysis and classification.
#'
#' This function investigates how each analysis differs from the general median cutoff of 0.5, therefore to see the sensitive nature of the survival result different quantiles ranging from 10th percentile to 90th percentiles were used. The sensitive nature of the quantile is investigated under \code{\link[MetabolicSurv]{SurvPlsClass}}, \code{\link[MetabolicSurv]{SurvPcaClass}} or \code{\link[MetabolicSurv]{Lasoelacox}} while relate to the 3 different Dimension method to select from.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Mdata A large or small metabolic profile matrix. A matrix with metabolic profiles where the number of rows should be equal to the number of metabolites and number of columns should be equal to number of patients.
#' @param Censor A vector of censoring indicator
#' @param Reduce A boolean parameter indicating if the metabolic profile matrix should be reduced, default is TRUE and larger metabolic profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of metabolites (default is 15) to be selected from supervised PCA. This is valid only if th argument Reduce=TRUE
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Plots A boolean parameter indicating if the graphical represenataion of the analysis should be shown. Default is FALSE and it is only valid for the PCA or PLS dimension method.
#' @param DimMethod The dimension method to be used. PCA implies using the \code{\link[MetabolicSurv]{SurvPcaClass}}, PLS uses \code{\link[MetabolicSurv]{SurvPcaClass}} while SM uses the \code{\link[MetabolicSurv]{Lasoelacox}} which ruses the shrinkage method techniques such as lasso and elastic net.
#' @param Alpha The mixing parameter for glmnet (see \code{\link[glmnet]{glmnet}}). The range is 0<= Alpha <= 1. The Default is 1
#' @return A Dataframe is returned depending on weather a data reduction method should be used or not. The dataframe contains the HR of the low risk group for each percentile.

#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}},\code{\link[MetabolicSurv]{EstimateHR}},
#' \code{\link[MetabolicSurv]{SurvPcaClass}},
#'  \code{\link[MetabolicSurv]{SurvPlsClass}},\code{\link[MetabolicSurv]{Lasoelacox}}
#' @references
#' \insertRef{ye1}{MetabolicSurv}
#' @examples
#' \donttest{
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#'Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## USING THE PCA METHOD
#'Result = QuantileAnalysis(Data$Survival,t(Data$Mdata),
#'Data$Censor,Reduce=FALSE, Select=150, Prognostic=Data$Prognostic,
#'Plots = TRUE,DimMethod="PCA",Alpha=1)
#'
#' ## USING THE PLS METHOD
#'Result = QuantileAnalysis(Data$Survival,t(Data$Mdata),
#'Data$Censor,Reduce=FALSE, Select=150, Prognostic=Data$Prognostic,
#'Plots = TRUE,DimMethod="PLS",Alpha=1)
#'
#' ## USING THE SM METHOD
#'Result = QuantileAnalysis(Data$Survival,t(Data$Mdata),
#'Data$Censor,Reduce=FALSE, Select=150, Prognostic=Data$Prognostic,
#'Plots = FALSE,DimMethod="SM",Alpha=1)
#'}

#' @export QuantileAnalysis

QuantileAnalysis<-function(Survival,Mdata,Censor,Reduce=TRUE, Select=150, Prognostic=NULL, Plots = FALSE,DimMethod=c("PLS","PCA","SM"), Alpha =1){

  DimMethod <- match.arg(DimMethod)

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


  n.metabolites<-nrow(ReduMdata)
  n.patients<-ncol(ReduMdata)

  cutpoint <- seq(0.10, 0.9, 0.05)

  Grinding <- matrix(NA,ncol=5,nrow=length(cutpoint))
  Grinding2 <- list()

  for (i in 1:length(cutpoint)) {
    message('Running analysis for Quantile = ',cutpoint[i])
    if (DimMethod=="PLS") {
      Temp<- SurvPlsClass(Survival, Mdata, Censor, Reduce = FALSE, Select = 150, Prognostic, Plots = FALSE, Quantile = cutpoint[i])
    } else if (DimMethod=="PCA") {
      Temp<-  SurvPcaClass(Survival, Mdata, Censor, Reduce = FALSE, Select = 150, Prognostic, Plots = FALSE, Quantile = cutpoint[i])
    } else{
      Temp <- Lasoelacox(Survival, Censor, Mdata, Prognostic, Quantile = cutpoint[i],
                         Metlist = NULL, Plots = FALSE, Standardize = TRUE, Alpha = Alpha)
    }
    #if (substring(DimMethod, 1, 1)=="P"){
      Grinding[i,]<-c(summary(Temp$SurvFit)[[8]][1,],cutpoint[i])
      colnames(Grinding)<-c("EstimatedHR","IHR","LowerCI","UpperCI","Quantile")
    #}
    #else{
      #Grinding2[i,] <- c(summary(Temp2$Result)[[8]][1,],cutpoint[i])
      #colnames(Grinding2)<-c("EstimatedHR","IHR","LowerCI","UpperCI","Quantile")
    #}
  }

  Grinding[is.infinite(Grinding[,3]), 3] <- Grinding[1,1]
  if (Plots) {
    #if (substring(DimMethod, 1, 1)=="P"){
      gplots::plotCI(x=Grinding[,1],xaxt="n", ui=Grinding[,4],li=Grinding[,3],  col="black", barcol="red", lwd=2,ylab="HR",xlab="CutOFF (Quantile)",main=paste("Dimension reduction using ",DimMethod,sep=""))
    axis(side=1,at=1:length(cutpoint),seq(0.10, 0.9, 0.05)*100,cex=0.7)
    #}
  }

  data1<-Grinding[,-2]
  colnames(data1)<-c("EstimatedHR","LowerCI","UpperCI","Quantile")
  #if (substring(DimMethod, 1, 1)=="P"){
    return(data1)
  #}else{
    #return(Result)
  #}
}
