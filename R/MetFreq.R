#'Frequency of Selected Metabolites from the LASSO, Elastic-net Cross-Validation
#'
#' The function selects the frquency of selection  from the shrinkage method (LASSO, Elastic-net) based on cross validation, that is the number of times each metabolite occur during the cross-validation process. In case of large metabolomic matrix then the N argument can be used to select metabolites occurence at a particular frequency.
#'
#' This function outputs the mostly selected metabolites during the LASSO and Elastic-net cross validation. Selected top metabolites are ranked based on frequency of selection and also a particular frequency cqn be selected. In addition, it visualizes the selected top metabolites based on the minimum frequency specified.
#' @param Object An object of class \code{\link[MetabolicSurv]{cvle}} returned from the function \code{\link[MetabolicSurv]{CVLasoelacox}}.
#' @param TopK The number of Top K metabolites (5 by default) to be displayed in the frequency of selection graph.
#' @param N  The metqbolites with the specified frequency should be displayed in the frequency of selection graph.
#' @return A vector of metabolites and their frequency of selection. Also, a graphical representation is displayed.
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{cvmm}}, \code{\link[survival]{coxph}},
#' \code{\link[MetabolicSurv]{EstimateHR}},\code{\link[MetabolicSurv]{CVLasoelacox}}
#' @examples
#' \donttest{
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#' Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## CROSS-VALIDATION FOR LASSO AND ELASTIC-NET
#' Result = CVLasoelacox(Survival = Data$Survival,
#' Censor = Data$Censor, Mdata = t(Data$Mdata),
#' Prognostic = Data$Prognostic, Quantile = 0.5,
#' Metlist = NULL,Standardize = TRUE, Reduce=FALSE, Select=15,
#' Alpha = 1,Fold = 4,Ncv = 10,nlambda = 100)
#'
#' ## CONFIRMING THE CLASS
#' class(Result)
#'
#' ## USING THE FUNCTION
#' MetFreq(Result,TopK = 5, N=5)
#' }
#' @export MetFreq

MetFreq<-function(Object,TopK=20,N=3){

  #Decrease=FALSE
  if (class(Object)!="cvle") stop("Invalid object class.")

  MFreq <- Object@Met.mat
  fr<-colSums(MFreq)
  names(fr)<-rownames(Object@Mdata)

  if(!is.null(TopK)){
    top <- sort(fr,decreasing=TRUE)
    topn <- top[1:TopK]
    barplot(topn,las=2,col=rainbow(length(topn)), ylim = c(0,max(topn) + 1),  ylab="",cex.names=0.6,main=paste( "Top ", TopK, " most selected metabolites ",sep=""),cex.lab=1,cex.main=1.5  )
    return(top)
  } else
  {
    top.fr<-fr[fr==N]
    #topn <- sort(top.fr,decreasing=TRUE)
    #topnn <- topn[1:N]
    barplot(top.fr,las=2,col=rainbow(length(top.fr)), ylim = c(0,max(top.fr)),  ylab="",cex.names=0.6,main=paste( "Metabolites ", "with selection frequency=", N ,sep=""),cex.lab=1,cex.main=1.5  )
    return(fr)
  }
}
