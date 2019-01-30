#'Frequency of Selected Metabolites from the Metabolite specific Cross Validation
#'
#' The function selects the top K metabolites from the metabolite by metabolite cross validation, that is the number of times each metabolite occur during the cross-validation process. In case of large metabolomic matrix then the minfreq is required to select metabolites occurence more than a certain frequency. If the number specified is greater than the maximum frequency of seection then max(frequency)-1 will be automatically used as MinFreq.
#'
#' This function outputs the mostly selected metabolites during the metabolite specific cross validation. Selected Top metabolites are ranked based on estimated HR for low risk group. Therefore top metabolites should have minimum HR estimates. Number of top K metabolites needed should be prespecified. In addition, it visualizes the selected top metabolites based on the minimum frequency specified..
#' @param Object An object of class \code{\link[MetabolicSurv]{cvmm}} returned from the function \code{\link[MetabolicSurv]{CVMetSpecificCoxPh}}.
#' @param TopK The number of Top K metabolites (20 by default) to be used in the analysis
#' @param Minfreq  The minimum frequency required for a particular meabolites to be included in the result.
#' @param N The number of Top N metabolites to be displayed in the frequency of selection graph
#' @return A vector of metabolites and their frequency of selection. Also, a graphical representation is displayed.
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{cvmm}}, \code{\link[survival]{coxph}},
#' \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[MetabolicSurv]{MSpecificCoxPh}},\code{\link[MetabolicSurv]{CVMetSpecificCoxPh}}
#' @examples
#' \donttest{
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#' Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## METABOLITE SPECIFIC CROSS VALIDATION
#' Result = CVMetSpecificCoxPh(Fold=3,Survival=Data$Survival,
#' Mdata=t(Data$Mdata),Censor= Data$Censor,Reduce=FALSE,
#' Select=15,Prognostic=Data$Prognostic,Quantile = 0.5,Ncv=30)
#'
#' ## CONFIRMING THE CLASS
#' class(Result)
#'
#' ## USING THE FUNCTION
#' MetFreq(Result,TopK = 5,Minfreq = 20, N=5)
#' }
#' @export MetFreq

MetFreq<-function(Object,TopK=20,Minfreq=5, N =10){

  Decrease=FALSE
  if (class(Object)!="cvmm") stop("Invalid object class.")

  Mdata <- Object@Rdata
  matrix.Mdata<-matrix(0,Object@Ncv,Object@n.mets)
  Mdatas.names <- rownames(Mdata)
  for (i in 1:Object@Ncv){ #i=1

    SetTE<-Object@HRTest[,,i] #test
    SetT<-Object@HRTrain[,,i] # train

    Names.KMdatas<-rownames(Mdata)

    Top.KMdatas.GSplus<-data.frame(rownames(Mdata),SetTE[,-2])
    colnames(Top.KMdatas.GSplus)<-c("Mdata","HR","LCI","UCI")


    Top.KMdatas.GSminus<-data.frame(rownames(Mdata),SetT[,-2])
    colnames(Top.KMdatas.GSminus)<-c("Mdata","HR","LCI","UCI")

    cilevel <- 1-0.05*nrow(Top.KMdatas.GSplus)/Object@n.mets

    HRpadj <- exp(log(SetTE[,1]) + log(SetTE[,c(3,4)])-log(SetTE[,1])*qnorm(cilevel)/1.96)  #
    res.topkMdatas<-data.frame(Top.KMdatas.GSplus,HRpadj,SetT[,1])

    colnames(res.topkMdatas)<-c("Mdata","HRTest","LCI","UCI","FDRLCI","FDRUCI","HRTrain")

    sort.topK<-res.topkMdatas[order(res.topkMdatas[,c("HRTest")],decreasing=Decrease),]

    topg20<-sort.topK[1:TopK,c("Mdata")]
    matrix.Mdata[i,is.element(Mdatas.names,topg20)]<-1
  }

  fr<-colSums(matrix.Mdata)
  names(fr)<-Mdatas.names

  if (max(fr)<Minfreq) Minfreq<-max(fr)-1

  top.fr<-fr[fr>Minfreq]
  top <- sort(top.fr,decreasing=TRUE)
  topn <- top[1:N]
  barplot(topn,las=2,col=rainbow(length(topn)), ylim = c(0,max(top.fr) + 3),  ylab="",cex.names=0.6,main=paste( "Top ", N, " most selected metabolites ",sep=""),cex.lab=1,cex.main=1.5  )

  return(top)

}
