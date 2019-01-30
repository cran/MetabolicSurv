#' The cvsim Class.
#'
#' Class of object returned by function \code{\link[MetabolicSurv]{CVSim}}.
#'
#' @name cvsim-class
#' @rdname cvsim-class
#' @exportClass cvsim
#' @param x	 A cvsim class object
#' @param y	 missing
#' @param type Plot type. 1 distribution of the HR under test For the Top K metabolites using PCA. 2 distribution of the HR under test For the Top K metabolites using PLS.
#' @param  object A cvsim class object
#' @param ...	 The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}
#' @slot HRpca A 3-way array in which first, second, and third dimensions correspond to number of metabolites, Hazard ratio infromation(Estimated HR, LowerCI and UpperCI), and number of cross validation respectively. This contains the estimated HR on test data and dimension reduction method is PCA.
#' @slot HRpls A 3-way array in which first, second, and third dimensions correspond to number of metabolites, Hazard ratio infromation(Estimated HR, LowerCI and UpperCI), and number of cross validation respectively. This contains the estimated HR on test data and dimension reduction method is PLS.
#' @slot Nmets The number of metabolites in the reduced matrix
#' @slot Ncv The number of cross validation done
#' @slot Top A sequence of top k metabolites considered. Default is Top=seq(5,100,by=5)
#'
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{CVPcaPls}}, \code{\link[MetabolicSurv]{SurvPcaClass}}, \code{\link[MetabolicSurv]{SurvPlsClass}}
#' @examples
#' \donttest{
#' ## GENERATE SOME METABOLIC SURVIVAL DATA WITH PROGNOSTIC FACTORS
#' Data<-MSData(nPatients=100,nMet=150,Prop=0.5)
#'
#' ## FIRST IS THE NETABOLITE BY METABOLITE ANALYSIS
#' w = CVMetSpecificCoxPh(Fold=3,Survival=Data$Survival,
#' Mdata=t(Data$Mdata),Censor= Data$Censor,Reduce=TRUE,
#' Select=150,Prognostic=Data$Prognostic,Quantile = 0.5,Ncv=3)
#'
#' ## USING THE FUNCTION
#' Result = CVSim(w, Top = seq(5, 100, by = 5), Survival=Data$Survival,
#'  Censor=Data$Censor, Prognostic = Data$Prognostic)
#'
#' ## GET THE CLASS OF THE OBJECT
#' class(Result)     # A "cvsim" Class
#'
#' ##  METHOD THAT CAN BE USED FOR THE RESULT
#' show(Result)
#' summary(Result)
#' plot(Result, type =2)
#' }
#' @importFrom methods setClass setGeneric setMethod setRefClass

setClass("cvsim",representation(HRpca="array",HRpls="array",Nmets="numeric",Ncv="numeric",Top="numeric"),
         prototype=list(HRpca=array(NA,dim=c(1,1,1)),HRpls=array(NA,dim=c(1,1,1)),Nmets=1,Ncv=3,Top=seq(5,100,by=5))
)

#' Method show.
#' @name cvsim
#' @rdname cvsim-class
#' @exportMethod show
#setGeneric("show", function(object) standardGeneric("show"))

#' @rdname cvsim-class
#' @aliases show,cvsim-method
setMethod("show",signature="cvsim"
          , function(object){
            cat("Cross Validation for sequentially increase metabolites Analysis\n")
            cat("Number of Top K matbolitea used: ", object@Top, "\n")
            cat("Number of cross valdiations used: ", object@Ncv, "\n")
          })

#' Method summary.
#' @name cvsim-class
#' @rdname cvsim-class
#' @exportMethod summary
#setGeneric("summary", function(object,...) standardGeneric("summary"))

#' @rdname cvsim-class
#' @aliases summary,cvsim-method
setMethod("summary",signature="cvsim", function(object){
  cat("Results Based on Test Data\n")
  cat("Summary of Cross Validated Top K metabolites  Analysis\n")
  cat("Estimated Median of the HR for the cross Validated HR for Top K metabolites \n")
  nn<-dim(object@HRpca)[3]
  mean.alpha<- sapply(1:nn,function(i) median(object@HRpca[,1,i],na.rm=T))
  se.alphal<- sapply(1:nn,function(i) quantile(object@HRpca[,1,i],na.rm=T,probs = c(0.025)))
  se.alphau<- sapply(1:nn,function(i) quantile(object@HRpca[,1,i],na.rm=T,probs = c(0.975)))
  mx1<-paste(round(mean.alpha,3),"(",round(se.alphal,3),"-",round(se.alphau,3),")",sep="")


  mean.alpha<- sapply(1:nn,function(i) median(object@HRpls[,1,i],na.rm=T))
  se.alphal<- sapply(1:nn,function(i) quantile(object@HRpls[,1,i],na.rm=T,probs = c(0.025)))
  se.alphau<- sapply(1:nn,function(i) quantile(object@HRpls[,1,i],na.rm=T,probs = c(0.975)))
  mx3<-paste(round(mean.alpha,3),"(",round(se.alphal,3),"-",round(se.alphau,3),")",sep="")


  HR<-data.frame(rbind(mx1,mx3))
  colnames(HR)<-object@Top
  rownames(HR)<-c("(PCA)","(PLS)")
  print(HR)
}
)

#' Method plot.
#' setGeneric("plot",function(x,y,...){standardGeneric("plot")})
#' @name cvsim-class
#' @rdname cvsim-class
#' @exportMethod plot

#' @rdname cvsim-class
#' @aliases plot,cvsim,missing-method
#' @aliases cvsim-method
setMethod("plot", signature(x="cvsim", y="missing"),
          function(x,  y, type=1, ...) {
            if (class(x)!="cvsim") stop("Invalid class object")
            if (type==1) {
            nn<-dim(x@HRpca)[3]
            PC.HRp<-x@HRpca[,1,1:nn]
            colnames(PC.HRp)<-x@Top

            dotsCall <- substitute(list(...))
            ll <- eval(dotsCall)
            if(!hasArg("xlab")) ll$xlab <- "Top K Metabolites"
            if(!hasArg("ylab")) ll$ylab <- "Cross Validated HR"
            ll$main <- "Estimated HR on Test Data \n for Top K Metabolites (PCA)"
            if(!hasArg("cex.lab")) ll$cex.lab <- 1.2
            if(!hasArg("cex.main")) ll$cex.main <- 1.3
            if(!hasArg("col")) ll$col <- 1:nn
            if(!hasArg("ylim")) ll$ylim <- c(0,max(PC.HRp))
            ll$x<-PC.HRp
            do.call(boxplot,args=ll)
            }

            if (type==2) {
            nn<-dim(x@HRpca)[3]
            PL.HRp<-x@HRpls[,1,1:nn]
            colnames(PL.HRp)<-x@Top

            dotsCall <- substitute(list(...))
            lll <- eval(dotsCall)
            if(!hasArg("xlab")) lll$xlab <- "Top K Metabolites"
            if(!hasArg("ylab")) lll$ylab <- "Cross Validated HR"
            lll$main <- "Estimated HR Test Data \n for Top K Metabolites (PLS)"
            if(!hasArg("cex.lab")) lll$cex.lab <- 1.2
            if(!hasArg("cex.main")) lll$cex.main <- 1.3
            if(!hasArg("col")) lll$col <- 1:nn
            if(!hasArg("ylim")) lll$ylim <- c(0,max(PL.HRp))
            lll$x<-PL.HRp
            do.call(boxplot,args=lll)
            return(invisible())
          }
}
)
