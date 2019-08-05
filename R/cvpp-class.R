#' The cvpp Class.
#'
#' Class of object returned by function \code{\link[MetabolicSurv]{CVPcaPls}}.
#'
#' @name cvpp-class
#' @rdname cvpp-class
#' @exportClass cvpp
#' @param x	 A cvpp class object
#' @param y	 missing
#' @param  object A cvpp class object
#' @param ...	 The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}
#' @slot Results A dataframe containg the estimated Hazard ratio of the test dataset and the training dataset
#' @slot Ncv The number of cross validation performed
#' @slot Method The dimesion reduction method used
#' @slot CVtrain The training dataset indices matrix used for the cross validation
#' @slot CVtest The test dataset indices matrix used for the cross validation
#' @slot Select The number of metabolite used for the dimesion reduction method used
#'
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{CVPcaPls}}, \code{\link[MetabolicSurv]{SurvPcaClass}}, \code{\link[MetabolicSurv]{SurvPlsClass}}
#' @examples
#' ## GENERATE SOME METABOLIC SURVIVAL DATA WITH PROGNOSTIC FACTORS
#' Data<-MSData(nPatients=100,nMet=150,Prop=0.5)
#'
#' ## USING THE FUNCTION
#' Result = CVPcaPls(Fold = 4, Survival = Data$Survival,
#' Mdata = t(Data$Mdata), Censor = Data$Censor, Reduce=TRUE,
#' Select=19, Prognostic= Data$Prognostic,Ncv=55,DR ="PLS")
#'
#' ## GET THE CLASS OF THE OBJECT
#' class(Result)     # A "cvpp" Class
#'
#' ##  METHOD THAT CAN BE USED FOR THE RESULT
#' show(Result)
#' summary(Result)
#' plot(Result)

setClass("cvpp",representation(Results="data.frame",Ncv="numeric",Method="vector",CVtrain="matrix",CVtest="matrix",Select="numeric"),
         prototype=list(Results=data.frame(1),Ncv=numeric(),Method="PCA",CVtrain=matrix(0,0,0),CVtest=matrix(0,0,0),Select=numeric()))

#' Method show.
#' @name cvpp
#' @rdname cvpp-class
#' @exportMethod show
#setGeneric("show", function(object) standardGeneric("show"))

#' @rdname cvpp-class
#' @aliases show,cvpp-method
setMethod("show",signature="cvpp"
          , function(object){
            cat("Cross Validation based HR estimation methods:",object@Method,"\n")
            cat("Dimension Reduction Method: ", object@Method,"\n", sep="")
            cat("Number of CVs : ", object@Ncv,"\n", sep="")
            cat("Number of Metabolites: ", object@Select, "\n")
          })




#' Method summary.
#' @name cvpp-class
#' @rdname cvpp-class
#' @exportMethod summary
#setGeneric("summary", function(object,...) standardGeneric("summary"))

#' @rdname cvpp-class
#' @aliases summary,cvpp-method
setMethod("summary",signature="cvpp", function(object){
  cat("Summary of computation Analysis\n")
  cat("Reduction method used :",object@Method,"\n")
  cat("Number of cross validation: ", object@Ncv, "\n")
  cat("Estimated  quantiles of the HR in the train dataset \n")
  print(quantile(object@Results[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
  cat("\n")
  cat("Estimated  quantiles of the HR in the test dataset \n")
  print(quantile(object@Results[,2],probs=c(0.05,0.25,0.5,0.75,0.95)))
})

#" setGeneric("plot", function(x,y, ...) standardGeneric("plot"))

#' Method plot.
#' setGeneric("plot",function(x,y,...){standardGeneric("plot")})
#' @name cvpp-class
#' @rdname cvpp-class
#' @exportMethod plot
#' @import graphics

#' @rdname cvpp-class
#' @aliases plot,cvpp,missing-method
#' @aliases cvpp-method
setMethod("plot", signature(x="cvpp", y="missing"),
          function(x,  y, ...) {
            dotsCall <- substitute(list(...))
            ll <- eval(dotsCall)
            if(!hasArg("xlab")) ll$xlab <- ""
            if(!hasArg("ylab")) ll$ylab <- "HR estimate"
            if(!hasArg("main")) ll$main <- paste("Distribution of HR on Training and Test Set \n for Low risk group using ", x@Method,sep="")
            if(!hasArg("cex.lab")) ll$cex.lab <- 1.5
            if(!hasArg("cex.main")) ll$cex.main <- 1
            if(!hasArg("ylim")) ll$ylim <- c(0,max(x@Results))
            if(!hasArg("col")) ll$col <- c(2,3)
            ll$x<-x@Results
            do.call(boxplot,args=ll)
            #boxplot(x@Results,ylim=c(0,max(x@Results,na.rm=T)),names=c("Train ","Test"),main=mtitle,ylab="HR",col=c("green","red"))
            return(invisible())
          }
)
