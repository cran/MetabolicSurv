#' The perm Class.
#'
#' Class of object returned by function \code{\link[MetabolicSurv]{DistHR}}.
#'
#' @name perm-class
#' @rdname perm-class
#' @exportClass perm
#' @param x	 A perm class object
#' @param y	 missing
#' @param  object A perm class object
#' @param ...	 The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}
#' @slot HRobs Estimated HR for low risk group on the original data.
#' @slot HRperm Estimated HR for low risk group on the permuted data
#' @slot nperm Number of permutations carried out.
#' @slot Validation The validation scheme that was used.
#'
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{DistHR}}, \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[MetabolicSurv]{SurvPcaClass}}, \code{\link[MetabolicSurv]{SurvPlsClass}}, \code{\link[MetabolicSurv]{Majorityvotes}}, \code{\link[MetabolicSurv]{Lasoelacox}}, \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[MetabolicSurv]{Lasoelacox}}
#' @examples
#' \donttest{
#' ## GENERATE SOME METABOLIC SURVIVAL DATA WITH PROGNOSTIC FACTORS
#' Data<-MSData(nPatients=100,nMet=150,Prop=0.5)
#'
#' ## USING THE FUNCTION
#' Example <- DistHR(Survival = Data$Survival,Mdata = t(Data$Mdata),
#' Censor = Data$Censor,Reduce=FALSE,Select=15,Prognostic=Data$Prognostic,
#' Quantile = 0.5, nperm=10, case=2, Validation=c("L1based"))
#'
#' ## GET THE CLASS OF THE OBJECT
#' class(Example)     # A "perm" Class
#'
#' ##  METHOD THAT CAN BE USED FOR THIS CLASS
#' show(Example)
#' summary(Example)
#' plot(Example)
#' }
#' @note The first, third and last vertical line on the plot are the lower,
#' median  and upper CI of the permuted data estimated HR while
#' the red line is the estimated HR of the original data

setClass("perm",representation(HRobs="vector",HRperm="matrix",nperm="numeric",Validation="vector"),
         prototype=list(HRobs=as.vector(rep(NA,3)),HRperm=matrix(1,1,1),nperm=100,Validation=c(NA))
)

#' Method show.
#' @name perm
#' @rdname perm-class
#' @exportMethod show
#setGeneric("show", function(object) standardGeneric("show"))

#' @rdname perm-class
#' @aliases show,perm-method
setMethod("show",signature="perm"
          , function(object){
            cat("Estimated Null Ditribution of the ",object@Validation,"\n")
            cat("Number of Permutations: ", object@nperm, "\n")
          })


#' Method summary.
#' @name perm-class
#' @rdname perm-class
#' @exportMethod summary
#setGeneric("summary", function(object,...) standardGeneric("summary"))

#' @rdname perm-class
#' @aliases summary,perm-method
setMethod("summary",signature="perm", function(object){
  cat("Summary of Permutation Analysis\n")
  cat("validation scheme used :",object@Validation,"\n")
  cat("Number of Permutations: ", object@nperm, "\n")
  cat("\n")
  cat("Estimated  quantiles of the null distribution of HR\n")
  print(quantile(object@HRperm[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
  cat("\n")
  cat("Estimated HR on original data\n")
  ttt<-object@HRobs
  names(ttt)<-c("Estimate","lower95CI","Upper95CI")
  print(ttt)
  })


#' Method plot.
#' setGeneric("plot",function(x,y,...){standardGeneric("plot")})
#' @name perm-class
#' @rdname perm-class
#' @exportMethod plot

#' @rdname perm-class
#' @aliases plot,perm,ANY-method
#' @aliases perm-method
setMethod("plot", signature("perm"),
          function(x, y, ...) {
            if (class(x)!="perm") stop("Invalid class object")

            HR<-x@HRperm[,1]
            HR<-na.exclude(HR)
            n<-x@nperm
            vv=x@HRobs[1]
            pvalue<-sum(vv>HR)/n

             dotsCall <- substitute(list(...))
            ll <- eval(dotsCall)
            if(!hasArg("xlab")) ll$xlab <- paste("Estimated HR: Emperical p-value = ", pvalue ,sep="")
            if(!hasArg("ylab")) ll$ylab <- ""
            ll$main <- "Null Distribution of HR on Permuted Data \n for low risk group"
            if(!hasArg("cex.lab")) ll$cex.lab <- 0.8
            if(!hasArg("cex.main")) ll$cex.main <- 1
            if(!hasArg("col")) ll$col <- 1
            if(!hasArg("ylim")) ll$ylim <- c(0,4)

            ll$x <- density(HR,from=0,to=(max(HR)+0.25))
            do.call(plot,args=ll)
            abline(v=vv,col=2)
            #CI for permuated cases
            qq<-quantile(sort(HR),prob=c(0.05,0.95))
            abline(v=qq[1],col=3)
            abline(v=qq[2],col=3)
            abline(v=median(HR),col=3,lwd=3)
            })

