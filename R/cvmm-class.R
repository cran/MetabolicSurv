#' The cvmm Class.
#'
#' Class of object returned by function \code{\link[MetabolicSurv]{CVMetSpecificCoxPh}}.
#'
#' plot signature(x = "cvmm"): Plots for \code{\link[MetabolicSurv]{CVMetSpecificCoxPh}} class analysis results.
#'
#' Any parameters of \code{\link[graphics]{plot.default}} may be passed on to this particular plot method.
#'
#' @rdname cvmm-class
#' @exportClass cvmm
#' @param x	 A CVMetSpecificCoxPh class object
#' @param y	 missing
#' @param object A CVMetSpecificCoxPh class object
#' @param which This specify which metabolite for which estimated HR information need to be visualized. By default results of the first metabolite is used.
#' @param ...	The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}
#' @slot HRTrain A 3-way array, The first dimension is the number of metabolites, the second dimension is the HR statistics for the low risk group in the train dataset (HR,1/HR LCI, UCI) while the third dimension is the number of cross validation performed.
#' @slot HRTest A 3-way array, The first dimension is the number of metabolites, the second dimension is the HR statistics for the low risk group in the test dataset (HR,1/HR LCI, UCI) while the third dimension is the number of cross validation performed.
#' @slot train The selected subjects for each CV in the train dataset
#' @slot test The selected subjects for each CV in the test dataset
#' @slot n.mets The number of metabolite used in the analysis
#' @slot Ncv The number of cross validation performed
#' @slot Rdata The Metabolite data matrix that was used for the analysis either same as Mdata or a reduced version
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{CVMetSpecificCoxPh}}
#' @examples
#' \donttest{
#' ## GENERATE SOME METABOLIC SURVIVAL DATA WITH PROGNOSTIC FACTORS
#' Data<-MSData(nPatients=100,nMet=150,Prop=0.5)
#'
#' ## USING THE FUNCTION
#' Result = CVMetSpecificCoxPh(Fold=3,Survival=Data$Survival,
#' Mdata=t(Data$Mdata),Censor= Data$Censor,Reduce=TRUE,
#' Select=150,Prognostic=Data$Prognostic,Quantile = 0.5,Ncv=3)
#'
#' ## GET THE CLASS OF THE OBJECT
#' class(Result)     # An "cvmm" Class
#'
#' ##  METHOD THAT CAN BE USED FOR THIS CLASS
#' show(Result)
#' summary(Result)
#' plot(Result)
#' }

setClass("cvmm",slots = representation(HRTrain="array",HRTest="array",train="matrix",test="matrix",n.mets="numeric",Ncv="numeric",Rdata="matrix"),
         prototype=list(HRTrain=array(NA,dim=c(1,1,1)),HRTest=array(NA,dim=c(1,1,1)),train=matrix(0,0,0), test=matrix(0,0,0),n.mets=1,Ncv=3,Rdata=matrix(0,0,0)))

#' Method show.
#' @name cvmm
#' @rdname cvmm-class
#' @exportMethod show
#setGeneric("show", function(object,...) standardGeneric("show"))

#' @rdname cvmm-class
#' @aliases show,cvmm-method
setMethod("show",signature="cvmm"
          , function(object){
            cat("Cross Validation for the Metabolite specific analysis\n")
            cat("Number of Metabolite used: ", object@n.mets, "\n")
            cat("Number of cross validations  performed: ", object@Ncv, "\n")
          })

#' Method summary.
#' @name cvmm-class
#' @rdname cvmm-class
#' @exportMethod summary
#' @aliases summary,cvmm-method
setMethod("summary",signature="cvmm"
          ,function(object,which=1){
            cat("Summary of Cross Validation for the Metabolite specific analysis\n")
            cat("Estimated Median of the cross Validated HR for Metabolite: ",which,"\n")
            HRTest  <- object@HRTest[which,,][1,]
            HRTrain <- object@HRTrain[which,,][1,]

            mean.alpha <- median(HRTrain,na.rm=T)
            se.alphal <- quantile(HRTrain,na.rm=T,probs = c(0.025))
            se.alphau <- quantile(HRTrain,na.rm=T,probs = c(0.975))
            cat("Estimated HR for Train Dataset \n")
            cat(paste(round(mean.alpha,4),"(",round(se.alphal,4)," , ",
                      round(se.alphau,4),")",sep=""))

            cat("\n")
            mean.alpha <- median(HRTest,na.rm=T)
            se.alphal <- quantile(HRTest,na.rm=T,probs = c(0.025))
            se.alphau <- quantile(HRTest,na.rm=T,probs = c(0.975))
            cat("Estimated HR for Test Dataset \n")
            cat(paste(round(mean.alpha,4),"(",round(se.alphal,4)," , ",
                      round(se.alphau,4),")",sep=""))
            cat("\n")
          }
)



#' Method plot.
#' @name cvmm-class
#' @rdname cvmm-class
#' @exportMethod plot

#' @rdname cvmm-class
#' setGeneric("plot",function(x,y,...){standardGeneric("plot")})
#' @aliases plot,cvmm,ANY-method
#' @aliases cvmm-method
setMethod(f="plot", signature = "cvmm",
          definition =  function(x,  y, which=1, ...) {
            HRTest  <- x@HRTest[which,,][1,]
            HRTrain <- x@HRTrain[which,,][1,]
            Results<-data.frame(HRTrain,HRTest)
            dotsCall <- substitute(list(...))
            ll <- eval(dotsCall)
            if(!hasArg("xlab")) ll$xlab <- ""

            if(!hasArg("ylab")) ll$ylab <- "HR estimate"
            if(!hasArg("main")) ll$main <- paste("Estimated HR of Low risk group for Metabolite ", which, "\n Number of CVs = ",x@Ncv,sep="")
            if(!hasArg("cex.lab")) ll$cex.lab <- 1.5
            if(!hasArg("cex.main")) ll$cex.main <- 1
            if(!hasArg("ylim")) ll$ylim <- c(0,15)
            if(!hasArg("col")) ll$col <- 2:3
            if(!hasArg("names"))  ll$names=c("Training","Test")
            ll$x<-Results
            do.call(boxplot,args=ll)
            return(invisible())
          }

)



