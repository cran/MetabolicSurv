#' The cvle Class.
#'
#' Class of object returned by function \code{\link[MetabolicSurv]{CVLasoelacox}}.
#'
#' @name cvle-class
#' @rdname cvle-class
#' @exportClass cvle
#' @param x	 A cvle class object
#' @param y	 missing
#' @param type Plot type. 1 distribution of the HR under training and test set. 2 HR vs number selected metabolites.
#' @param  object A cvle class object
#' @param ...	 The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}
#' @slot Coef.mat A matrix of coefficients with rows equals to number of cross validations and columns equals to number of metabolites.
#' @slot Runtime A vector of runtime for each iteration measured in seconds.
#' @slot lambda A vector of estimated optimum lambda for each iterations.
#' @slot n A vector of the number of selected metabolites
#' @slot Met.mat A matrix with 0 and 1. Number of rows equals to number of iterations and number of columns equals to number of metabolites. 1 indicates that the particular metabolite was selected or had nonzero coefficient and otherwise it is zero.
#' @slot HRTrain A matrix of survival information for the training dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.
#' @slot HRTest A matrix of survival information for the test dataset. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.
#' @slot pld A vector of partial likelihood deviance at each cross validations.
#' @slot Mdata The metabolite matrix that was used for the analysis which can either be the full the full data or a reduced supervised PCA version.
#'
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[glmnet]{glmnet}}, \code{\link[MetabolicSurv]{Lasoelacox}}
#' @examples
#' \donttest{
#' ## GENERATE SOME METABOLIC SURVIVAL DATA WITH PROGNOSTIC FACTORS
#' Data<-MSData(nPatients=100,nMet=150,Prop=0.5)
#'
#' ## USE THE FUNCTION
#' Eg = CVLasoelacox(Survival = Data$Survival,Censor = Data$Censor,
#' Mdata = t(Data$Mdata),Prognostic = Data$Prognostic, Quantile = 0.5,
#' Metlist = NULL,Standardize = TRUE, Reduce=FALSE, Select=15,
#' Alpha = 1,Fold = 4,Ncv = 10,nlambda = 100)
#'
#' ## GET THE CLASS OF THE OBJECT
#' class(Eg)     # An "cvle" Class
#'
#' ##  METHOD THAT CAN BE USED FOR THIS CLASS
#' show(Eg)
#' summary(Eg)
#' plot(Eg, type =3)
#' }
#' @docType class

setClass("cvle",representation(Coef.mat="matrix",Runtime="vector",lambda="vector",n="vector",Met.mat="matrix",HRTrain="matrix",HRTest="matrix",pld="vector",Mdata="matrix"), prototype=list(Coef.mat=matrix(1,1,1),Runtime=c(NA),lambda=c(NA), n=c(NA), Met.mat=matrix(1,1,1),HRTrain=matrix(1,1,1) , HRTest=matrix(1,1,1) , pld=c(NA),Mdata=matrix(1,1,1))
)


#' Method show.
#' @name cvle
#' @rdname cvle-class
#' @exportMethod show
#setGeneric("show", function(object) standardGeneric("show"))

#' @rdname cvle-class
#' @aliases show,cvle-method
setMethod("show",signature="cvle"
          , function(object){
            cat("Cross Valdiated Results for Lasso and Elastic Net based Predictive Metabolite signature\n")
            cat("Number of Metabolite used: ", length(rownames(object@Mdata)), "\n")
            cat("Number of CV: ", length(object@lambda), "\n")
          })


#' Method summary.
#' @name cvle-class
#' @rdname cvle-class
#' @exportMethod summary
#setGeneric("summary", function(object,...) standardGeneric("summary"))

#' @rdname cvle-class
#' @aliases summary,cvle-method
setMethod("summary",signature="cvle"
          , function(object){
            cat("Summary of cross valdiateion for Lasso and Elastic Net based Predictive Metabolite signature\n")
            cat("Number of Metabolite used: ", length(rownames(object@Mdata)), "\n")
            cat("Estimated  quantiles of HR on test data\n")
            print(quantile(object@HRTest[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
            cat("\n")
            cat("Estimated quantiles of HR on train data\n")
            print(quantile(object@HRTrain[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
            cat("Mostly selected 30 metabolites:\n")
            Freq=colSums(object@Met.mat)
            names(Freq)<-rownames(object@Mdata)
            sFreq<-sort(Freq,decreasing = TRUE)
            sFreq<-sFreq[sFreq>0]
            maxG<-length(sFreq)
            if (maxG>30) maxG<-30
print(names(sFreq)[1:maxG])
          })


#' Method plot.
#' setGeneric("plot",function(x,y,...){standardGeneric("plot")})
#' @name cvle-class
#' @rdname cvle-class
#' @exportMethod plot

#' @rdname cvle-class
#' @aliases plot,cvle,missing-method
#' @aliases cvle-method
          setMethod(f ="plot", signature(x="cvle", y="missing"),
                    function(x,  y, type=1, ...) {
                      if (class(x)!="cvle") stop("Invalid class object")
                      if (type==1) {

                        DistHR<-data.frame(HRTrain=x@HRTrain[,1],HRTest=x@HRTest[,1])

                        colnames(DistHR)<-c("Training","Test")
                        dotsCall <- substitute(list(...))
                        ll <- eval(dotsCall)
                        if(!hasArg("xlab")) ll$xlab <- ""
                        if(!hasArg("ylab")) ll$ylab <- "HR estimate"
                        ll$main <- "Distribution of HR on Training and Test Set \n for Low risk group"
                        if(!hasArg("cex.lab")) ll$cex.lab <- 1.5
                        if(!hasArg("cex.main")) ll$cex.main <- 1
                        if(!hasArg("col")) ll$col <- 2:3
                        ll$x<-DistHR
                        do.call(boxplot,args=ll)

                      }


                      if (type==2) {

                        HRTest=x@HRTest[,1]
                        dotsCall <- substitute(list(...))
                        ll <- eval(dotsCall)
                        if(!hasArg("xlab")) ll$xlab <- "Estimated HR on Test Data"
                        if(!hasArg("ylab")) ll$ylab <- "Number of non zero coef."
                        ll$main <- "HR vs number of metabolites"
                        if(!hasArg("cex.lab")) ll$cex.lab <- 0.8
                        if(!hasArg("cex.main")) ll$cex.main <- 1
                        if(!hasArg("col")) ll$col <- 2
                        ll$x<-HRTest
                        ll$y<-x@n
                        do.call(plot,args=ll)

                      }

                      # if (type==3) {
                      #
                      #   Freq=colSums(x@Met.mat)
                      #   names(Freq)<-rownames(x@Mdata)
                      #   sFreq<-sort(Freq,decreasing = TRUE)
                      #   sFreq<-sFreq[sFreq>0]
                      #
                      #   maxG<-length(sFreq)
                      #   if (maxG>30) maxG<-30
                      #
                      #
                      #   dotsCall <- substitute(list(...))
                      #   lll <- eval(dotsCall)
                      #
                      #   lll$height<-sFreq[1:maxG]
                      #   if(!hasArg("xlab" ))   lll$xlab<-""
                      #   if(!hasArg("ylab" ))   lll$ylab<-"Frequency"
                      #   if(!hasArg("main"))    lll$main<-"Mostly Selected Metabolites"
                      #   lll$col<-rainbow(maxG)
                      #   lll$names<-names(sFreq)[1:maxG]
                      #   if(!hasArg("cex.lab")) lll$cex.lab <- 1
                      #   if(!hasArg("las"))   lll$las<-2
                      #   lll$cex.names<-0.65
                      #   do.call(barplot,args=lll)
                      #
                      # }

                    }
          )
