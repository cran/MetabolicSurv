#' The fcv Class.
#'
#' Class of object returned by function \code{\link[MetabolicSurv]{Icvlasoel}}.
#'
#' @name fcv-class
#' @rdname fcv-class
#' @exportClass fcv
#' @param x	 A fcv class object
#' @param y	 missing
#' @param type Plot type. 1 is the distribution of the inner cross validated HR under test data for each outer iterations and estimated HR on the out of bag data are superimposed. 2 Estimated HR Density for low Risk Group .
#' @param  object A fcv class object
#' @param ...	 The usual extra arguments to generic functions — see \code{\link[graphics]{plot}}, \code{\link[graphics]{plot.default}}
#' @slot Runtime A vector of runtime for each iteration measured in seconds.
#' @slot Fold Number of folds used.
#' @slot Ncv Number of outer cross validations used.
#' @slot Nicv Number of inner cross validations used.
#' @slot TopK The Top metabolites used
#' @slot HRInner A 3-way array in which first, second, and third dimensions correspond to Nicv, 1, and Ncv respectively. This contains estimated HR for low risk group on the out of bag data.
#' @slot HRTest A matrix of survival information for the test dataset based on the out of bag data. It has three columns representing the estimated HR, the 95\% lower confidence interval and the 95\% upper confidence interval.
#' @slot Weight A matrix with columns equals number of TopK metabolites and rows Ncv. Note that Weights are estimated as colMeans of coefficients matrix return from the inner cross validations.
#'
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{CVLasoelacox}},
#' \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[glmnet]{glmnet}}, \code{\link[MetabolicSurv]{Lasoelacox}}
#' @examples
#' \donttest{
#' ## GENERATE SOME METABOLIC SURVIVAL DATA WITH PROGNOSTIC FACTORS
#' Data<-MSData(nPatients=100,nMet=150,Prop=0.5)
#'
#' ## USE THE FUNCTION
#' Eg = Icvlasoel(Data$Survival, Data$Censor, Data$Prognostic,
#' t(Data$Mdata), Fold = 3,Ncv = 5, Nicv = 7, Alpha = 1,
#' TopK = colnames(Data$Mdata[,80:100]), Weights = FALSE)
#'
#' ## GET THE CLASS OF THE OBJECT
#' class(Eg)     # An "fcv" Class
#'
#' ##  METHOD THAT CAN BE USED FOR THIS CLASS
#' show(Eg)
#' summary(Eg)
#' plot(Eg, type =1)
#' }
#' @docType class

setClass("fcv",representation(Runtime="vector",Fold="numeric", Ncv="numeric",Nicv="numeric",TopK="vector",
                                           HRInner="array",HRTest="matrix",Weight="matrix"),
        prototype=list(Runtime=c(NA),Fold=3,Ncv=10, Nicv=100, TopK=c(NA),HRInner=array(NA,c(1,1,1)) , HRTest=matrix(1,1,1), Weight=matrix(1,1,1))
)


#' Method show.
#' @name fcv
#' @rdname fcv-class
#' @exportMethod show
#setGeneric("show", function(object) standardGeneric("show"))

#' @rdname fcv-class
#' @aliases show,fcv-method
setMethod("show",signature="fcv"
          , function(object){
            cat("Further Cross Valdiated Results for  \n Lasso and Elastic Net based Predictive Metabolite signature\n")
            cat("Number of Outer CV used: ", object@Ncv, "\n")
            cat("Number of Inner CV used: ", object@Nicv, "\n")
            cat("Number of metabolites subject : ", length(object@TopK), "\n")
          })

#' Method summary.
#' @name fcv-class
#' @rdname fcv-class
#' @exportMethod summary
#setGeneric("summary", function(object,...) standardGeneric("summary"))

#' @rdname fcv-class
#' @aliases summary,fcv-method
setMethod("summary",signature="fcv", function(object){
  cat("Summary of Further Cross Validations\n")
  cat("Estimated Mean Weights for the Classifier\n")
  wmat<-colSums(object@Weight)/nrow(object@Weight)
  names(wmat)<-object@TopK
  print(wmat)
})


#' Method plot.
#' setGeneric("plot",function(x,y,...){standardGeneric("plot")})
#' @name fcv-class
#' @rdname fcv-class
#' @exportMethod plot

#' @rdname fcv-class
#' @aliases plot,fcv,missing-method
#' @aliases fcv-method
setMethod("plot", signature(x="fcv", y="missing"),
          function(x,  y, type=1, ...) {
            if (class(x)!="fcv") stop("Invalid class object")
            if (type==1) {
              dotsCall <- substitute(list(...))
              ll <- eval(dotsCall)
              if(!hasArg("xlab")) ll$xlab <- "Outer Iteration \n The red dot is the median of the HR for the out of bag dataset"
              if(!hasArg("ylab")) ll$ylab <- "HR estimate"
              ll$main <- "Distribution of HR on Test Data \n for low risk group"
              if(!hasArg("cex.lab")) ll$cex.lab <- 0.8
              if(!hasArg("cex.main")) ll$cex.main <- 1
              if(!hasArg("col")) ll$col <- 3
              ll$x<-sapply(1:x@Ncv,function(j) x@HRInner[,1,j])
              do.call(boxplot,args=ll)
              points(1: x@Ncv,x@HRTest[,1],pch=19,col="red",cex = 1.5)
              #for (i in 1: x@Ncv) lines(c(i+0.05,i+0.05),
                                     #c(x@HRTest[i,2],x@HRTest[i,3]),col="red")

            }

            if (type==2) {

              dotsCall <- substitute(list(...))
              ll <- eval(dotsCall)
              if(!hasArg("xlab")) ll$xlab <- "Estimated HR \n The thick vetical line is the median(HR) of the out of bag dataset"
              if(!hasArg("ylab")) ll$ylab <- ""
              ll$main <- "Estimated HR Density for low Risk Group "
              if(!hasArg("cex.lab")) ll$cex.lab <- 0.8
              if(!hasArg("cex.main")) ll$cex.main <- 1
              if(!hasArg("col")) ll$col <- 2

              ll$x <- density(x@HRTest[,1])

              do.call(plot,args=ll)
              #qq<-quantile(sort(x@HRTest[,1]),prob=c(0.05,0.5,0.95))
              #abline(v=qq[1],col=3)
              #abline(v=qq[2],col=3)
              abline(v=median(x@HRTest[,1]),col=3,lwd=5)
            }

          })

