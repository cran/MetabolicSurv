#' The ms class
#' @slot Result A list of dataframes of each output object of coxph for the metabolites.
#' @slot HRRG A dataframe with estimated metabolite-specific HR for low risk group and 95 percent CI.
#' @slot Group A matrix of the classification group a subject belongs to for each of the metabolite analysis. The metabolites are on the rows and the subjects are the columns
#' @slot Metnames The names of the metabolites for the analysis
#'
#' @author Olajumoke Evangelina Owokotomo,
#' \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{MSpecificCoxPh}}
#' @examples
#' ## GENERATE SOME METABOLIC SURVIVAL DATA WITH PROGNOSTIC FACTORS
#' Data<-MSData(nPatients=100,nMet=150,Prop=0.5)
#'
#' ## DO THE METABOLITE BY METABOLITE ANALYSIS
#' Eg = MSpecificCoxPh(Survival=Data$Survival, Mdata=t(Data$Mdata),
#' Censor=Data$Censor, Reduce = FALSE, Select = 15,
#' Prognostic=Data$Prognostic, Quantile = 0.5)
#'
#' ## GET THE CLASS OF THE OBJECT
#' class(Eg)     # An "ms" Class
#'
#' ##  METHOD THAT CAN BE USED FOR THIS CLASS
#' show(Eg)
#' summary(Eg)
#' plot(Eg)
#' @import survival methods
#' @export


ms <- setClass("ms",
               slots = list(
                 Result="list",
                 HRRG="matrix",Group="matrix",Metnames="vector"),
               prototype=list(Result=list(1),
                              HRRG=matrix(0,0,0),
                              Group=matrix(0,0,0),
                              Metnames = vector()
                              )
               )
#' Constructor for the ms class
#' @param  Result A list of dataframes of each output object of coxph for the metabolites.
#' @param HRRG A dataframe with estimated metabolite-specific HR for low risk group and 95 percent CI.
#' @param Group A matrix of the classification group a subject belongs to for each of the metabolite analysis. The metabolites are on the rows and the subjects are the columns
#' @param Metnames The names of the metabolites for the analysis
#'
#' @return  object of class ms

ms <- function(Result, HRRG, Group, Metnames){
  new("ms",
      Result = Result,
      HRRG = HRRG,
      Group = Group,
      Metnames = Metnames)
}

#' Show method for ms class
#' @param object An object of class ms
#' @return Metabolite by Metabolite CoxPh Model and the number of
#' metabolites used.
#' @export
setMethod("show",signature="ms"
          , function(object){
            cat("Metabolite by Metabolite CoxPh Model\n")
            cat("Number of Metabolite used: ", length(object@Metnames), "\n")
          })

#' Summary method for ms class
#' @param object An object of class ms
#' @return Metabolite by Metabolite CoxPh summary information
#' @export
setMethod("summary",signature="ms",function(object){
  cat("Summary of Metabolite by Metabolite CoxPh Models\n")
  cat("Number of Metabolites used: ", length(object@Metnames), "\n")
  cat("Top 15 Metabolites out of ", length(object@Metnames), "\n")
  cat("Estimated HR for the low risk group\n")
  # top Metabolites based on upper CI HR GS+
  Names.KMetabolites<-object@Metnames
  index.Top.KMetabolites <- order(object@HRRG[,1],decreasing =FALSE)
  index.Top.KMetabolites<-index.Top.KMetabolites[1:15]
  Top.KMetabolites.GSplus<-data.frame(Metnames=Names.KMetabolites[index.Top.KMetabolites],object@HRRG[index.Top.KMetabolites,])

  colnames(Top.KMetabolites.GSplus)<-c("Metnames","HR GS+","LowerCI","UpperCI")
  # FDR corrected CI for top k  metabolites
  cilevel <- 1-0.05*15/nrow(object@HRRG)

  HRpadj <- exp(log(object@HRRG[index.Top.KMetabolites,1]) + log(object@HRRG[index.Top.KMetabolites,c(2,3)])-log(object@HRRG[index.Top.KMetabolites,1])*qnorm(cilevel)/1.96)  #
  res.topkMetabolites<-data.frame(Top.KMetabolites.GSplus,HRpadj)

  colnames(res.topkMetabolites)<-c("Metnames","HR","LCI","UCI","FDRLCI","FDRUCI")
  print(res.topkMetabolites)
})


#' Plot method for ms class
#' @param x	 A ms class object
#' @param y	 missing
#' @param ...	 The usual extra arguments to generic functions — see  \code{\link[graphics]{plot.default}}
#' @return Metabolite by Metabolite CoxPh plots
#' @export
setMethod(f="plot", signature = "ms",
          definition = function(x,y,...){
            object <-  x
            Names.KMetabolites<-object@Metnames
            index.Top.KMetabolites <- order(object@HRRG[,1],decreasing =FALSE)
            Top.KMetabolites.GSplus<-data.frame(Metnames=Names.KMetabolites,object@HRRG)

            colnames(Top.KMetabolites.GSplus)<-c("Metnames","HR GS+","LowerCI","UpperCI")
            # FDR corrected CI for top k  metabolites
            cilevel <- 1-0.05*15/nrow(object@HRRG)

            HRpadj <- exp(log(object@HRRG[index.Top.KMetabolites,1]) + log(object@HRRG[index.Top.KMetabolites,c(2,3)])-log(object@HRRG[index.Top.KMetabolites,1])*qnorm(cilevel)/1.96)  #
            res.topkMetabolites<-data.frame(Top.KMetabolites.GSplus,HRpadj)

            colnames(res.topkMetabolites)<-c("Metnames","HR","LCI","UCI","FDRLCI","FDRUCI")

            x= 1:length(object@Metnames)
            par(mfrow=c(2,2))
            boxplot(object@HRRG[,1],col="blue",ylab = "Estimated HR",main = "Distribution of all the HR for the metabolite specific analysis")

            Group2 = object@Group
            prop = sapply(1:ncol(Group2),function(k) sum(Group2[,k]=="Low risk")/nrow(Group2))
            barplot(prop, col=grDevices::rainbow(length(Group2[1,])), xlab = "Patient index", ylab = "Proportion",main = "Proportion of being classified as low risk group",ylim = c(0,max(prop)))

            plot(x=x, y = res.topkMetabolites[,2],
                 ylim= c(0,max(res.topkMetabolites[,4])),
                 pch=19, xlab="Metabolites", ylab="Hazard ratio",
                 main="Hazard ratio plot with unadjusted confidence interval")
            arrows(x, res.topkMetabolites[,3], x, res.topkMetabolites[,4], length=0.05, angle=90, code=3)
            abline(h=1,col="red2",lwd=2.0)

            plot(x=x, y = res.topkMetabolites[,2],
                 ylim= c(0,max(res.topkMetabolites[,6])),
                 pch=19, xlab="Metabolites", ylab="Hazard ratio",
                 main="Hazard ratio plot with adjusted confidence interval")
            arrows(x, res.topkMetabolites[,5], x, res.topkMetabolites[,6], length=0.05, angle=90, code=3)
            abline(h=1,col="red2",lwd=2.0)

            return(invisible())
          })


######################### END OF MS CLASS ###############################

#' The cvle Class.
#' @slot Coef.mat A matrix of coefficients with rows equals to number of
#' cross validations and columns equals to number of metabolites.
#' @slot Runtime A vector of runtime for each iteration measured in seconds.
#' @slot lambda A vector of estimated optimum lambda for each iterations.
#' @slot n A vector of the number of selected metabolites
#' @slot Met.mat A matrix with 0 and 1. Number of rows equals to number of
#' iterations and number of columns equals to number of metabolites.
#' 1 indicates that the particular metabolite was selected or had
#' nonzero coefficient and otherwise it is zero.
#' @slot HRTrain A matrix of survival information for the training
#' dataset. It has three columns representing the estimated HR,
#' the 95% lower confidence interval and the 95% upper confidence interval.
#' @slot HRTest A matrix of survival information for the test dataset.
#' It has three columns representing the estimated HR, the 95% lower
#' confidence interval and the 95% upper confidence interval.
#' @slot pld A vector of partial likelihood deviance at each
#' cross validations.
#' @slot Mdata The metabolite matrix that was used for the analysis
#'  which can either be the full the full data or a reduced
#'  supervised PCA version.
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
#' @export


setClass("cvle",representation(Coef.mat="matrix",Runtime="vector",
      lambda="vector",n="vector",Met.mat="matrix",HRTrain="matrix",
      HRTest="matrix",pld="vector",Mdata="matrix"),
      prototype=list(Coef.mat=matrix(1,1,1),Runtime=c(NA),
      lambda=c(NA), n=c(NA), Met.mat=matrix(1,1,1),HRTrain=matrix(1,1,1),
      HRTest=matrix(1,1,1) , pld=c(NA),Mdata=matrix(1,1,1))
)

#' Constructor for the cvle class
#' @param  Coef.mat A matrix of coefficients with rows equals to number of cross
#' validations and columns equals to number of metabolites.
#' @param  Runtime A vector of runtime for each iteration measured in seconds.
#' @param  lambda A vector of estimated optimum lambda for each iterations.
#' @param  n A vector of the number of selected metabolites
#' @param  Met.mat A matrix with 0 and 1. Number of rows equals to number of
#' iterations and number of columns equals to number of metabolites. 1 indicates
#'  that the particular metabolite was selected or had nonzero coefficient and
#'  otherwise it is zero.
#' @param  HRTrain A matrix of survival information for the training dataset.
#' It has three columns representing the estimated HR, the 95% lower confidence
#'  interval and the 95% upper confidence interval.
#' @param  HRTest A matrix of survival information for the test dataset. It has
#' three columns representing the estimated HR, the 95% lower confidence
#' interval and the 95% upper confidence interval.
#' @param  pld A vector of partial likelihood deviance at each cross
#' validations.
#' @param  Mdata The metabolite matrix that was used for the analysis which can
#' either be the full the full data or a reduced supervised PCA version.
#'
#' @return  object of class cvle

cvle <- function(Coef.mat, Runtime, lambda, n,Met.mat,
                 HRTrain,HRTest, pld, Mdata){
  new("cvle",
      Coef.mat = Coef.mat,
      Runtime = Runtime,
      lambda = lambda,
      n = n,
      Met.mat = Met.mat,
      HRTrain =  HRTrain,
      HRTest = HRTest,
      pld = pld,
      Mdata = Mdata)
}

#' Show method for cvle class
#' @param object An object of class cvle
#' @return Cross Valdiated Results for Lasso and Elastic Net based Predictive Metabolite signature.
#' @export
setMethod("show",signature="cvle"
          , function(object){
            cat("Cross Valdiated Results for Lasso and Elastic Net based Predictive Metabolite signature\n")
            cat("Number of Metabolite used: ", length(rownames(object@Mdata)), "\n")
            cat("Number of CV: ", length(object@lambda), "\n")
          })



#' Summary method for cvle class
#' @param object An object of class cvle
#' @return Cross Valdiated Results for Lasso and Elastic Net based Predictive Metabolite signature summary information
#' @export
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



#' Plot method for cvle class
#' @param x	 A cvle class object
#' @param y	 missing
#' @param type Plot type. 1 distribution of the HR under training and test set. 2 HR vs number selected metabolites.
#' @param ...	 The usual extra arguments to generic functions — see  \code{\link[graphics]{plot.default}}
#' @return Cross Valdiated Results for Lasso and Elastic Net based Predictive Metabolite plots
#' @export
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
          }
)

######################### END OF CVLE CLASS ###############################

#' The cvmm Class.
#' @slot HRTrain A 3-way array, The first dimension is the number of metabolites, the second dimension is the HR statistics for the low risk group in the train dataset (HR,1/HR LCI, UCI) while the third dimension is the number of cross validation performed.
#' @slot HRTest A 3-way array, The first dimension is the number of metabolites,
#' the second dimension is the HR statistics for the low risk group in the test
#' dataset (HR,1/HR LCI, UCI) while the third dimension is the number of
#' cross validation performed.
#' @slot train The selected subjects for each CV in the train dataset
#' @slot test The selected subjects for each CV in the test dataset
#' @slot n.mets The number of metabolite used in the analysis
#' @slot Ncv The number of cross validation performed
#' @slot Rdata The Metabolite data matrix that was used for the analysis either
#' same as Mdata or a reduced version
#'
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
#' @export


setClass("cvmm",slots = representation(HRTrain="array",HRTest="array",
        train="matrix",test="matrix",n.mets="numeric",Ncv="numeric",
        Rdata="matrix"),
       prototype=list(HRTrain=array(NA,dim=c(1,1,1)),
        HRTest=array(NA,dim=c(1,1,1)),train=matrix(0,0,0), test=matrix(0,0,0),
        n.mets=1,Ncv=3,Rdata=matrix(0,0,0)))


#' Constructor for the cvmm class
#' @param HRTrain A 3-way array, The first dimension is the number of metabolites, the second dimension is the HR statistics for the low risk group in the train dataset (HR,1/HR LCI, UCI) while the third dimension is the number of cross validation performed.
#' @param HRTest A 3-way array, The first dimension is the number of metabolites,
#' the second dimension is the HR statistics for the low risk group in the test
#' dataset (HR,1/HR LCI, UCI) while the third dimension is the number of
#' cross validation performed.
#' @param train The selected subjects for each CV in the train dataset
#' @param test The selected subjects for each CV in the test dataset
#' @param n.mets The number of metabolite used in the analysis
#' @param Ncv The number of cross validation performed
#' @param Rdata The Metabolite data matrix that was used for the analysis either
#' same as Mdata or a reduced version
#'
#' @return  object of class cvmm

cvmm <- function(HRTrain, HRTest, train, test,n.mets,Ncv,Rdata)
                 {
  new("cvmm",
      HRTrain = HRTrain,
      HRTest = HRTest,
      train = train,
      test = test,
      n.mets =n.mets,
      Ncv =  Ncv,
      Rdata = Rdata)
}

#' Show method for cvmm class
#' @param object An object of class cvmm
#' @return Cross Valdiated Metabolic Specific CoxPh information
#' @export
setMethod("show",signature="cvmm"
          , function(object){
            cat("Cross Validation for the Metabolite specific analysis\n")
            cat("Number of Metabolite used: ", object@n.mets, "\n")
            cat("Number of cross validations  performed: ", object@Ncv, "\n")
          })




#' Summary method for cvmm class
#' @param object An object of class cvmm
#' @param which This specify which metabolite for which estimated HR
#' @return Cross Valdiated Metabolic Specific CoxPh summary
#' @export
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


#' Plot method for cvmm class
#' @param x	 A cvmm class object
#' @param y	 missing
#' @param which This specify which metabolite for which estimated HR
#' information need to be visualized. By default results of the first metabolite
#' is used.
#' @param ...	 The usual extra arguments to generic functions — see  \code{\link[graphics]{plot.default}}
#' @return Cross Valdiated Metabolic Specific CoxPh plots
#' @export
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

######################### END OF CVMM CLASS ###############################

#' The cvmv Class.
#' @slot HRTrain A matrix of survival information for the training dataset.
#' It has three columns representing the estimated HR, the 95% lower
#' confidence interval and the 95% upper confidence interval.
#' @slot HRTest A matrix of survival information for the test dataset.
#' It has three columns representing the estimated HR, the 95% lower
#' confidence interval and the 95% upper confidence interval.
#' @slot Ncv The number of cross validation used
#' @slot Mdata The Metabolite data matrix that was used for the analysis either
#' same as Mdata or a reduced version.
#' @slot Progfact The names of prognostic factors used
#'
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{Majorityvotes}}, \code{\link[MetabolicSurv]{CVPcaPls}}, \code{\link[MetabolicSurv]{SurvPcaClass}}, \code{\link[MetabolicSurv]{SurvPlsClass}}
#' @examples
#' \donttest{
#' ## GENERATE SOME METABOLIC SURVIVAL DATA WITH PROGNOSTIC FACTORS
#' Data<-MSData(nPatients=100,nMet=150,Prop=0.5)
#'
#' ## USING THE FUNCTION
#' Result = CVMajorityvotes(Survival=Data$Survival,Censor=Data$Censor,
#' Prognostic=Data$Prognostic, Mdata=t(Data$Mdata), Reduce=FALSE,
#' Select=15, Fold=3, Ncv=10)
#'
#' ## GET THE CLASS OF THE OBJECT
#' class(Result)     # A "cvmv" Class
#'
#' ##  METHOD THAT CAN BE USED FOR THE RESULT
#' show(Result)
#' summary(Result)
#'
#' }
#' @export


setClass("cvmv",representation(HRTrain="matrix",HRTest="matrix",Ncv="numeric",
                               Mdata="matrix",Progfact="vector"),
         prototype=list(HRTrain=matrix(1,1,1),HRTest=matrix(1,1,1),Ncv=100,
                        Mdata=matrix(1,1,1),Progfact=c(NA))
)

#' Constructor for the cvmv class
#' @param HRTrain A matrix of survival information for the training dataset.
#' It has three columns representing the estimated HR, the 95% lower
#' confidence interval and the 95% upper confidence interval.
#' @param HRTest A matrix of survival information for the test dataset.
#' It has three columns representing the estimated HR, the 95%
#' lower confidence interval and the 95% upper confidence interval.
#' @param Ncv The number of cross validation used
#' @param Mdata The Metabolite data matrix that was used for the analysis
#' either same as Mdata or a reduced version.
#' @param Progfact The names of prognostic factors used
#'
#' @return  object of class cvmv

cvmv <- function(HRTrain, HRTest, Ncv,Mdata, Progfact)
{
  new("cvmv",
      HRTrain = HRTrain,
      HRTest = HRTest,
      Ncv =  Ncv,
      Mdata = Mdata,
      Progfact = Progfact)
}

#' Show method for cvmv class
#' @param object An object of class cvmv
#' @return Cross validation for Majority Votes Based Classification Analysis information
#' @export
setMethod("show",signature="cvmv"
          , function(object){
            cat("Cross validation for Majority Votes Based Classification Analysis\n")
            cat("Number of cross valdiations used: ", object@Ncv, "\n")
            if(!is.null(object@Progfact)) cat("Prognostic factors used: ",object@Progfact,"\n")
          })




#' Summary method for cvmv class
#' @param object An object of class cvmv
#' @return Cross validation for Majority Votes Based Classification Analysis summary
#' @export
setMethod("summary",signature="cvmv", function(object){
  cat("Summary of majority votes cross validation analysis\n")
  cat("Number of prognostic factor used :",length(object@Progfact),"\n")
  cat("Number of cross validation: ", object@Ncv, "\n")
  cat("Estimated  quantiles of the HR in the train dataset \n")
  print(quantile(object@HRTrain[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
  cat("\n")
  cat("Estimated  quantiles of the HR in the test dataset \n")
  print(quantile(object@HRTest[,1],probs=c(0.05,0.25,0.5,0.75,0.95)))
})




#' Plot method for cvmv class
#' @param x	 A cvmv class object
#' @param y	 missing
#' @param ...	 The usual extra arguments to generic functions — see  \code{\link[graphics]{plot.default}}
#' @return Cross validation for Majority Votes Based Classification Analysis plots
#' @export
setMethod("plot", signature(x="cvmv"),
          function(x,  y, ...) {
            if (class(x)!="cvmv") stop("Invalid class object")
            HRTest<-x@HRTest
            HRTrain<-x@HRTrain
            nCV<-x@Ncv
            dotsCall <- substitute(list(...))
            ll <- eval(dotsCall)
            if(!hasArg("xlab")) ll$xlab <- "Cross Validation index"
            if(!hasArg("ylab")) ll$ylab <- "HR estimate"
            ll$main <- "Estimated HR on Test Set \n for Low risk group"
            if(!hasArg("cex.lab")) ll$cex.lab <- 1.5
            if(!hasArg("cex.main")) ll$cex.main <- 1
            if(!hasArg("col")) ll$col <- 2

            ll$x<-HRTest[,1]
            if(!hasArg("ylim")) ll$ylim <- c(0,max(x@HRTrain,x@HRTest))


            par(mfrow=c(1,2))
            t1 <- which(HRTest[,1]<1)
            do.call(plot,args=ll)
            #plot(HRp.test[,1],ylim=c(0,2),ylab="HR",main="")
            for(i in 1:nCV){
              lines(c(i,i),HRTest[i,2:3])
            }
            for(i in t1){
              lines(c(i,i),HRTest[i,2:3],col=2)
            }
            abline(h=1)


            Results<-data.frame(HRTrain=HRTrain[,1],HRTest=as.numeric(HRTest[,1]))
            ll$x<-Results
            ll$names<-c("Training ","Test ")
            ll$main <- "Estimated HR on Training and Test Set \n for Low risk group"
            if(!hasArg("xlab")) ll$xlab <- ""
            if(!hasArg("ylab")) ll$ylab <- "HR estimate"
            if(!hasArg("cex.lab")) ll$cex.lab <- 1.5
            if(!hasArg("cex.main")) ll$cex.main <- 1
            if(!hasArg("col")) ll$col <- 2:3
            do.call(boxplot,args=ll)
          })


######################### END OF CVMV CLASS ###############################


#' The cvpp Class.
#' @slot Results A dataframe containg the estimated Hazard ratio of the test
#' dataset and the training dataset
#' @slot Ncv The number of cross validation performed
#' @slot Method The dimesion reduction method used
#' @slot CVtrain The training dataset indices matrix used for the cross
#' validation
#' @slot CVtest The test dataset indices matrix used for the cross validation
#' @slot Nmet The number of metabolite used for the dimesion reduction method
#' used
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
#' @export


setClass("cvpp",representation(Results="data.frame",Ncv="numeric",
                               Method="vector",CVtrain="matrix",CVtest="matrix",
                               Nmet="numeric"),
         prototype=list(Results=data.frame(1),Ncv=numeric(),Method="PCA",
                        CVtrain=matrix(0,0,0),CVtest=matrix(0,0,0),
                        Nmet=numeric()))


#' Constructor for the cvpp class
#' @param Results A dataframe containg the estimated Hazard ratio of the
#' test dataset and the training dataset
#' @param Ncv The number of cross validation performed
#' @param Method The dimesion reduction method used
#' @param CVtrain The training dataset indices matrix used for the cross
#' validation
#' @param CVtest The test dataset indices matrix used for the cross
#' validation
#' @param Nmet The number of metabolite used for the dimesion reduction
#' method used
#'
#' @return  object of class cvpp

cvpp <- function(Results, Ncv, Method, CVtrain, CVtest, Nmet)
{
  new("cvpp",
      Results = Results,
      Ncv =  Ncv,
      Method = Method,
      CVtrain = CVtrain,
      CVtest = CVtest,
      Nmet = Nmet)
}

#' Show method for cvpp class
#' @param object An object of class cvpp
#' @return CCross Validations for PCA and PLS based information
#' @export
setMethod("show",signature="cvpp"
          , function(object){
            cat("Cross Validations for PCA and PLS based:",object@Method,"\n")
            cat("Dimension Reduction Method: ", object@Method,"\n", sep="")
            cat("Number of CVs : ", object@Ncv,"\n", sep="")
            cat("Number of Metabolites: ", object@Nmet, "\n")
          })




#' Summary method for cvpp class
#' @param object An object of class cvpp
#' @return Cross Validations for PCA and PLS based summary
#' @export
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

#' Plot method for cvpp class
#' @param x	 A cvpp class object
#' @param y	 missing
#' @param ...	 The usual extra arguments to generic functions — see  \code{\link[graphics]{plot.default}}
#' @return Cross Validations for PCA and PLS based plots
#' @export
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



######################### END OF CVPP CLASS ###############################


#' The cvsim Class.
#' @slot HRpca A 3-way array in which first, second, and third
#' dimensions correspond to number of metabolites, Hazard ratio
#' infromation(Estimated HR, LowerCI and UpperCI),
#' and number of cross validation respectively.
#' This contains the estimated HR on test data and
#' dimension reduction method is PCA.
#' @slot HRpls A 3-way array in which first, second, and third
#' dimensions correspond to number of metabolites,
#' Hazard ratio infromation(Estimated HR, LowerCI and UpperCI),
#' and number of cross validation respectively.
#' This contains the estimated HR on test data and dimension
#' reduction method is PLS.
#' @slot Nmets The number of metabolites in the reduced matrix
#' @slot Ncv The number of cross validation done
#' @slot Top A sequence of top k metabolites considered.
#' Default is Top=seq(5,100,by=5)
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
#' Result = CVSimet(w, Top = seq(5, 100, by = 5), Survival=Data$Survival,
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
#' @export


setClass("cvsim",representation(HRpca="array",HRpls="array",
            Nmets="numeric",Ncv="numeric",Top="numeric"),
         prototype=list(HRpca=array(NA,dim=c(1,1,1)),
      HRpls=array(NA,dim=c(1,1,1)),Nmets=1,Ncv=3,
      Top=seq(5,100,by=5))
)


#' Constructor for the cvsim class
#' @param HRpca A 3-way array in which first, second, and third
#' dimensions correspond to number of metabolites, Hazard ratio
#' infromation(Estimated HR, LowerCI and UpperCI),
#' and number of cross validation respectively.
#' This contains the estimated HR on test data and
#' dimension reduction method is PCA.
#' @param HRpls A 3-way array in which first, second, and
#' third dimensions correspond to number of metabolites,
#' Hazard ratio infromation(Estimated HR, LowerCI and UpperCI),
#' and number of cross validation respectively. This contains the
#'  estimated HR on test data and dimension reduction method is
#'   PLS.
#' @param Nmets The number of metabolites in the reduced matrix
#' @param Ncv The number of cross validation done
#' @param Top A sequence of top k metabolites considered.
#' Default is Top=seq(5,100,by=5)

#'
#' @return  object of class cvsim

cvsim <- function(HRpca, HRpls, Nmets, Ncv, Top)
{
  new("cvsim",
      HRpca = HRpca,
      HRpls = HRpls,
      Nmets = Nmets,
      Ncv = Ncv,
      Top = Top)
}

#' Show method for cvsim class
#' @param object An object of class cvsim
#' @return Cross validation for sequentially increases metabolites information
#' @export
setMethod("show",signature="cvsim"
          , function(object){
            cat("Cross Validation for sequentially increase metabolites Analysis\n")
            cat("Number of Top K matbolitea used: ", object@Top, "\n")
            cat("Number of cross valdiations used: ", object@Ncv, "\n")
          })





#' Summary method for cvsim class
#' @param object An object of class cvsim
#' @return Cross validation for sequentially increases metabolites summary
#' @export
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



#' Plot method for cvsim class
#' @param x	 A cvsim class object
#' @param y	 missing
#' @param type Plot type. 1 distribution of the HR under test For the Top K
#' metabolites using PCA. 2 distribution of the HR under test For the Top K
#' metabolites using PLS.
#' @param ...	 The usual extra arguments to generic functions — see
#'  \code{\link[graphics]{plot.default}}
#' @return Cross validation for sequentially increases metabolites plots
#' @export
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


######################### END OF CVSIM CLASS ###############################

#' The fcv Class.
#' @slot Runtime A vector of runtime for each iteration measured in seconds.
#' @slot Fold Number of folds used.
#' @slot Ncv Number of outer cross validations used.
#' @slot Nicv Number of inner cross validations used.
#' @slot TopK The Top metabolites used
#' @slot HRInner A 3-way array in which first, second, and third dimensions
#' correspond to Nicv, 1, and Ncv respectively. This contains estimated HR
#' for low risk group on the out of bag data.
#' @slot HRTest A matrix of survival information for the test dataset based
#' on the out of bag data. It has three columns representing the estimated HR,
#' the 95% lower confidence interval and the 95% upper confidence interval.
#' @slot Weight A matrix with columns equals number of TopK metabolites and
#' rows Ncv. Note that Weights are estimated as colMeans of coefficients
#' matrix return from the inner cross validations.
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
#' @export


setClass("fcv",representation(Runtime="vector",Fold="numeric", Ncv="numeric",
                Nicv="numeric",TopK="vector",
            HRInner="array",HRTest="matrix",Weight="matrix"),
         prototype=list(Runtime=c(NA),Fold=3,Ncv=10, Nicv=100, TopK=c(NA),
      HRInner=array(NA,c(1,1,1)) , HRTest=matrix(1,1,1), Weight=matrix(1,1,1))
)



#' Constructor for the fcv class
#' @param  Runtime A vector of runtime for each iteration measured in seconds.
#' @param Fold Number of folds used.
#' @param Ncv Number of outer cross validations used.
#' @param Nicv Number of inner cross validations used.
#' @param TopK The Top metabolites used
#' @param HRInner A 3-way array in which first, second, and third dimensions
#' correspond to Nicv, 1, and Ncv respectively. This contains estimated HR
#' for low risk group on the out of bag data.
#' @param HRTest A matrix of survival information for the test dataset based on
#' the out of bag data. It has three columns representing the estimated HR,
#' the 95% lower confidence interval and the 95% upper confidence interval.
#' @param Weight A matrix with columns equals number of TopK metabolites and
#' rows Ncv. Note that Weights are estimated as colMeans of
#' coefficients matrix return from the inner cross validations.
#'
#' @return  object of class fcv

fcv <- function(Runtime, Fold, Ncv, Nicv, TopK, HRInner, HRTest, Weight)
{
  new("fcv",
      Runtime = Runtime,
      Fold = Fold,
      Ncv = Ncv,
      Nicv = Nicv,
      TopK = TopK,
      HRInner = HRInner,
      HRTest = HRTest,
      Weight = Weight)
}

#' Show method for fcv class
#' @param object An object of class fcv
#' @return Inner and Outer Cross Validations for Lasso Elastic Net Survival predictive models and Classification information
#' @export
setMethod("show",signature="fcv"
          , function(object){
            cat("Further Cross Valdiated Results for  \n Lasso and Elastic Net based Predictive Metabolite signature\n")
            cat("Number of Outer CV used: ", object@Ncv, "\n")
            cat("Number of Inner CV used: ", object@Nicv, "\n")
            cat("Number of metabolites subject : ", length(object@TopK), "\n")
          })





#' Summary method for fcv class
#' @param object An object of class fcv
#' @return Inner and Outer Cross Validations for Lasso Elastic Net Survival predictive models and Classification summary
#' @export
setMethod("summary",signature="fcv", function(object){
  cat("Summary of Further Cross Validations\n")
  cat("Estimated Mean Weights for the Classifier\n")
  wmat<-colSums(object@Weight)/nrow(object@Weight)
  names(wmat)<-object@TopK
  print(wmat)
})


#' Plot method for fcv class
#' @param x	 A fcv class object
#' @param y	 missing
#' @param type Plot type. 1 is the distribution of the inner cross validated
#' HR under test data for each outer iterations and estimated HR on the out
#' of bag data are superimposed. 2 Estimated HR Density for low Risk Group.
#' @param ...	 The usual extra arguments to generic functions — see
#'  \code{\link[graphics]{plot.default}}
#' @return Inner and Outer Cross Validations for Lasso Elastic Net Survival predictive models and Classification plots
#' @export
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


######################### END OF FCV CLASS ###############################


#' The perm Class.
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
#' @export


setClass("perm",representation(HRobs="vector",HRperm="matrix",
            nperm="numeric",Validation="vector"),
         prototype=list(HRobs=as.vector(rep(NA,3)),HRperm=matrix(1,1,1),
                        nperm=100,Validation=c(NA))
)


#' Constructor for the perm class
#' @param HRobs Estimated HR for low risk group on the original data.
#' @param HRperm Estimated HR for low risk group on the permuted data
#' @param nperm Number of permutations carried out.
#' @param Validation The validation scheme that was used.
#'
#' @return  object of class perm

perm <- function(HRobs, HRperm, nperm,Validation)
{
  new("perm",
      HRobs = HRobs,
      HRperm = HRperm,
      nperm = nperm,
      Validation = Validation)
}

#' Show method for perm class
#' @param object An object of class perm
#' @return Null Distribution of the Estimated HR information
#' @export
setMethod("show",signature="perm"
          , function(object){
            cat("Estimated Null Ditribution of the ",object@Validation,"\n")
            cat("Number of Permutations: ", object@nperm, "\n")
          })



#' Summary method for perm class
#' @param object An object of class perm
#' @return Null Distribution of the Estimated HR summary
#' @export
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


#' Plot method for perm class
#' @param x	 A perm class object
#' @param y	 missing
#' @param ...	 The usual extra arguments to generic functions — see
#'  \code{\link[graphics]{plot.default}}
#' @return Null Distribution of the Estimated HR plots
#' @export
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


######################### END OF PERM CLASS ###############################
