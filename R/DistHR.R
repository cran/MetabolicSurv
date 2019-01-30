#'Null Distribution of the Estimated HR
#'
#' This function generates the null distribution of the HR by permutation approach. Several ways of permutation setting can be implemented. That is, function can be used to generate null distributions for four different validation schemes, PLS based, PCA based, Majority votes based and Lasso based.
#'
#' This function generates the null distribution of the HR by permutation approach either using a large metabolite matrix or a reduced version by supervised pca approach. Several ways of permutation setting can be implemented. That is, the function can be used to generate null distributions for four different validation schemes which are PLS based, PCA based, Majority votes based and Lasso based. Note this function internally calls function  \code{\link[MetabolicSurv]{SurvPcaClass}}, \code{\link[MetabolicSurv]{SurvPlsClass}}, \code{\link[MetabolicSurv]{Majorityvotes}}, and \code{\link[MetabolicSurv]{Lasoelacox}}.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Censor A vector of censoring indicator
#' @param Mdata A large or small metabolic profile matrix. A matrix with metabolic profiles where the number of rows should be equal to the number of metabolites and number of columns should be equal to number of patients.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Quantile The cut off value for the classifier, default is the median cutoff
#' @param Reduce A boolean parameter indicating if the metabolic profile matrix should be reduced, default is TRUE and larger metabolic profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of metabolites (default is 15) to be selected from supervised PCA. This is valid only if the argument Reduce=TRUE
#' @param nperm Number of permutations to be used and default 100
#' @param case There are seven different ways on how to call this argument:
#' \enumerate{
#' \item{Permute survival only.}
#' \item{Permute survival and rows of data frame of the prognostic factors.}
#' \item{Permute survival, rows of data frame of the prognostic factors, columns of metabolite matrix independently.}
#' \item{Permute metabolite matrix only.}
#' }
#' @param Validation There are four different validation schemes where the null distribution can be estimated. That is c("PLSbased","PCAbased","L1based","MVbased").
#' @return A object of class \code{\link[MetabolicSurv]{perm}} is returned with the following values
#'   \item{HRobs}{Estimated HR for low risk group on the original data}
#'   \item{HRperm}{Estimated HR for low risk group on the permuted data}
#'   \item{nperm}{Number of permutations carried out}
#'   \item{Validation}{The validation scheme that was used}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}}, \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[MetabolicSurv]{SurvPcaClass}}, \code{\link[MetabolicSurv]{SurvPlsClass}}, \code{\link[MetabolicSurv]{Majorityvotes}}, \code{\link[MetabolicSurv]{Lasoelacox}}, \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[MetabolicSurv]{Lasoelacox}}
#' @examples
#' \donttest{
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#' Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## USING THE FUNCTION
#' Example <- DistHR(Survival = Data$Survival,Mdata = t(Data$Mdata),
#' Censor = Data$Censor,Reduce=FALSE,Select=15,Prognostic=Data$Prognostic,
#' Quantile = 0.5, nperm=10, case=2, Validation=c("L1based"))
#' }
#' @export DistHR

DistHR<-function(Survival,
                 Censor,
                 Mdata,
                 Prognostic=NULL,
                 Quantile=0.5,
                 Reduce= FALSE,
                 Select = 15,
                 nperm=100,
                 case=2,
                 Validation=c("PLSbased","PCAbased","L1based","MVbased")

){

  options(warn=-1)
  Validation <- match.arg(Validation)
  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Mdata)) stop("Argument 'Mdata' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")


  if (Reduce) {
    DataForReduction<-list(x=Mdata,y=Survival, censoring.status=Censor, featurenames=rownames(Mdata))
    TentativeList<-names(sort(abs(superpc.train(DataForReduction, type="survival")$feature.scores),decreasing =TRUE))[1:Select]

    Mdata<-Mdata[TentativeList,]
  } else {
    Mdata<-Mdata
  }


  n.mets<-nrow(Mdata)
  n.patients<-ncol(Mdata)


  HRlowPerm<-matrix(NA,nrow=nperm,ncol=3)
  HRlowObs<-as.vector(rep(NA,3))


  #set.seed(123)
  ind.s<-ind.met<-ind.w<-matrix(NA,nrow=nperm,ncol=n.patients)

  for (i in 1:nperm) {
    ind.s[i,]<-1:n.patients
    ind.w[i,]<-1:n.patients
    ind.met[i,]<-1:n.patients
  }

  switch(case,

         {#case 2: permute survival
           for (i in 1:nperm) {
             ind.s[i,]<-sample(c(1:n.patients),replace=F)
           }

         },

         {#case 3: permute survival, prognostic
           for (i in 1:nperm) {
             ind.s[i,]<-sample(c(1:n.patients),replace=F)
           }
           ind.w<-ind.s

         },


         {#case 4 b: permute survival, prognostic and permute metabolites independently
           for (i in 1:nperm) {
             ind.s[i,]<-sample(c(1:n.patients),replace=F)
           }

           ind.w<-ind.s
           for (i in 1:nperm) {
             ind.met[i,]<-sample(c(1:n.patients),replace=F)
           }

         },

         {#case 4 : permute  metabolites only
           for (i in 1:nperm) {
             ind.met[i,]<-sample(c(1:n.patients),replace=F)
           }
         }
  )

  perPrognostic<-NULL

  if (Validation=="PLSbased") {
    for (i in 1:nperm) {
      message('Permutation loop ',i)
      if (!is.null(Prognostic)) perPrognostic<-Prognostic[ind.w[i,],]

      Temp<-SurvPlsClass(Survival=Survival[ind.s[i,]],Mdata= Mdata[,ind.met[i,]],Censor= Censor[ind.s[i,]], Reduce=FALSE,Select=15, Prognostic=perPrognostic, Plots = FALSE,  Quantile = Quantile )

      if (!is.null(Prognostic))  HRlowPerm[i,]<-summary(Temp$SurvFit)[[8]][1,][-2]
      if ( is.null(Prognostic))  HRlowPerm[i,]<-summary(Temp$SurvFit)[[8]][-2]
    }


    TempObs<-SurvPlsClass(Survival, Mdata, Censor, Reduce=FALSE,Select=15, Prognostic=Prognostic, Plots = FALSE,  Quantile = Quantile )
    if (!is.null(Prognostic)) HRlowObs<-summary(TempObs$SurvFit)[[8]][1,]
    if ( is.null(Prognostic)) HRlowObs<-summary(TempObs$Survfit)[[8]][-2]
  }



  if (Validation=="PCAbased") {
    for (i in 1:nperm) {
      message('Permutation loop ',i)
      if (!is.null(Prognostic)) perPrognostic<-Prognostic[ind.w[i,],]

      Temp<-SurvPcaClass(Survival=Survival[ind.s[i,]], Mdata=Mdata[,ind.met[i,]], Censor=Censor[ind.s[i,]], Reduce=FALSE,Select=15, Prognostic=perPrognostic, Plots = FALSE,  Quantile = Quantile )

      if (!is.null(Prognostic))  HRlowPerm[i,]<-summary(Temp$SurvFit)[[8]][1,][-2]
      if ( is.null(Prognostic))  HRlowPerm[i,]<-summary(Temp$SurvFit)[[8]][-2]

    }
    TempObs<-SurvPcaClass(Survival, Mdata, Censor, Reduce=FALSE,Select=15, Prognostic=Prognostic, Plots = FALSE,  Quantile = Quantile )

    if (!is.null(Prognostic)) HRlowObs<-summary(Temp$SurvFit)[[8]][1,][-2]
    if ( is.null(Prognostic)) HRlowObs<-summary(Temp$SurvFit)[[8]][-2]
  }



  if (Validation=="MVbased") {
    for (i in 1:nperm) {
      message('Permutation loop ',i)
      if (!is.null(Prognostic)) perPrognostic<-Prognostic[ind.w[i,],]

      Ana1<-MSpecificCoxPh( Survival=Survival[ind.s[i,]], Mdata=Mdata[,ind.met[i,]], Censor=Censor[ind.s[i,]], Reduce=FALSE,Select=15, Prognostic=perPrognostic, Quantile = Quantile)

      Temp<-Majorityvotes(Ana1,Prognostic=Prognostic[ind.w[i,],], Survival=Survival[ind.s[i,]],Censor=Censor[ind.s[i,]],J=1)

      if (!is.null(Prognostic)) HRlowPerm[i,]<-summary(Temp$Model.result)[[8]][1,][-2]
      if ( is.null(Prognostic)) HRlowPerm[i,]<-summary(Temp$Model.result)[[8]][-2]
    }

    Ana2<-MSpecificCoxPh( Survival,
                          Mdata,
                          Censor,
                          Reduce=FALSE,
                          Select=15,
                          Prognostic=Prognostic,
                          Quantile = Quantile)

    TempObs<-Majorityvotes(Ana2,Prognostic, Survival,Censor,J=1)

    if (!is.null(Prognostic)) HRlowObs<-(summary(TempObs$Model.result)[[8]])[1,][-2]
    if ( is.null(Prognostic)) HRlowObs<-(summary(TempObs$Model.result)[[8]])[-2]
  }



  if (Validation=="L1based") {
    for (i in 1:nperm) {
      message('Permutation loop ',i)
      Temp<-NA
      if (!is.null(Prognostic)) perPrognostic<-Prognostic[ind.w[i,],]

      try(Temp<-Lasoelacox(Survival=Survival[ind.s[i,]],Censor=Censor[ind.s[i,]],Mdata[,ind.met[i,]], Prognostic=perPrognostic, Plots = FALSE,
                           Quantile = Quantile , Metlist=NULL, Standardize = TRUE,  Alpha=1,Fold = 4, nlambda = 100), silent = TRUE)

      if ((!is.na(Temp))[1]) {
        if (!is.null(Prognostic)) HRlowPerm[i,]<-summary(Temp$SurvFit)[[8]][1,][-2]
        if ( is.null(Prognostic)) HRlowPerm[i,]<-summary(Temp$SurvFit)[[8]][-2]
      }
      if ((is.na(Temp))[1])  HRlowPerm[i,]<-NA
    }
    TempObs<-Lasoelacox(Survival,Censor,Mdata, Prognostic=Prognostic, Plots = FALSE, Quantile = Quantile , Metlist=NULL, Standardize = TRUE, Alpha = 1,
                        Fold = 4, nlambda = 100)

    if (!is.null(Prognostic)) HRlowObs<-summary(TempObs$SurvFit)[[8]][1,][-2]
    if ( is.null(Prognostic)) HRlowObs<-summary(TempObs$SurvFit)[[8]][-2]
  }

  return(new("perm",HRobs=HRlowObs,HRperm=HRlowPerm,nperm=nperm,Validation=Validation))

}
