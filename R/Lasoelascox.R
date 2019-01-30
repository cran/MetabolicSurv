#'   Wapper function for glmnet
#'
#' The function uses the glmnet function to firstly do the variable selection either with Lasso, Elastic net or ridge regressions before the survial analysis. The survival analysis is based on the selected metabolites in the presence or absene of prognostic factors.
#'
#' This is a wrapper function for glmnet and it fits models using either Lasso, Elastic net and Ridge regressions. This is done in the presence or absene of prognostic factors. The prognostic factor when avaialable will always be forced to be in the model so no penalty for it. Optimum lambda will be used to select the non-zero shrinkage coefficients, the nonzero selceted metabolites will thus be used in the survival analysis and in calculation of the risk scores.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Censor A vector of censoring indicator
#' @param Mdata A large or small metabolic profile matrix. A matrix with metabolic profiles where the number of rows should be equal to the number of metabolites and number of columns should be equal to number of patients.
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Quantile The cut off value for the classifier, default is the median cutoff
#' @param Metlist A list of metabolites to be considered in the model usually smaller than the metabolites in the Mdata . Default is to use all metabolites available
#' @param Plots A boolean parameter indicating if plots should be shown. Default is FALSE. If TRUE, the first plot is the partial likelihood deviance against the logarithmn of each lambda while the second is the coefficients versus the lamdas
#' @param Standardize A Logical flag for the standardization of the metabolite matrix, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is standardize=TRUE.
#' @param Alpha The mixing parameter for glmnet (see \code{\link[glmnet]{glmnet}}). The range is 0<= Alpha <= 1. The Default is 1
#' @param Fold number of folds to be used for the cross validation. Its value ranges between 3 and the numbe rof subjects in the dataset
#' @param nlambda The number of lambda values - default is 100 as in glmnet.
#' @return A object is returned with the following values
#'   \item{Coefficients.NonZero}{The coefficients of the selected metabolites}
#'   \item{Selected.Mets}{The selected metabolites}
#'   \item{n}{The number of selected metabolites}
#'   \item{Risk.scores}{The risk scores of the subjects}
#'   \item{Risk.group}{The risk classification of the subjects based on the specified quantile}
#'   \item{SurvFit}{The cox analysis of the riskgroup based on the selected metabolites and the prognostic factors}
#'   \item{Select}{A Boolean argument indicating if there was selection or not}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[survival]{coxph}},
#' \code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[glmnet]{glmnet}},
#' @examples
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#' Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## USING THE FUNCTION
#' Results = Lasoelacox(Survival=Data$Survival, Censor=Data$Censor,
#' Mdata=t(Data$Mdata), Prognostic = Data$Prognostic, Quantile = 0.5,
#' Metlist = NULL, Plots = FALSE, Standardize = TRUE, Alpha = 1)
#'
#' ## VIEW THE SELECTED METABOLITES
#' Results$Selected.mets
#' ## NUMBER OF SELECTED METABOLITES
#' Results$n
#'
#' ## VIEW THE CLASSIFICATION GROUP OF EACH SUBJECT
#' Results$Risk.Group
#'
#' ## VIEW THE SURVIVAL ANALYSIS RESULT
#' Results$SurvFit
#'
#' ## TO CHECK IF THERE WAS ANY SELECTION
#' Results$Select
#' @import utils
#' @import stats
# #' @import Biobase


#' @export Lasoelacox


Lasoelacox <- function (Survival,
          Censor,
          Mdata,
          Prognostic,
          Quantile = 0.5,
          Metlist = NULL,
          Plots = FALSE,
          Standardize = TRUE,
          Alpha = 1,
          Fold = 4,
          nlambda = 100)
{

  if (missing(Survival)) stop("Argument 'Survival' is missing...")
  if (missing(Mdata)) stop("Argument 'Mdata' is missing...")
  if (missing(Censor)) stop("Argument 'Censor' is missing...")

  n.mets<-nrow(Mdata)
  n.patients<-ncol(Mdata)

  if(is.null(Metlist)){
    Mdata <- Mdata
  } else {
    Mdata <- Mdata[is.element(rownames(Mdata),Metlist),]
  }

  metnames <- rownames(Mdata)

  if(is.null(Prognostic)){
    Data <- Mdata
    Data.Full <- t(Mdata)
    Penalty <- rep(1,row(Data.Full))
  } else{
    Data <- t(Mdata)
    Data.Full <- cbind(Data,Prognostic)
    Penalty <- c(rep(1,ncol(Data)),rep(0,ncol(Prognostic)))
  }

  # Survival times must be larger than 0

  Survival[Survival <= 0] <- quantile(Survival, probs = 0.01)
  Lasso.Cox.CV <- glmnet::cv.glmnet(x = as.matrix(Data.Full),
                            y = survival::Surv(as.vector(Survival),as.vector(Censor) == 1),
                            family = 'cox',
                            alpha = Alpha,
                            nfolds= Fold,
                            nlambda = nlambda,
                            penalty.factor = Penalty,
                            standardize = Standardize)


  # Results of the cv.glmnet procedure

  Lambda <- Lasso.Cox.CV$lambda.min
  Alllamda <- Lasso.Cox.CV$lambda
  Coefficients <- coef(Lasso.Cox.CV, s =Lambda)
  Coefficients.NonZero <- Coefficients[Coefficients[, 1] != 0, ]

  if (!is.null(dim(Coefficients.NonZero)))
  {
    Selected.mets <- setdiff(colnames(Data.Full),colnames(Prognostic))
    Coefficients.NonZero <- Coefficients[c(Selected.mets,colnames(Prognostic)),]
    Select <- "F"
  }else{
    Select <- "T"
  }

  if (is.null(Prognostic))
  {
    Selected.mets <- names(Coefficients.NonZero)
  }

  else
  {
    Selected.mets <- setdiff(names(Coefficients.NonZero), colnames(Prognostic))
  }



  # What to do if no metabolite is selected?
  # Decrease lambda until a metabolite is selected
  # Going through the list of lambda values and repeat the above procedure

  n <- length(Selected.mets)

  Lambda <- Lasso.Cox.CV$lambda.min
  Lambda.Sequence <- Lasso.Cox.CV$lambda
  Lambda.Index <- which(Lambda.Sequence == Lasso.Cox.CV$lambda.min)

  Lambda.Index.Add = 0

  while (n <= 0) {
    Lambda.Index.Add <- Lambda.Index.Add + 1
    Coefficients <- coef(Lasso.Cox.CV, s = 0.12)
    Coefficients.NonZero <- Coefficients[Coefficients[, 1] != 0,]

    if (!is.null(dim(Coefficients.NonZero)))
    {
      Selected.mets <- setdiff(colnames(Data.Full),colnames(Prognostic))
      Coefficients.NonZero <- Coefficients[c(Selected.mets,colnames(Prognostic)),]
      Select <- "F"
    }else{
      Select <- "T"
    }

    if (is.null(Prognostic))
    {
      Selected.mets <- names(Coefficients.NonZero)
    }
    else
    {
      Selected.mets <- setdiff(names(Coefficients.NonZero), colnames(Prognostic))
    }
    Lambda <- Lambda.Sequence[Lambda.Index + Lambda.Index.Add]
    n <- length(Selected.mets)
  }

  Risk.Scores <- as.vector(Coefficients.NonZero[Selected.mets] %*% t(Data[,Selected.mets]))

  # Estimate the HR by running the appropriate function

  Data.Survival <- cbind(Survival,Censor)
  Estimation.HR <- EstimateHR(Risk.Scores = Risk.Scores,Data.Survival = Data.Survival,
                       Prognostic=Prognostic, Plots = FALSE, Quantile =Quantile)
  Risk.Group <- Estimation.HR$Riskgroup
  Cox.Fit.Risk.Group <- Estimation.HR$SurvResult

  # Produce plots if requested

  if (Plots)
  {
    par(mfrow=c(1,2))
    # Plot 1

    plot(log(Lasso.Cox.CV$lambda), Lasso.Cox.CV$cvm,
         main = expression(paste(alpha,'=', 1 , sep=" ")),
         ylim = c(min(Lasso.Cox.CV$cvlo), max(Lasso.Cox.CV$cvup)),
         xlab = expression(log(lambda)), ylab = 'Partial likelihood deviance',
         pch = 19, col = 'red')

    for (i in 1:length(Lasso.Cox.CV$cvm))
      lines(log(c(Lasso.Cox.CV$lambda[i], Lasso.Cox.CV$lambda[i])), c(Lasso.Cox.CV$cvlo[i], Lasso.Cox.CV$cvup[i]))

    Lasso.Cox <- glmnet::glmnet(x = as.matrix(Data.Full),
                        y = survival::Surv(as.vector(Survival),as.vector(Censor) == 1),
                        family = 'cox',
                        alpha = Alpha,
                        nlambda = 100,
                        penalty.factor = Penalty,
                        standardize = Standardize)

    # Plot 2

    plot(Lasso.Cox, xvar = 'lambda', label = TRUE,xlab=expression(lambda))
    abline(v = log(Lambda), lwd = 2, lty = 2, col = 'red')

  }
  return(list(Coefficients.NonZero=Coefficients.NonZero,
              Selected.mets = Selected.mets,
              n = n,
              Risk.Scores=Risk.Scores,
              Risk.Group=Risk.Group,
              SurvFit=Cox.Fit.Risk.Group,
              Select = Select))
}
