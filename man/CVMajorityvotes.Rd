% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CVMajorityVotes.R
\name{CVMajorityvotes}
\alias{CVMajorityvotes}
\title{Cross validation for majority votes}
\usage{
CVMajorityvotes(
  Survival,
  Censor,
  Prognostic = NULL,
  Mdata,
  Reduce = TRUE,
  Select = 15,
  Fold = 3,
  Ncv = 100
)
}
\arguments{
\item{Survival}{A vector of survival time with length equals to number
of subjects}

\item{Censor}{A vector of censoring indicator}

\item{Prognostic}{A dataframe containing possible prognostic(s) factor
and/or treatment effect to be used in the model.}

\item{Mdata}{A large or small metabolic profile matrix. A matrix with
metabolic profiles where the number of rows should be equal to the number
of metabolites and number of columns should be equal to number of patients.}

\item{Reduce}{A boolean parameter indicating if the metabolic profile matrix
should be reduced, default is TRUE and larger metabolic profile matrix is
reduced by supervised pca approach and first pca is extracted from the
reduced matrix to be used in the classifier.}

\item{Select}{Number of metabolites (default is 15) to be selected from
supervised PCA. This is valid only if th argument Reduce=TRUE}

\item{Fold}{Number of times in which the dataset is divided. Default is
3 which implies dataset will be divided into three groups and 2/3 of the
dataset will be the train datset and 1/3 will be to train the results.}

\item{Ncv}{The Number of cross validation loop. Default is 50 but it is
recommended to have at least 100.}
}
\value{
A object of class \code{\link[MetabolicSurv]{cvmv}} is returned
with the following values
\item{HRTrain}{A matrix of survival information for the training dataset.
It has three columns representing the estimated HR, the 95\% lower
confidence interval and the 95\% upper confidence interval.}
\item{HRTest}{A matrix of survival information for the test dataset.
It has three columns representing the estimated HR, the 95\% lower
confidence interval and the 95\% upper confidence interval.}
\item{Ncv}{The number of cross validation used}
\item{Mdata}{The Metabolite data matrix that was used for the analysis either same as Mdata or a reduced version.}
\item{Progfact}{The names of prognostic factors used}
}
\description{
This function does cross validation for the Majority votes based
classification.
}
\details{
This function does cross validation for the Majority votes based
classification which is a cross validated approach to
\code{\link[MetabolicSurv]{Majorityvotes}}.
}
\examples{
\donttest{
## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)

## USING THE FUNCTION
Result = CVMajorityvotes(Survival=Data$Survival,Censor=Data$Censor,
Prognostic=Data$Prognostic, Mdata=t(Data$Mdata), Reduce=FALSE,
Select=15, Fold=3, Ncv=10)

## GET THE CLASS OF THE OBJECT
class(Result)     # An "cvmv" Class

##  METHOD THAT CAN BE USED FOR THE RESULT
show(Result)
summary(Result)

}
}
\seealso{
\code{\link[MetabolicSurv]{Majorityvotes}}
}
\author{
Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@uhasselt.be}

Ziv Shkedy
}
