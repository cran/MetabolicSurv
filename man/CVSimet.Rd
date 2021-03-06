% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CVSimet.R
\name{CVSimet}
\alias{CVSimet}
\title{Cross validation for sequentially increases metabolites}
\usage{
CVSimet(Object, Top = seq(5, 100, by = 5), Survival, Censor, Prognostic = NULL)
}
\arguments{
\item{Object}{An object of class \code{\link[MetabolicSurv]{cvmm}}}

\item{Top}{The Top k number of metabolites to be used}

\item{Survival}{A vector of survival time with length equals to number of subjects}

\item{Censor}{A vector of censoring indicator}

\item{Prognostic}{A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.}
}
\value{
A object of class \code{\link[MetabolicSurv]{cvsim}} is returned with the following values
\item{HRpca}{A 3-way array in which first, second, and third dimensions correspond to number of metabolites, Hazard ratio infromation(Estimated HR, LowerCI and UpperCI), and number of cross validation respectively. This contains the estimated HR on test data and dimension reduction method is PCA.}
\item{HRpls}{A 3-way array in which first, second, and third dimensions correspond to number of metabolites, Hazard ratio infromation(Estimated HR, LowerCI and UpperCI), and number of cross validation respectively. This contains the estimated HR on test data and dimension reduction method is PLS.}
\item{Nmets}{The number of metabolites in the reduced matrix}
\item{Ncv}{The number of cross validation done}
\item{Top}{A sequence of top k metabolites considered. Default is Top=seq(5,100,by=5)}
}
\description{
This function does cross validation for the metabolite by metabolite analysis while sequentially increasing the number of metabolites as specified.
}
\details{
This function firstly processes the cross validation for the metabolite by metabolite  analysis results, and then sequentially considers top k metabolites. The function recompute first PCA or PLS on train data and estimate risk scores on both test and train data only on the metabolite matrix with top k metabolites. Patients are then classified as having low or high risk based on the test data where the cutoff used is median of the risk score. The process is repeated for each top K metabolite sets.
}
\examples{
\donttest{
## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)

## GETTING THE cvmm OBJECT
Result = CVMetSpecificCoxPh(Fold=3,Survival=Data$Survival,
Mdata=t(Data$Mdata),Censor= Data$Censor,Reduce=TRUE,Select=150,
Prognostic=Data$Prognostic,Quantile = 0.5,Ncv=3)

## USING THE FUNCTION
 Result2 = CVSimet(Result, Top = seq(5, 100, by = 5), Data$Survival,
 Data$Censor,Prognostic = Data$Prognostic)

## GET THE CLASS OF THE OBJECT
class(Result2)     # An "cvsim" Class
}
}
\seealso{
\code{\link[MetabolicSurv]{MSpecificCoxPh}}
}
\author{
Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@uhasselt.be}

Ziv Shkedy
}
