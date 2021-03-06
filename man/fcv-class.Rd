% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MetabolicSurvClass.R
\docType{class}
\name{fcv-class}
\alias{fcv-class}
\title{The fcv Class.}
\description{
The fcv Class.
}
\section{Slots}{

\describe{
\item{\code{Runtime}}{A vector of runtime for each iteration measured in seconds.}

\item{\code{Fold}}{Number of folds used.}

\item{\code{Ncv}}{Number of outer cross validations used.}

\item{\code{Nicv}}{Number of inner cross validations used.}

\item{\code{TopK}}{The Top metabolites used}

\item{\code{HRInner}}{A 3-way array in which first, second, and third dimensions
correspond to Nicv, 1, and Ncv respectively. This contains estimated HR
for low risk group on the out of bag data.}

\item{\code{HRTest}}{A matrix of survival information for the test dataset based
on the out of bag data. It has three columns representing the estimated HR,
the 95\% lower confidence interval and the 95\% upper confidence interval.}

\item{\code{Weight}}{A matrix with columns equals number of TopK metabolites and
rows Ncv. Note that Weights are estimated as colMeans of coefficients
matrix return from the inner cross validations.}
}}

\examples{
\donttest{
## GENERATE SOME METABOLIC SURVIVAL DATA WITH PROGNOSTIC FACTORS
Data<-MSData(nPatients=100,nMet=150,Prop=0.5)

## USE THE FUNCTION
Eg = Icvlasoel(Data$Survival, Data$Censor, Data$Prognostic,
t(Data$Mdata), Fold = 3,Ncv = 5, Nicv = 7, Alpha = 1,
TopK = colnames(Data$Mdata[,80:100]), Weights = FALSE)

## GET THE CLASS OF THE OBJECT
class(Eg)     # An "fcv" Class

##  METHOD THAT CAN BE USED FOR THIS CLASS
show(Eg)
summary(Eg)
plot(Eg, type =1)
}
}
\seealso{
\code{\link[MetabolicSurv]{CVLasoelacox}},
\code{\link[MetabolicSurv]{EstimateHR}}, \code{\link[glmnet]{glmnet}}, \code{\link[MetabolicSurv]{Lasoelacox}}
}
\author{
Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@uhasselt.be}

Ziv Shkedy
}
