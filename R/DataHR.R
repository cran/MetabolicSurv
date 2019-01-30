#' Survival and Prognostic Data .
#'
#' A dataset containing the riskscore, survival parameters (Overall survival and censoring indicator) and other pronostic factors of 149 subjects.
#'
#' @format A data frame with 149 rows and 5 variables:
#' \describe{
#' \item{Riskscore}{Riskscores of the subjects}
#'   \item{Survival}{Overall survival of the subjects}
#'   \item{Censor}{Censoring indicator for all the patients; 1= Dead and 0 = Alive}
#'   \item{Gender}{The first prognostic factor which is the gender of all the patients; 1=Male and 0 = Female}
#'   \item{Stage}{The second prognostic factor which is the cancer stage of all the patients; 1= Early stage and 0= Advanced stage}
#'   ...
#' }
#' @keywords datasets
#' @usage data(DataHR)
#' @source \url{https://bmccancer.biomedcentral.com/articles/10.1186/s12885-018-4755-1}
#' @examples
#' data(DataHR)
#' summary(DataHR[,1:2])
"DataHR"
