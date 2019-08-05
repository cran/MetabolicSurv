#'Cross Validations for PCA and PLS based methods
#'
#' This function does cross validation  for the analysis performs by \code{\link[MetabolicSurv]{SurvPcaClass}} and \code{\link[MetabolicSurv]{SurvPlsClass}} functions where the dimension reduction methods can either be PCA and PLS.
#'
#' This function does cross validation for the analysis using two reduction method. The reduction method can be PCA or PLS. If it is PCA then the \code{\link[MetabolicSurv]{SurvPcaClass}} is internally used for the cross validation and \code{\link[MetabolicSurv]{SurvPlsClass}} otherwise.
#' @param Fold Number of times in which the dataset is divided. Default is 3 which implies dataset will be divided into three groups and 2/3 of the dataset will be the train datset and 1/3 will be to train the results.
#' @param Survival A vector of survival time with length equals to number of subjects
#' @param Mdata A large or small metabolic profile matrix. A matrix with metabolic profiles where the number of rows should be equal to the number of metabolites and number of columns should be equal to number of patients.
#' @param Censor A vector of censoring indicator
#' @param Reduce A boolean parameter indicating if the metabolic profile matrix should be reduced, default is TRUE and larger metabolic profile matrix is reduced by supervised pca approach and first pca is extracted from the reduced matrix to be used in the classifier.
#' @param Select Number of metabolites (default is 15) to be selected from supervised PCA. This is valid only if th argument Reduce=TRUE
#' @param Prognostic A dataframe containing possible prognostic(s) factor and/or treatment effect to be used in the model.
#' @param Ncv The Number of cross validation loop. Default is 50 but it is recommended to have at least 100.
#' @param DR The dimension reduction method. It can be either "PCA" for Principle components analysis or "PLS" for Partial least squares.
#' @return A object of class \code{\link[MetabolicSurv]{cvpp}} is returned with the following values
#'   \item{Result}{A dataframe containg the estimated Hazard ratio of the test dataset and the training dataset}
#'  \item{Ncv}{The number of cross validation performed}
#'   \item{Method}{The dimesion reduction method used}
#'   \item{CVtrain}{The training dataset indices matrix used for the cross validation}
#'   \item{CVtest}{The test dataset indices matrix used for the cross validation}
#' \item{Select}{The number of metabolite used for the dimesion reduction method used}
#' @author Olajumoke Evangelina Owokotomo, \email{olajumoke.owokotomo@@uhasselt.be}
#' @author Ziv Shkedy
#' @seealso \code{\link[MetabolicSurv]{SurvPlsClass}},
#' \code{\link[MetabolicSurv]{SurvPcaClass}}
#' @references
#' \insertRef{ye1}{MetabolicSurv}
#'
#' \insertRef{ye2}{MetabolicSurv}
#' @examples
#' ## FIRSTLY SIMULATING A METABOLIC SURVIVAL DATA
#' Data = MSData(nPatients = 100, nMet = 150, Prop = 0.5)
#'
#' ## USING THE FUNCTION
#' Result = CVPcaPls(Fold = 4, Survival = Data$Survival,
#' Mdata = t(Data$Mdata), Censor = Data$Censor, Reduce=TRUE,
#' Select=19, Prognostic= Data$Prognostic,Ncv=55,DR ="PLS")
#'
#' ## GET THE CLASS OF THE OBJECT
#' class(Result)     # An "cvpp" Class
#'
#' ##  METHOD THAT CAN BE USED FOR THE RESULT
#' show(Result)
#' summary(Result)
#' plot(Result)

#' @export CVPcaPls

CVPcaPls<-function(Fold = 3,
                   Survival,
                         Mdata,
                         Censor,
                         Reduce=TRUE,
                         Select=15,
                         Prognostic=NULL,
                         Ncv=5,
                         DR ="PCA"){
  options( warn = -1)

  if (Reduce) {
    DataForReduction<-list(x=Mdata,y=Survival, censoring.status=Censor, featurenames=rownames(Mdata))
    TentativeList<-names(sort(abs(superpc.train(DataForReduction, type="survival")$feature.scores),decreasing =TRUE))[1:Select]
    Mdata<-Mdata[TentativeList,]
  } else {
    Mdata<-Mdata
  }

  n.mets<-nrow(Mdata)
  n.patients<-ncol(Mdata)
  sen<-Censor               # Censoring indicator
  dfi <- Survival

  if (is.data.frame(Prognostic)) {   data1<-data.frame(sen,dfi,Prognostic)
  } else {
    data1<-data.frame(sen,dfi)
  }


  HRp.train <- HRn.train <- matrix(0,Ncv,4)  # Training
  HRp.test <- HRn.test<- matrix(0,Ncv,4)     # Testing

  n.train<-(n.patients-floor(n.patients/Fold))
  n.test<-floor(n.patients/Fold)
  cv.train <-matrix(0,Ncv,n.train)
  cv.test  <-matrix(0,Ncv,n.test)


  set.seed(123)
  pIndex <- c(1:n.patients)
  res <-res1<-res2<-res3<- vector("list", Ncv)
  #------------------------  Begin FOR LOOP :1  --------------------------- i=1
  for (i in 1:Ncv){
    message('Cross validation loop ',i)
    p1<-NA
    p2<-NA

    cv.train[i,] <-sort(sample(pIndex,n.train,replace=F) )
    cv.test[i,] <-c(1:n.patients)[-c(intersect(cv.train[i,] ,c(1:n.patients)))]

    if (DR=="PCA") { #-------------------------------------------------------------------------------------------

  Temp1<-IntermediatePCA(Mdata,Prognostic,Survival,Censor,as.vector(cv.train[i,]))
      pc1  <- Temp1$pc1
      cdata <-Temp1$cdata
      m2 <- Temp1$m0

  Temp2<-IntermediatePCA(Mdata,Prognostic,Survival,Censor,cv.test[i,])
      pc1test <- Temp2$pc1
      ctestdata <- Temp2$cdata


      # classification method A1
      TrtandPC1<-summary(m2)[[7]][c("pc1"),1]
      p1 <-  TrtandPC1[1]*pc1
      p2 <-  TrtandPC1[1]*pc1test
    }#-------------------------------------------------------------------------------------------


    if (DR=="PLS") { #-------------------------------------------------------------------------------------------

      Temp3<-IntermediatePLS(Mdata,Prognostic,Survival,Censor,cv.train[i,])
      pls.comp1  <- Temp3$pc1
      cdata <-Temp3$cdata
      m3 <- Temp3$m0

      Temp4<-IntermediatePLS(Mdata,Prognostic,Survival,Censor,cv.test[i,])
      pls.comp1.test <- Temp4$pc1
      ctestdata <- Temp4$cdata

      # classification method A1
      TrtandPLSc1<-summary(m3)[[7]][c("pc1"),1]
      p1 <- TrtandPLSc1[1]*pls.comp1
      p2 <- TrtandPLSc1[1]*pls.comp1.test

    }#-------------------------------------------------------------------------------------------


    #######################
    ## training set ###########

    TrainSet<-EstimateHR(p1,Data.Survival=cdata, Prognostic=Prognostic[cv.train[i,],],Plots = FALSE, Quantile = 0.5 )
  HRp.train[i,] <- summary(TrainSet$SurvResult)[[8]][1,]

    #######################
    ## test set ###########

    mp1 <- median(p1)
    TestSet<-EstimateHR(p2,Data.Survival=ctestdata, Prognostic=Prognostic[cv.test[i,],],Plots = FALSE, Quantile = 0.5 )

    HRp.test[i,] <- summary( TestSet$SurvResult)[[8]][1,]


  }
  #------------------------  END of FOR LOOP :1


  Results<-data.frame(Training=HRp.train[,1],Test=as.numeric(HRp.test[,1]))

  return(new("cvpp",Results=Results,Ncv=Ncv,Method=DR,CVtrain=cv.train,CVtest=cv.test,Select=n.mets))
}
