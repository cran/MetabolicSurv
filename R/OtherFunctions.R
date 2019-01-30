#### mean for each column of x
average=function(x){
  transposex = t(x)
  nn = ncol(transposex)
  c(transposex %*% rep(1, nn))/nn
}


#-----------------------------------END OF function for the mean ------------------------------------


hrDens<-function(HR,vv=1,mymain=""){
  r <- density(HR[,1],from=0,to=3)
  plot(r,main=mymain,ylim=c(0,3),xlab="Estimated HR")

  abline(v=vv,col=2)
  #CI for permuated cases
  qq<-quantile(sort(HR[,1]),prob=c(0.05,0.5,0.95))
  abline(v=qq[1],col=3)
  abline(v=qq[2],col=3)
  abline(v=median(HR[,1]),col=3,lwd=5)

  pvalue<-sum(vv>HR[,1])/nrow(HR)
  return(list(pvalue,qq))
}


#---------------------------------------------------------------------------------------------------
#----------------------------------- function for PCA ----------------------------------------------
f.pca=function (x)
{
  ca <- match.call()
  if (ncol(x) > nrow(x)) {
    u = princomp(t(x))
    u$call = ca
    return(u)
  }
  xb <- x - (mu <- average(x))
  xb.svd <- svd(xb)
  u <- t(xb) %*% xb.svd$u
  dimnames(u)[[2]] <- paste("PC", 1:ncol(u), sep = "")
  l <- xb.svd$u
  dimnames(l) <- list(paste("V", 1:nrow(l), sep = ""), paste("Comp.",
                                                             1:ncol(l), sep = ""))
  class(l) <- "loadings"
  sd = xb.svd$d/sqrt(ncol(x))
  names(sd) <- paste("Comp.", 1:length(sd), sep = "")
  u <- list(sdev = sd, loadings = l, center = mu, scale = rep(1,
                                                              length(mu)), n.obs = ncol(x), scores = u, call = ca)
  class(u) <- "princomp"
  return(u)
}


#---------------------------------------------------------------------------------------------------
#----------------------------------- Variance ----------------------------------------------
#### variance for each row or column of x
#### 1=row and 2 = columns
variance <- function(x, type =1){
  if (type==1){
    apply(x,1,var)
  }else {
    apply(x,2,var)
  }
}




#---------------------------------------------------------------------------------------------------
#----------------------------------- Intermediate PCA ----------------------------------------------

IntermediatePCA<-function(Mdata,Prognostic,Survival,Censor,index){
  if (is.matrix(Mdata)) {
    pc1 <- f.pca(as.matrix(Mdata[,index]))[[6]][,1]
  } else {
    pc1<-Mdata[,index]
  }

  if (is.null(Prognostic)) {

    cdata <- data.frame(Survival=Survival[index],Censor=Censor[index],pc1)
    m0 <- survival::coxph(Surv(Survival, Censor==1) ~ pc1,data=cdata)
  }

  if (!is.null(Prognostic)) {
    if (is.data.frame(Prognostic)) {
      nPrgFac<-ncol(Prognostic)
      cdata <- data.frame(Survival=Survival[index],Censor=Censor[index],pc1,Prognostic[index,])
      NameProg<-colnames(Prognostic)
      eval(parse(text=paste( "m0 <-survival::coxph(Surv(Survival, Censor==1) ~ pc1",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
    } else {
      stop(" Argument 'Prognostic' is NOT a data frame ")
    }

  }

  return(list(m0=m0,pc1=pc1,cdata=cdata))
}


#---------------------------------------------------------------------------------------------------
#----------------------------------- Intermediate PLS ----------------------------------------------


IntermediatePLS<-function(Mdata,Prognostic,Survival,Censor,index){
  if (is.matrix(Mdata)) {

    PLSforGSK<-data.frame(1:length(index))
    PLSforGSK$g<-as.matrix(t(Mdata[,index]))
    colnames(PLSforGSK)[1]<-c("Survival")
    PLSforGSK[,1]<-Survival[index]
    plsr.1 <- pls::plsr(Survival ~ g, method="simpls",ncomp = 2, scale =TRUE,data = PLSforGSK, validation =  "CV")
    pc1<-pls::scores(plsr.1)[,1] # extract the first com
  } else {
    pc1<-Mdata[,index]
  }

  if (is.null(Prognostic)) {
    cdata <- data.frame(Survival=Survival[index],Censor=Censor[index],pc1)
    m0 <- survival::coxph(Surv(Survival, Censor==1) ~ pc1,data=cdata)
  }
  if (!is.null(Prognostic)) {
    if (is.data.frame(Prognostic)) {
      nPrgFac<-ncol(Prognostic)
      cdata <- data.frame(Survival=Survival[index],Censor=Censor[index],pc1,Prognostic[index,])
      NameProg<-colnames(Prognostic)
      eval(parse(text=paste( "m0 <-survival::coxph(Surv(Survival, Censor==1) ~ pc1",paste("+",NameProg[1:nPrgFac],sep="",collapse =""),",data=cdata)" ,sep="")))
    } else {
      stop(" Argument 'Prognostic' is NOT a data frame ")
    }
  }
  return(list(m0=m0,pc1=pc1,cdata=cdata))
}
