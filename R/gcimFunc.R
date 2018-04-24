gcimFunc <- function(mxmp,galaxyy1,res11,chr_name,legend_size,mainline_size,backline_size,margin_space,axis_space,logPCoff,color1,color2,lodthred)
{
  chr_pos <- mxmp[,1:2]
  chr_num <- length(chr_name)
  chr <- matrix(0,chr_num,1)
  pos <- matrix(0,chr_num,1)
  for(i in 1:chr_num)
  {
    temp <- numeric()
    temp <- length(which(chr_pos[,1]==i))
    if(i==1)
    {
      pos[i] <- temp
      chr[i] <- chr_pos[pos[i],2]
    }else{
      pos[i] <- pos[i-1] + temp
      chr[i] <- chr_pos[pos[i],2]
    }
  }
  
  pos_acc <- matrix(0,chr_num,1)
  for(i in 1:chr_num)
  {
    if(i==1){
      pos_acc[i] <- chr[i]
    }else{
      pos_acc[i] <- pos_acc[i-1] + chr[i]
    }
  }
  
  firFil <- res11[,1:2]
  newposadd <- as.matrix(firFil[,2])
  for(i in 1:chr_num)
  {
    temp1 <- numeric()
    temp1 <- which(firFil[,1]==i)
    if(i>1)
    {
      newposadd[temp1] <- newposadd[temp1]+pos_acc[i-1]
    }
  }
  if(is.null(galaxyy1)==FALSE){
    if(is.null(dim(galaxyy1))==TRUE){
      galaxyy1<-matrix(galaxyy1,1,3)
    }
    newres_pos <- galaxyy1[,2]
    res_sumpos <- pos_acc[galaxyy1[which(galaxyy1[,1]>1),1]-1] + galaxyy1[which(galaxyy1[,1]>1),2]
    newres_pos[which(galaxyy1[,1]>1)] <- res_sumpos
    pospic<-c(newres_pos)
    lodpic<-c(galaxyy1[,3])
    resdf <- data.frame(pospic,lodpic)
  }
  
  resp<-as.matrix(res11[,3])
  pmin<-min(resp[resp!=0])
  locsub<-which(resp==0)
  if(length(locsub)!=0){
    subvalue<-10^(1.1*log10(pmin))
    res11[locsub,3]<-subvalue
  }else{
    res11<-res11
  }
  
  negloP <- -log10(as.matrix(res11[,3]))
  if(is.null(galaxyy1)==FALSE){
    par(mar=c(2*margin_space,2*margin_space,2*margin_space,2*margin_space)+margin_space,mgp=c(3*axis_space,axis_space,0))
    plot(pospic,lodpic,type="h",col=color1,xlab="",ylab="Logarithm of odds (LOD)",cex.axis=legend_size,cex.lab=legend_size,lwd=mainline_size,xlim=c(0,max(newposadd)),ylim=c(0,max(lodpic)))
    abline(h=lodthred)
    par(new=TRUE)
    plot(newposadd,negloP,type="l",col=color2,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,xlim=c(0,max(newposadd)),ylim=c(0,logPCoff*max(negloP)))
    axis(side=4,cex.axis=legend_size)
    mtext(expression('-log'[10]*'(P)'),side=4,line=3*axis_space,cex=legend_size)
    abline(v=pos_acc,lty=2,col="gray")
  }else{
    plot(newposadd,negloP,type="l",col=color1,xlab="Genome position (cM)",ylab=expression('Expected -log'[10]*'(P)'),cex.axis=legend_size,cex.lab=legend_size,lwd=mainline_size,xlim=c(0,max(newposadd)),ylim=c(0,logPCoff*max(negloP)))
  }
}

gcimFuncF2 <- function(mxmp,galaxyy1,res1a,res1d,chr_name,legend_size,mainline_size,backline_size,margin_space,axis_space,logPCoff,color1,color2,color3,lodthred)
{
  chr_pos <- mxmp[,1:2]
  chr_num <- length(chr_name)
  chr <- matrix(0,chr_num,1)
  pos <- matrix(0,chr_num,1)
  for(i in 1:chr_num)
  {
    temp <- numeric()
    temp <- length(which(chr_pos[,1]==i))
    if(i==1)
    {
      pos[i] <- temp
      chr[i] <- chr_pos[pos[i],2]
    }else{
      pos[i] <- pos[i-1] + temp
      chr[i] <- chr_pos[pos[i],2]
    }
  }
  
  pos_acc <- matrix(0,chr_num,1)
  for(i in 1:chr_num)
  {
    if(i==1){
      pos_acc[i] <- chr[i]
    }else{
      pos_acc[i] <- pos_acc[i-1] + chr[i]
    }
  }
  firFila <- res1a[,1:2]
  newposadda <- as.matrix(firFila[,2])
  for(i in 1:chr_num)
  {
    temp1a <- numeric()
    temp1a <- which(firFila[,1]==i)
    if(i>1)
    {
      newposadda[temp1a] <- newposadda[temp1a]+pos_acc[i-1]
    }
  }
  firFild <- res1d[,1:2]
  newposaddd <- as.matrix(firFild[,2])
  for(i in 1:chr_num)
  {
    temp1d <- numeric()
    temp1d <- which(firFild[,1]==i)
    if(i>1)
    {
      newposaddd[temp1d] <- newposaddd[temp1d]+pos_acc[i-1]
    }
  }
  if(is.null(galaxyy1)==FALSE){
    if(is.null(dim(galaxyy1))==TRUE){
      galaxyy1<-matrix(galaxyy1,1,3)
    }
    newres_pos <- galaxyy1[,2]
    res_sumpos <- pos_acc[galaxyy1[which(galaxyy1[,1]>1),1]-1] + galaxyy1[which(galaxyy1[,1]>1),2]
    newres_pos[which(galaxyy1[,1]>1)] <- res_sumpos
    pospic<-c(newres_pos)
    lodpic<-c(galaxyy1[,3])
    resdf <- data.frame(pospic,lodpic)
  }
  negloPa <- as.matrix(res1a[,3])
  negloPd <- as.matrix(res1d[,3])
  if(is.null(galaxyy1)==FALSE){
    par(mar=c(2*margin_space,2*margin_space,2*margin_space,2*margin_space)+margin_space,mgp=c(3*axis_space,axis_space,0))
    plot(pospic,lodpic,type="h",col=color1,xlab="",ylab="Logarithm of odds (LOD)",cex.axis=legend_size,cex.lab=legend_size,lwd=mainline_size,xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,max(lodpic)))
    abline(h=lodthred)
    par(new=TRUE)
    plot(newposadda,negloPa,type="l",col=color3,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPa)))
    par(new=TRUE)
    plot(newposaddd,negloPd,type="l",col=color2,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPd)))
    axis(side=4,cex.axis=legend_size)
    mtext(expression('-log'[10]*'(P)'),side=4,line=3*axis_space,cex=legend_size)
    abline(v=pos_acc,lty=2,col="gray")
  }else{
    plot(newposadda,negloPa,type="l",col=color3,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPa)))
    par(new=TRUE)
    plot(newposaddd,negloPd,type="l",col=color2,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPd)))
  }
}

#########
#-1,0,1 change a and d matrix
fix<-function(x,gen,y,kk){
  
  loglike<-function(theta){
    lambda<-exp(theta)
    logdt<-sum(log(lambda*delta+1))
    h<-1/(lambda*delta+1)
    yy<-sum(yu*h*yu)
    yx<-matrix(0,q,1)
    xx<-matrix(0,q,q)
    for(i in 1:q){
      yx[i]<-sum(yu*h*xu[,i])
      for(j in 1:q){
        xx[i,j]<-sum(xu[,i]*h*xu[,j])
      }
    }
    if(abs(min(eigen(xx)$values))<1e-6)
      loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx+diag(ncol(xx))*0.01)%*%yx)-0.5*log(det(xx+diag(ncol(xx))*0.01))
    else
      loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
    return(-loglike)
  }
  
  fixed<-function(lambda){
    h<-1/(lambda*delta+1)
    yy<-sum(yu*h*yu)
    yx<-matrix(0,q,1)
    xx<-matrix(0,q,q)
    for(i in 1:q){
      yx[i]<-sum(yu*h*xu[,i])
      for(j in 1:q){
        xx[i,j]<-sum(xu[,i]*h*xu[,j])
      }
    }
    if(abs(min(eigen(xx)$values))<1e-6)
      beta<-solve(xx+diag(ncol(xx))*0.01,yx)
    else
      beta<-solve(xx,yx)
    if(abs(min(eigen(xx)$values))<1e-6)
      sigma2<-(yy-t(yx)%*%solve(xx+diag(ncol(xx))*0.01)%*%yx)/(n-q)
    else
      sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
    sigma2<-drop(sigma2)
    if(abs(min(eigen(xx)$values))<1e-6)
      vertue<-solve(xx+diag(ncol(xx))*0.01)
    else
      vertue<-solve(xx)
    var<-diag(vertue*sigma2)
    stderr<-sqrt(var)
    return(c(beta,stderr,sigma2))
  }
  qq<-eigen(as.matrix(kk))
  delta<-qq[[1]]
  uu<-qq[[2]]
  qx<-ncol(x)
  n<-length(y)
  yu<-t(uu)%*%y
  tempx<-x
  
  cl.cores <- detectCores()
  if (cl.cores<=2) {
    cl.cores<-1
  }else if(cl.cores>2){
    cl.cores <- detectCores()-1
  }
  cll <- makeCluster(cl.cores)
  registerDoParallel(cll)
  
  i<-numeric()
  parmm<-foreach(i=1:nrow(gen),.combine=rbind)%dopar%{
    x<-tempx
    z<-gen[i,3:(ncol(gen))]
    qz<-ncol(z)
    x<-cbind(x,z)
    q<-ncol(x)
    xu<-t(uu)%*%x
    theta<-0
    parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-50,upper=10)
    lambda<-exp(parm$par)
    conv<-parm$convergence
    fn1<-parm$value
    fn0<-loglike(-Inf)
    lrt<-2*(fn0-fn1)
    hess<-parm$hessian
    parmfix<-fixed(lambda)
    beta<-parmfix[1:q]
    stderr<-parmfix[(q+1):(2*q)]
    sigma2<-parmfix[2*q+1]
    poly.lod<-lrt/4.61
    poly.p<-1-pchisq(lrt,1)
    sigma2g<-lambda*sigma2
    g<-beta[-c(1:qx)]
    g.err<-stderr[-c(1:qx)]
    b<-beta[c(1:qx)]
    b.err<-stderr[c(1:qx)]
    wald<-g^2/g.err^2
    p<-1-pchisq(wald,1)
    parmm<-(c(b[1],sigma2,lambda,sigma2g,poly.lod,poly.p,g,g.err,wald,p))
  }
  stopCluster(cll)
  return(parmm)
}


random<-function(fx,gen,phe,kk)
{
  mixed<-function(x,y,kk){
    
    loglike<-function(theta){
      lambda<-exp(theta)
      logdt<-sum(log(lambda*delta+1))
      h<-1/(lambda*delta+1)
      yy<-sum(yu*h*yu)
      yx<-matrix(0,q,1)
      xx<-matrix(0,q,q)
      for(i in 1:q){
        yx[i]<-sum(yu*h*xu[,i])
        for(j in 1:q){
          xx[i,j]<-sum(xu[,i]*h*xu[,j])
        }
      }
      loglike<- -0.5*logdt-0.5*(n-q)*log(yy-t(yx)%*%solve(xx)%*%yx)-0.5*log(det(xx))
      return(-loglike)
    }
    
    fixed<-function(lambda){
      h<-1/(lambda*delta+1)
      yy<-sum(yu*h*yu)
      yx<-matrix(0,q,1)
      xx<-matrix(0,q,q)
      for(i in 1:q){
        yx[i]<-sum(yu*h*xu[,i])
        for(j in 1:q){
          xx[i,j]<-sum(xu[,i]*h*xu[,j])
        }
      }
      beta<-solve(xx,yx)
      sigma2<-(yy-t(yx)%*%solve(xx)%*%yx)/(n-q)
      sigma2<-drop(sigma2)
      var<-diag(solve(xx)*sigma2)
      stderr<-sqrt(var)
      return(c(beta,stderr,sigma2))
    }
    
    qq<-eigen(kk)
    delta<-qq[[1]]
    uu<-qq[[2]]
    q<-ncol(x)
    n<-ncol(kk)
    vp<-var(y)
    yu<-t(uu)%*%y
    xu<-t(uu)%*%x
    theta<-0
    parm<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-50,upper=10)
    lambda<-exp(parm$par)
    conv<-parm$convergence
    fn1<-parm$value
    fn0<-loglike(-Inf)
    lrt<-2*(fn0-fn1)
    hess<-parm$hessian
    parmfix<-fixed(lambda)
    beta<-parmfix[1:q]
    stderr<-parmfix[(q+1):(2*q)]
    sigma2<-parmfix[2*q+1]
    lod<-lrt/4.61
    p_value<-1-pchisq(lrt,1)
    sigma2g<-lambda*sigma2
    goodness<-(vp-sigma2)/vp
    par<-data.frame(lrt,beta,stderr,sigma2,lambda,sigma2g,lod,p_value)
    return(par)
  }
  
  
  loglike<-function(theta){
    xi<-exp(theta)
    tmp0<-zz*xi+1
    tmp<-xi*solve(tmp0)
    yHy<-yy-t(zy)%*%tmp%*%zy
    yHx<-yx-zx%*%tmp%*%zy
    xHx<-xx-zx%*%tmp%*%t(zx)
    logdt2<-log(det(tmp0))
    loglike<- -0.5*logdt2-0.5*(n-s)*log(yHy-t(yHx)%*%solve(xHx)%*%yHx)-0.5*log(det(xHx))
    return(-loglike)
  }
  
  fixed<-function(xi){
    tmp0<-zz*xi+diag(1)
    tmp<-xi*solve(tmp0)
    yHy<-yy-t(zy)%*%tmp%*%zy
    yHx<-yx-zx%*%tmp%*%zy
    xHx<-xx-zx%*%tmp%*%t(zx)
    zHy<-zy-zz%*%tmp%*%zy
    zHx<-zx-zx%*%tmp%*%zz
    zHz<-zz-zz%*%tmp%*%zz
    beta<-solve(xHx,yHx)
    tmp2<-solve(xHx)
    sigma2<-(yHy-t(yHx)%*%tmp2%*%yHx)/(n-s)
    gamma<-xi*zHy-xi*t(zHx)%*%tmp2%*%yHx
    var<-abs((xi*diag(1)-xi*zHz*xi)*as.numeric(sigma2))
    stderr<-sqrt(diag(var))
    result<-list(gamma,stderr,beta,sigma2)
    return(result)
  }
  name<-gen[,1:2]
  gen<-gen[,3:(ncol(gen))]
  gen<-t(gen)
  n<-nrow(gen)
  m<-ncol(gen)
  x<-fx
  
  s<-ncol(x)
  kk<-as.matrix(kk)
  qq<-eigen(kk)
  delta<-qq[[1]]
  uu<-qq[[2]]
  xu<-t(uu)%*%x
  y<-as.matrix(phe)
  parm<-mixed(x=x,y=y,kk=kk)
  lambda<-parm$lambda[1]
  h<-1/(delta*lambda+1)
  yu<-t(uu)%*%y
  xx<-matrix(0,s,s)
  for(i in 1:s){
    for(j in 1:s){
      xx[i,j]<-sum(xu[,i]*h*xu[,j])
    }
  }
  yy<-sum(yu*h*yu)
  yx<-matrix(0,s,1)
  for(i in 1:s){
    yx[i]<-sum(yu*h*xu[,i])
  }
  
  cl.cores <- detectCores()
  if (cl.cores<=2) {
    cl.cores<-1
  }else if(cl.cores>2){
    cl.cores <- detectCores()-1
  }
  cll <- makeCluster(cl.cores)
  registerDoParallel(cll)
  
  k<-numeric()
  parms<-foreach(k=1:m,.combine=rbind)%dopar%{
    z<-as.matrix(gen[,k])
    zu<-t(uu)%*%z
    zy<-as.matrix(sum(yu*h*zu))
    zz<-as.matrix(sum(zu*h*zu))
    zx<-matrix(0,s,1)
    for(i in 1:s){
      zx[i]<-sum(xu[,i]*h*zu)
    }
    theta<-c(0)
    par<-optim(par=theta,fn=loglike,hessian = TRUE,method="L-BFGS-B",lower=-10,upper=10)
    xi<-exp(par$par)
    conv<-par$convergence
    fn1<-par$value
    hess<-par$hessian
    parmfix<-fixed(xi)
    gamma<-parmfix[[1]]
    stderr<-parmfix[[2]]
    beta<-parmfix[[3]]
    sigma2<-parmfix[[4]]
    lambda<-xi
    sigma2g<-lambda*sigma2
    fn0<-loglike(-Inf)
    lrt<-2*(fn0-fn1)
    p_lrt<-1-pchisq(lrt,1)
    wald<-(gamma/stderr)^2
    p_wald<-1-pchisq(wald,1)
    parm0<-c(k,name[k,1],name[k,2],beta[1],sigma2,sigma2g,gamma,stderr,wald,p_wald)
  }
  stopCluster(cll)
  return(parms)
}


multinormal<-function(y,mean,sigma)
{
  pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
  return (pdf_value)
}
#LOD value test
likelihood<-function(xxn,xxx,yn)
{
  nq<-ncol(xxx)
  ns<-nrow(yn)
  at1<-0
  ww1<-1:ncol(xxx)
  ww1<-as.matrix(ww1)
  at1<-dim(ww1)[1]
  lod<-matrix(rep(0,nq),nq,1)
  if(at1>0.5)
    ad<-cbind(xxn,xxx[,ww1])
  else
    ad<-xxn
  #if(abs(det(crossprod(ad,ad)))<1e-6)
  if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6)
    bb<-solve(crossprod(ad,ad)+diag(ncol(ad))*0.01)%*%crossprod(ad,yn)
  else
    bb<-solve(crossprod(ad,ad))%*%crossprod(ad,yn)
  vv1<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns);
  ll1<-sum(log(abs(multinormal(yn,ad%*%bb,vv1))))
  
  sub<-1:ncol(ad);
  if(at1>0.5)
  {
    for(i in 1:at1)
    {
      ij<-which(sub!=sub[(i+ncol(xxn))])
      ad1<-ad[,ij]
      #if(abs(det(crossprod(ad1,ad1)))<1e-6)
      if(abs(min(eigen(crossprod(ad1,ad1))$values))<1e-6)
        bb1<-solve(crossprod(ad1,ad1)+diag(ncol(ad1))*0.01)%*%crossprod(ad1,yn)
      else
        bb1<-solve(crossprod(ad1,ad1))%*%crossprod(ad1,yn)
      vv0<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns);
      ll0<-sum(log(abs(multinormal(yn,ad1%*%bb1,vv0))))
      lod[ww1[i]]<--2.0*(ll0-ll1)/(2.0*log(10))
    }
  }
  return (lod)
}

#2010 EM_Bayes
ebayes_EM<-function(x,z,y)
{
  n<-nrow(z);k<-ncol(z)
  if(abs(min(eigen(crossprod(x,x))$values))<1e-6)
    b<-solve(crossprod(x,x)+diag(ncol(x))*1e-8)%*%crossprod(x,y)
  else
    b<-solve(crossprod(x,x))%*%crossprod(x,y)
  v0<-as.numeric(crossprod((y-x%*%b),(y-x%*%b))/n)
  u<-matrix(rep(0,k),k,1)
  v<-matrix(rep(0,k),k,1)
  s<-matrix(rep(0,k),k,1)
  for(i in 1:k)
  {
    zz<-z[,i]
    s[i]<-((crossprod(zz,zz))^(-1))*v0
    u[i]<-s[i]*crossprod(zz,(y-x%*%b))/v0
    v[i]<-u[i]^2+s[i]
  }
  vv<-matrix(rep(0,n*n),n,n);
  for(i in 1:k)
  {
    zz<-z[,i]
    vv=vv+tcrossprod(zz,zz)*v[i]
  }
  vv<-vv+diag(n)*v0
  iter<-0;err<-1000;iter_max<-100;err_max<-1e-8
  tau<-0;omega<-0
  while((iter<iter_max)&&(err>err_max))
  {
    iter<-iter+1
    v01<-v0
    v1<-v
    b1<-b
    vi<-solve(vv)
    xtv<-crossprod(x,vi)
    if(ncol(x)==1)
    {
      b<-((xtv%*%x)^(-1))*(xtv%*%y)
    }else
    {
      if(abs(min(eigen(xtv%*%x)$values))<1e-6){
        b<-solve((xtv%*%x)+diag(ncol(x))*1e-8)%*%(xtv%*%y)
      }
      else{
        b<-solve(xtv%*%x)%*%(xtv%*%y)
      }
    }
    r<-y-x%*%b
    ss<-matrix(rep(0,n),n,1)
    for(i in 1:k)
    {
      zz<-z[,i]
      zztvi<-crossprod(zz,vi)
      u[i]<-v[i]*zztvi%*%r
      s[i]<-v[i]*(1-zztvi%*%zz*v[i])
      v[i]<-(u[i]^2+s[i]+omega)/(tau+3)
      ss<-ss+zz*u[i]
    }
    v0<-as.numeric(crossprod(r,(r-ss))/n)
    vv<-matrix(rep(0,n*n),n,n)
    for(i in 1:k)
    {
      zz<-z[,i]
      vv<-vv+tcrossprod(zz,zz)*v[i]
    }
    vv<-vv+diag(n)*v0
    err<-(crossprod((b1-b),(b1-b))+(v01-v0)^2+crossprod((v1-v),(v1-v)))/(2+k)
    beta<-t(b)
    sigma2<-v0
  }
  wang<-matrix(rep(0,k),k,1)
  for (i in 1:k){
    stderr<-sqrt(s[i]+1e-20)
    t<-abs(u[i])/stderr
    f<-t*t
    p<-1-pchisq(f,1)
    wang[i]<-p
  }
  return (wang)
}


#####k kinship start################
kinship.F2 <- function(gen){
  X1<-as.matrix(gen)
  kk<-X1%*%t(X1)
  cc<-mean(diag(kk))
  kk<-kk/cc
  return(kk)
}
###################################################
#K kinship end
#############################################
###mixed.vars start
#Mixednew can fit the models involving mutiple variace matrixs. 

#input
##dataframe<-cbind(y,w)
#y:n*1
#w:n*c,include intercept,fixed effect
#kinship<-list(K.a,K.d,...)
#K.a,K.d,...,:n*n,matrix

#output
#RR$tau.kk:sigma.a2,sigma.d2,...,sigma.e2

mixed.vars<-function(dataframe,kinship,optim.speed=TRUE){
  #data
  d<-dataframe
  y<-d[,1]
  x<-as.matrix(d[,-1,drop=FALSE])
  kk<-kinship
  indi<-nrow(d)
  r<-ncol(x)
  num.kk<-length(kk)
  
  cal.xx<-function(x,y,H=NA){
    x<-as.matrix(x)
    y<-as.matrix(y)
    nr<-ncol(x)
    nc<-ncol(y)
    n0<-nrow(x)
    if (is.na(H)) H<-diag(n0)
    xx<-t(x)%*%H%*%y
    xx<-as.matrix(xx)
    return(xx)
  }
  
  
  ##likelihood function
  ##loglike
  loglike<-function(theta){
    i0<-num.kk+1
    V<-diag(indi)*as.numeric(theta[i0])
    for ( i in 1:num.kk){
      lambda<-as.numeric(theta[i])
      V<-V+kk[[i]]*lambda
    }
    
    v.inv<-solve(V,tol=1e-50)
    xx<-cal.xx(x=x,y=x,H=v.inv)
    xy<-cal.xx(x=x,y=y,H=v.inv)
    yy<-cal.xx(x=y,y=y,H=v.inv)
    
    d1<-unlist(determinant(V))
    d1<-d1[1]
    d2<-unlist(determinant(xx))
    d2<-d2[1]
    p<-yy-t(xy)%*%solve(xx)%*%xy
    
    
    logvalue<- -0.5*(d1+d2+p)
    return(-logvalue)
  }
  
  ## 
  ##parms involved in each kinship,sg2.1,sg2.2,...,se2. 
  if (optim.speed){
    theta<-rep(1,num.kk+1)
    parm<-optim(par=theta,fn=loglike,hessian=TRUE,method="L-BFGS-B",lower=1e-10,upper=1e+10)
    tau.kk<-parm$par
    conv<-parm$convergence
    fn1<-parm$value
    hess<-parm$hessian
  }else{
    intervals<-c(1e-10,10,5e+02,1e+05)
    intials<-rep(1,num.kk+1)
    ni<-length(intervals)-1
    thetas<-list()
    for ( i in 1:ni){
      theta<-intials
      parms<-optim(par=theta,fn=loglike,hessian=TRUE,method="L-BFGS-B",lower=intervals[1],upper=intervals[i+1])
      theta<-parms$par
      parms<-optim(par=theta,fn=loglike,hessian=TRUE,method="L-BFGS-B",lower=1e-10,upper=1e+10)
      fi<-parms$par
      thetas[[i]]<-fi
    }
    fns<-sapply(1:ni,function(i){
      theta<-thetas[[i]]
      fn<-loglike(theta)
      return(fn)
    } )
    ii<-which.min(fns)
    tau.kk<-thetas[[ii]]
    fn1<-fns[ii]
  }
  
  ###  
  ###beta,var(beta),residual variance and blup estimation 
  i0<-num.kk+1
  V<-diag(indi)*as.numeric(tau.kk[i0])
  G<-matrix(0,indi,indi)
  for ( i in 1:num.kk){
    sg2<-as.numeric(tau.kk[i])
    V<-V+kk[[i]]*sg2
    G<-G+kk[[i]]*sg2
  }
  v.inv<-solve(V)
  xx<-cal.xx(x=x,y=x,H=v.inv)
  xy<-cal.xx(x=x,y=y,H=v.inv)
  yy<-cal.xx(x=y,y=y,H=v.inv)
  
  beta<-solve(xx,xy,tol=1e-50)
  v.beta<-solve(xx,tol=1e-50)      
  blup.sum<-G%*%v.inv%*%(y-x%*%beta)
  #Estimated weights
  yhat<-y-x%*%as.matrix(beta)
  yrandom<-NULL
  for (i in 1:num.kk){
    V0<-kk[[i]]*as.numeric(tau.kk[i])
    y0<-V0%*%v.inv%*%yhat
    yrandom<-cbind(yrandom,y0)
  }
  #
  xx<-t(yrandom)%*%yrandom
  xy<-t(yrandom)%*%yhat
  ww<-solve(xx,xy,tol=1e-50)
  ww<-as.vector(ww)
  #   
  
  RR<-list(beta=beta,v.beta=v.beta,random=blup.sum,se2=tau.kk[i0],tau.kk=tau.kk,weights=ww,lrt=fn1)
  ###
  return(RR)
}
#mixed.vars end
################################################
#Y=w*alpha+Z*u+e start
################################################
#likelihood
emma.eigen.L <- function(Z,K,complete=TRUE) {
  if ( is.null(Z) ) {
    return(emma.eigen.L.wo.Z(K))
  }
  else {
    return(emma.eigen.L.w.Z(Z,K,complete))
  }
}
#likelihood
emma.eigen.L.wo.Z <- function(K) {
  eig <- eigen(K,symmetric=TRUE)
  return(list(values=eig$values,vectors=eig$vectors))
}
#likelihood
emma.eigen.L.w.Z <- function(Z,K,complete=TRUE) {
  if ( complete == FALSE ) {
    vids <- colSums(Z)>0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  eig <- eigen(K%*%crossprod(Z,Z),symmetric=FALSE,EISPACK=TRUE)
  return(list(values=eig$values,vectors=qr.Q(qr(Z%*%eig$vectors),complete=TRUE)))
}

#restricted likelihood
emma.eigen.R <- function(Z,K,X,complete=TRUE) {
  if ( ncol(X) == 0 ) {
    return(emma.eigen.L(Z,K))
  }
  else if ( is.null(Z) ) {
    return(emma.eigen.R.wo.Z(K,X))
  }
  else {
    return(emma.eigen.R.w.Z(Z,K,X,complete))
  }
}
#restricted likelihood
emma.eigen.R.wo.Z <- function(K, X) {
  n <- nrow(X)
  q <- ncol(X)
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}
#restricted likelihood
emma.eigen.R.w.Z <- function(Z, K, X, complete = TRUE) {
  if ( complete == FALSE ) {
    vids <-  colSums(Z) > 0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  n <- nrow(Z)
  t <- ncol(Z)
  q <- ncol(X)
  
  SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
  eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE,EISPACK=TRUE)
  if ( is.complex(eig$values) ) {
    eig$values <- Re(eig$values)
    eig$vectors <- Re(eig$vectors)    
  }
  qr.X <- qr.Q(qr(X))
  return(list(values=eig$values[1:(t-q)],
              vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:(t-q)],qr.X)),
                           complete=TRUE)[,c(1:(t-q),(t+1):n)]))   
}

emma.delta.ML.LL.wo.Z <- function(logdelta, lambda, etas, xi) {
  n <- length(xi)
  delta <- exp(logdelta)
  return( 0.5*(n*(log(n/(2*pi))-1-log(sum((etas*etas)/(delta*lambda+1))))-sum(log(delta*xi+1))) )  
}

emma.delta.ML.LL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
  delta <- exp(logdelta)
  return( 0.5*(n*(log(n/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*xi.1+1)) ))
}

emma.delta.ML.dLL.wo.Z <- function(logdelta, lambda, etas, xi) {
  n <- length(xi)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- delta*lambda+1
  return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(xi/(delta*xi+1))) )
}

emma.delta.ML.dLL.w.Z <- function(logdelta, lambda, etas.1, xi.1, n, etas.2.sq ) {
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- delta*lambda+1
  return( 0.5*(n*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(xi.1/(delta*xi.1+1))) )
}

emma.delta.REML.LL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas*etas/(delta*lambda+1))))-sum(log(delta*lambda+1))) )
}

emma.delta.REML.LL.w.Z <- function(logdelta, lambda, etas.1, n, t, etas.2.sq ) {
  tq <- length(etas.1)
  nq <- n - t + tq#(n-t):the number of eigvalue 1;tq:the number of components of etas.1
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*lambda+1))) ) 
}

emma.delta.REML.dLL.wo.Z <- function(logdelta, lambda, etas) {
  nq <- length(etas)
  delta <- exp(logdelta)
  etasq <- etas*etas
  ldelta <- delta*lambda+1
  return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/sum(etasq/ldelta)-sum(lambda/ldelta)) )
}

emma.delta.REML.dLL.w.Z <- function(logdelta, lambda, etas.1, n, t1, etas.2.sq ) {
  t <- t1
  tq <- length(etas.1)
  nq <- n - t + tq
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- delta*lambda+1
  return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(lambda/ldelta) ))
}

emma.MLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                     esp=1e-10, eig.L = NULL, eig.R = NULL)
{
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)

  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)
  
  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(ML=0,delta=0,ve=0,vg=0))
  }
  
  if ( is.null(Z) ) {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.wo.Z(K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    
    Lambdas.1<-matrix(eig.R$values,n-q,m)    
    Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE)+1
    Xis.1<-matrix(eig.L$values,n,m)
    Xis <- Xis.1* matrix(delta,n,m,byrow=TRUE)+1
 
    Etasq <- matrix(etas*etas,n-q,m)
   dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Xis.1/Xis))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.R$values,etas,eig.L$values))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas, xi=eig.L$values)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
      }
    }
    #    optdelta <- exp(optlogdelta)
  }
  else {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.w.Z(Z,K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    
    m <- length(logdelta)
    delta <- exp(logdelta)
    
    Lambdas.1<-matrix(eig.R$values,t-q,m)
    Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
    
    Xis.1<-matrix(eig.L$values,t,m)
    Xis <- Xis.1 * matrix(delta,t,m,byrow=TRUE) + 1
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    
    dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Xis.1/Xis))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
      }
    }
    #    optdelta <- exp(optlogdelta)
  }
  #handler of grids with NaN log
  optLL=replaceNaN(optLL)  #20160728
  
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/n    
  }
  else {
    maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
  }
  maxvg <- maxve*maxdelta
  
  return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
  
}


emma.REMLE <- function(y, X, K, Z=NULL, ngrids=100, llim=-10, ulim=10,
                       esp=1e-10, eig.L = NULL, eig.R = NULL) {
  n <- length(y)
  t <- nrow(K)
  q <- ncol(X)

  stopifnot(ncol(K) == t)
  stopifnot(nrow(X) == n)
  
  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(REML=0,delta=0,ve=0,vg=0))
  }
  
  if ( is.null(Z) ) {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    
    Lambdas.1<-matrix(eig.R$values,n-q,m)
    Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE) + 1
    Etasq <- matrix(etas*etas,n-q,m)
    dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Lambdas.1/Lambdas))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
      }
    }
    #    optdelta <- exp(optlogdelta)
  }
  else {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    etas.1 <- etas[1:(t-q)]
    etas.2 <- etas[(t-q+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    
    Lambdas.1 <- matrix(eig.R$values,t-q,m) 
    Lambdas <- Lambdas.1 * matrix(delta,t-q,m,byrow=TRUE) + 1
    Etasq <- matrix(etas.1*etas.1,t-q,m)
    
    dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Lambdas.1/Lambdas))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(llim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z(ulim,eig.R$values,etas.1,n,t,etas.2.sq))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, t1=t, etas.2.sq = etas.2.sq )
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z(r$root,eig.R$values, etas.1, n, t, etas.2.sq ))
      }
    }
    #    optdelta <- exp(optlogdelta)
  }  
  #handler of grids with NaN log
  optLL=replaceNaN(optLL)  #20160728
  
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  
  maxLL <- max(optLL)
  
  if ( is.null(Z) ) {
    maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/(n-q)    
  }
  else {
    maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/(n-q)
  }
  maxvg <- maxve*maxdelta
  return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg))
}

################################################
#Y=w*alpha+Z*u+e end

#change to Y_c=W_c*alpha+X_c*beta+e_c,K=1 start
##################################################################################
#change,k=1,Z=X_c,X=W_c:


################################################
#likelihood: 
emma.eigen.L.c <- function(Z,K,complete=TRUE) {
  if ( is.null(Z) ) {
    return(emma.eigen.L.wo.Z.c(K))
  }
  else {
    return(emma.eigen.L.w.Z.c(Z,K,complete))
  }
}
#likelihood
emma.eigen.L.wo.Z.c <- function(K) {
  eig <- eigen(K,symmetric=TRUE)
  return(list(values=eig$values,vectors=eig$vectors))
}

#likelihood
emma.eigen.L.w.Z.c <- function(Z,K,complete=TRUE) {
  if ( complete == FALSE ) {
    vids <- colSums(Z)>0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  eig <- eigen(K%*%crossprod(Z,Z),symmetric=TRUE,EISPACK=TRUE)
  return(list(values=eig$values,vectors=qr.Q(qr(Z%*%eig$vectors),complete=TRUE)))
}

#restricted likelihood
emma.eigen.R.c <- function(Z,K,X,complete=TRUE) {
  if ( ncol(X) == 0 ) {
    return(emma.eigen.L.c(Z,K))
  }
  else if ( is.null(Z) ) {
    return(emma.eigen.R.wo.Z.c(K,X))
  }
  else {
    return(emma.eigen.R.w.Z.c(Z,K,X,complete))
  }
}
#restricted likelihood
emma.eigen.R.wo.Z.c <- function(K, X) {
 
  if(is.matrix(X)) {    
    n<-nrow(X)   
    q<-ncol(X) 
  } 
  else{   
    n<-length(X) 
    q<-1 
  }
  S <- diag(n)-X%*%solve(crossprod(X,X))%*%t(X)
  eig <- eigen(S%*%(K+diag(1,n))%*%S,symmetric=TRUE)
  stopifnot(!is.complex(eig$values))
  return(list(values=eig$values[1:(n-q)]-1,vectors=eig$vectors[,1:(n-q)]))
}

#restricted likelihood
emma.eigen.R.w.Z.c <- function(Z, K, X, complete = TRUE) {
  if ( complete == FALSE ) {
    vids <-  colSums(Z) > 0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  if(!is.matrix(Z)) n<-length(Z) 
  t <-1 
  if(is.matrix(X)) {
    q<-ncol(X)
  }    
  else  {
    q<-1
  }

  SZ <- Z - X%*%solve(crossprod(X,X))%*%crossprod(X,Z)
  eig <- eigen(K%*%crossprod(Z,SZ),symmetric=FALSE,EISPACK=TRUE)
  if ( is.complex(eig$values) ) {
    eig$values <- Re(eig$values)
    eig$vectors <- Re(eig$vectors)    
  }
  qr.X <- qr.Q(qr(X))
  return(list(values=eig$values[1],
              vectors=qr.Q(qr(cbind(SZ%*%eig$vectors[,1:t],qr.X)),
                           complete=TRUE)[,c(1:t,(t+q+1):n)]))  

}
######################################################################################
emma.delta.REML.LL.w.Z.c <- function(logdelta, lambda, etas.1, n, q, etas.2.sq ) {
   nq<-n-q 
  delta <-  exp(logdelta)
  return( 0.5*(nq*(log(nq/(2*pi))-1-log(sum(etas.1*etas.1/(delta*lambda+1))+etas.2.sq))-sum(log(delta*lambda+1))) ) 
}

emma.delta.REML.dLL.w.Z.c <- function(logdelta, lambda, etas.1, n, q1, etas.2.sq ) {
  q<-q1
  nq<-n-q
  delta <- exp(logdelta)
  etasq <- etas.1*etas.1
  ldelta <- delta*lambda+1
  return( 0.5*(nq*sum(etasq*lambda/(ldelta*ldelta))/(sum(etasq/ldelta)+etas.2.sq)-sum(lambda/ldelta) ))
}
###########################

emma.MLE.c <- function(y, X, K=1, Z=NULL, ngrids=100, llim=-10, ulim=10,
                       esp=1e-10, eig.L = NULL, eig.R = NULL)#Z=X_c,K=1,X=W_c
{
  if(is.matrix(y)){
    n<-nrow(y)
  }
  else{
    n <- length(y)
  }
  
  t <- 1#,K=1,t=1
  if ( is.matrix(X) ){
    q <- ncol(X)
    stopifnot(nrow(X)==n)
  }
  else{#X:n*1 vector
    q <- 1
    stopifnot(length(X)==n)
  }
  
  stopifnot(K == 1)
  
  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(ML=0,delta=0,ve=0,vg=0))
  }
  
  if ( is.null(Z) ) {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.wo.Z.c(K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z.c(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    
    Lambdas.1<-matrix(eig.R$values,n-q,m)    
    Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE)+1
    Xis.1<-matrix(eig.L$values,n,m)
    Xis <- Xis.1* matrix(delta,n,m,byrow=TRUE)+1
    Etasq <- matrix(etas*etas,n-q,m)
    dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Xis.1/Xis))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(llim,eig.R$values,etas,eig.L$values))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.wo.Z(ulim,eig.R$values,etas,eig.L$values))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.ML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas, xi=eig.L$values)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.ML.LL.wo.Z(r$root,eig.R$values, etas, eig.L$values))
      }
    }
    #    optdelta <- exp(optlogdelta)
  }
  else {
    if ( is.null(eig.L) ) {
      eig.L <- emma.eigen.L.w.Z.c(Z,K)
    }
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z.c(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
   etas.1<-etas[1:t]
    etas.2<-etas[(t+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    
    m <- length(logdelta)
    delta <- exp(logdelta)
    
    Lambdas.1<-matrix(eig.R$values,t,m)
    Lambdas <- Lambdas.1 * matrix(delta,t,m,byrow=TRUE) + 1
    
    Xis.1<-matrix(eig.L$values,t,m)
    Xis <- Xis.1 * matrix(delta,t,m,byrow=TRUE) + 1
    Etasq <- matrix(etas.1*etas.1,t,m)
     dLL <- 0.5*delta*(n*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Xis.1/Xis))
     optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(llim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.ML.LL.w.Z(ulim,eig.R$values,etas.1,eig.L$values,n,etas.2.sq))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.ML.dLL.w.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, xi.1=eig.L$values, n=n, etas.2.sq = etas.2.sq )
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.ML.LL.w.Z(r$root,eig.R$values, etas.1, eig.L$values, n, etas.2.sq ))
      }
    }
    #    optdelta <- exp(optlogdelta)
  }
  
  #handler of grids with NaN log
  optLL=replaceNaN(optLL)  #20160728
  
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  
  maxLL <- max(optLL)
  if ( is.null(Z) ) {
    maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/n    
  }
  else {
    maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/n
  }
  maxvg <- maxve*maxdelta
  
  return (list(ML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg,U_R=eig.R$vectors,etas.1=etas.1,etas=etas,lambda=eig.R$values))
}


emma.REMLE.c <- function(y, X, K=1, Z=NULL, ngrids=100, llim=-10, ulim=10,
                         esp=1e-10, eig.L = NULL, eig.R = NULL) #K=1,X=W_c,Z=X_c
{
  if(is.matrix(y)){
    n<-nrow(y)
  }
  else{
    n <- length(y)
  }
  
  t <- 1#,K=1,t=1
  if ( is.matrix(X) ){
    q <- ncol(X)
  }
  else{#X:n*1 vector
    q <- 1
  }
  stopifnot(K == 1)
  stopifnot(nrow(X) == n)
  
  if ( det(crossprod(X,X)) == 0 ) {
    warning("X is singular")
    return (list(REML=0,delta=0,ve=0,vg=0))
  }
  
  if ( is.null(Z) ) {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.wo.Z.c(K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)
    
    Lambdas.1<-matrix(eig.R$values,n-q,m)
    Lambdas <- Lambdas.1 * matrix(delta,n-q,m,byrow=TRUE) + 1
     Etasq <- matrix(etas*etas,n-q,m)
    
    dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/colSums(Etasq/Lambdas)-colSums(Lambdas.1/Lambdas))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(llim,eig.R$values,etas))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.wo.Z(ulim,eig.R$values,etas))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.wo.Z, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas=etas)
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.wo.Z(r$root,eig.R$values, etas))
      }
    }
    #    optdelta <- exp(optlogdelta)
  }
  else {
    if ( is.null(eig.R) ) {
      eig.R <- emma.eigen.R.w.Z.c(Z,K,X)
    }
    etas <- crossprod(eig.R$vectors,y)
 
    etas.1 <- etas[1:t]
    etas.2 <- etas[(t+1):(n-q)]
    etas.2.sq <- sum(etas.2*etas.2)
    
    logdelta <- (0:ngrids)/ngrids*(ulim-llim)+llim
    m <- length(logdelta)
    delta <- exp(logdelta)

    Lambdas.1 <- matrix(eig.R$values,t,m) 
    Lambdas <- Lambdas.1 * matrix(delta,t,m,byrow=TRUE) + 1
    Etasq <- matrix(etas.1*etas.1,t,m)
    
    dLL <- 0.5*delta*((n-q)*colSums(Etasq*Lambdas.1/(Lambdas*Lambdas))/(colSums(Etasq/Lambdas)+etas.2.sq)-colSums(Lambdas.1/Lambdas))
    optlogdelta <- vector(length=0)
    optLL <- vector(length=0)
    if ( dLL[1] < esp ) {
      optlogdelta <- append(optlogdelta, llim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z.c(llim,eig.R$values,etas.1,n,q,etas.2.sq))
    }
    if ( dLL[m-1] > 0-esp ) {
      optlogdelta <- append(optlogdelta, ulim)
      optLL <- append(optLL, emma.delta.REML.LL.w.Z.c(ulim,eig.R$values,etas.1,n,q,etas.2.sq))
    }
    
    for( i in 1:(m-1) )
    {
      if ( ( dLL[i]*dLL[i+1] < 0-esp*esp ) && ( dLL[i] > 0 ) && ( dLL[i+1] < 0 ) ) 
      {
        r <- uniroot(emma.delta.REML.dLL.w.Z.c, lower=logdelta[i], upper=logdelta[i+1], lambda=eig.R$values, etas.1=etas.1, n=n, q1=q, etas.2.sq = etas.2.sq )#t1=t change to q1=q#have revised 20140830
        optlogdelta <- append(optlogdelta, r$root)
        optLL <- append(optLL, emma.delta.REML.LL.w.Z.c(r$root,eig.R$values, etas.1, n, q, etas.2.sq ))
      }
    }
    #    optdelta <- exp(optlogdelta)
  }  
  
  #handler of grids with NaN log
  optLL=replaceNaN(optLL)  #20160728
  
  maxdelta <- exp(optlogdelta[which.max(optLL)])
  
  maxLL <- max(optLL)
  
  if ( is.null(Z) ) {
    maxve <- sum(etas*etas/(maxdelta*eig.R$values+1))/(n-q)    
  }
  else {
    maxve <- (sum(etas.1*etas.1/(maxdelta*eig.R$values+1))+etas.2.sq)/(n-q)
  }
  maxvg <- maxve*maxdelta
   return (list(REML=maxLL,delta=maxdelta,ve=maxve,vg=maxvg,U_R=eig.R$vectors,etas.1=etas.1,etas=etas,lambda=eig.R$values))
}

###################################################################
#################################################
emma.maineffects.B<-function(Z=NULL,K,deltahat.g,complete=TRUE){
  if( is.null(Z) ){
    return(emma.maineffects.B.Zo(K,deltahat.g))
  }
  else{
    return(emma.maineffects.B.Z(Z,K,deltahat.g,complete))
  }
}

#####
emma.maineffects.B.Zo <-function(K,deltahat.g){
  t <- nrow(K)
  stopifnot(ncol(K) == t)
  
  B<-deltahat.g*K+diag(1,t)
  eig<-eigen(B,symmetric=TRUE)
  qr.B<-qr(B)
  q<-qr.B$rank
  
  stopifnot(!is.complex(eig$values))
  
  A<-diag(1/sqrt(eig$values[1:q]))
  Q<-eig$vectors[,1:q]
  C<-Q%*%A%*%t(Q)
  return(list(mC=C,Q=Q,A=A))
}

emma.maineffects.B.Z <- function(Z,K,deltahat.g,complete=TRUE){
  if ( complete == FALSE ) {
    vids <- colSums(Z)>0
    Z <- Z[,vids]
    K <- K[vids,vids]
  }
  
  n <- nrow(Z)  
  B <- deltahat.g*Z%*%K%*%t(Z)+diag(1,n)
  eig <- eigen(B,symmetric=TRUE,EISPACK=TRUE)
  qr.B<-qr(B)
  q<-qr.B$rank
  
  stopifnot(!is.complex(eig$values))
  
  A<-diag(1/sqrt(eig$values[1:q]))
  Q<-eig$vectors[,1:q]
  C<-Q%*%A%*%t(Q)
  return(list(mC=C,Q=Q,A=A,complete=TRUE))
}
##########################################################
IRMMA.aK.dK.effects.B<-function(Z=NULL,aK,dK,deltahat.aK,deltahat.dK,complete=TRUE){
  if( is.null(Z)){
    return(IRMMA.aK.dK.effects.B.Zo(aK,dK,deltahat.aK,deltahat.dK))
  }
  else{
    return(IRMMA.aK.dK.effects.B.Z(Z,aK,dK,deltahat.aK,deltahat.dK,complete))
  }
}
####
IRMMA.aK.dK.effects.B.Zo<- function(aK,dK,deltahat.aK,deltahat.dK){
  n<- nrow(aK)
  stopifnot(nrow(dK)==n)
  
  stopifnot(ncol(dK)==n)
  stopifnot(ncol(dK)==n)
  
  B <- deltahat.aK*aK+deltahat.dK*dK+diag(1,n)
  eig <- eigen(B,symmetric=TRUE,EISPACK=TRUE)
  qr.B<-qr(B)
  q<-qr.B$rank
  
  stopifnot(!is.complex(eig$values))
  
  A<-diag(1/sqrt(eig$values[1:q]))
  Q<-eig$vectors[,1:q]
  C<-Q%*%A%*%t(Q)
  return(list(C.g=C,Q=Q,A=A))
}

####
IRMMA.aK.dK.effects.B.Z<-function(Z,aK,dK,deltahat.aK,deltahat.dK,complete){
  if ( complete == FALSE ) {
    vids <- colSums(Z)>0
    Z <- Z[,vids]
    aK <- aK[vids,vids]
    dK <- dK[vids,vids]
  }
  
  n <- nrow(Z)  
  t <- nrow(aK)
  
  stopifnot(ncol(aK)==t)
  stopifnot(nrow(dK)==t)
  
  stopifnot(ncol(dK)==t)
  
  B <- Z%*%(deltahat.aK*aK+deltahat.dK*dK)%*%t(Z)+diag(1,n)
  eig <- eigen(B,symmetric=TRUE,EISPACK=TRUE)
  qr.B<-qr(B)
  q<-qr.B$rank
  
  stopifnot(!is.complex(eig$values))
  
  A<-diag(1/sqrt(eig$values[1:q]))
  Q<-eig$vectors[,1:q]
  C<-Q%*%A%*%t(Q)
  return(list(C.g=C,Q=Q,A=A,complete=TRUE))
}



################################################################################################
emma.MLE0.c <- function(Y_c,W_c){
  
  n <- length(Y_c)
  stopifnot(nrow(W_c)==n)
  M_c<-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
  etas<-crossprod(M_c,Y_c)
  LL <- 0.5*n*(log(n/(2*pi))-1-log(sum(etas*etas)))
  return(list(ML=LL))
}

emma.ML.LRT.c.noalpha <- function(ys, xs, K=1, Z, X0, ngrids=100, llim=-10, ulim=10, esp=1e-10) {
  #K=1,Z=C:n*n,X0=W:n*c(or,1_n*1),ys:n*1
  
  stopifnot(K == 1)
  
  
  ys <- Z%*%ys  
  xs <- Z%*%xs
  X0 <- Z%*%X0
  
  ys<-as.matrix(ys)
  xs<-as.matrix(xs)
  X0<-as.matrix(X0)
  n <- nrow(ys)
  m <- nrow(xs)
  t <- ncol(xs)
  q0 <- ncol(X0)

  MLE0<-emma.MLE0.c(ys,X0)

  ML1s <- vector(length=t)
  ML0s <- vector(length=t)
  vgs <- vector(length=t)
  ves <- vector(length=t)
  deltas<-vector(length=t)
  bhats<-vector(length=t)
  
  var.bhats.ratio <- vector(length=t)
  d <- vector(length=t)
  stats <- vector(length=t)
  ps <- vector(length=t)
  
  for (i in 1:t){
    
    vids <- !is.na(xs[,i])
    xv <- xs[vids,i]
    
    yv <- ys[vids]
    x0v<-X0[vids,]
    
    MLE1 <- emma.MLE.c (yv, x0v, K=1, xv, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
    if(length(MLE1$vg)!=0){
      ML1s[i]<-MLE1$ML
      ML0s[i]<-MLE0$ML
      vgs[i]<-MLE1$vg
      ves[i]<-MLE1$ve
      deltas[i]<-MLE1$delta
      
      
      nv<-length(MLE1$etas)
      Lam<-diag(c(1/(MLE1$delta*MLE1$lambda+1),rep(1,nv-1)))
      Lam1 <- diag(c(1/(MLE1$delta*MLE1$lambda+1)^2,rep(1,nv-1)))
      
      temp <- crossprod(xv,MLE1$U_R)
      bhats[i] <- MLE1$delta*temp%*%Lam%*%MLE1$etas
      var.bhats.ratio[i] <- MLE1$delta^2*temp%*%Lam%*%t(temp)%*%temp%*%Lam%*%t(temp)+MLE1$delta*temp%*%Lam1%*%t(temp)
      
      d[i] <- (1-var.bhats.ratio[i])*((1-var.bhats.ratio[i])>=0) #to record me=sum(d)
      stats[i]<- 2*(MLE1$ML-MLE0$ML)
      ps[i]<-if(stats[i]<=1e-100) 1 else pchisq(stats[i],1,lower.tail=F)/2#20160619
      
    }else{ps[i]<-1}
  }
  return(list(ID=1:t,ps=ps,bhats=bhats,deltas=deltas,d=d,ML1s=ML1s,ML0s=ML0s,stats=stats,vgs=vgs,ves=ves))
  
} 
######################################
emma.REMLE0.c <- function(Y_c,W_c){
  
  n <- length(Y_c)
  stopifnot(nrow(W_c)==n)
  M_c <-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
  eig <-eigen(M_c)
  t <-qr(W_c)$rank
  v <-n-t
  U_R <-eig$vector[,1:v]
  etas<-crossprod(U_R,Y_c)
  LL <- 0.5*v*(log(v/(2*pi))-1-log(sum(etas*etas)))
  return(list(REML=LL))
  
}

emma.REML.LRT.c.noalpha <- function(ys, xs, K=1, Z, X0, ngrids=100, llim=-10, ulim=10, esp=1e-10) {
  #K=1,Z=C:n*n,X0=W:n*c(or,1_n*1),ys:n*1
  
  stopifnot(K == 1)
 
  ys <- Z%*%ys   
  xs <- Z%*%xs
  X0 <- Z%*%X0
  
  ys<-as.matrix(ys)
  xs<-as.matrix(xs)
  X0<-as.matrix(X0)
  
  
  n <- nrow(ys)
  m <- nrow(xs)
  t <- ncol(xs)
  q0 <- ncol(X0)

  MLE0<-emma.REMLE0.c(ys,X0)

  ML1s <- vector(length=t)
  ML0s <- vector(length=t)
  vgs <- vector(length=t)
  ves <- vector(length=t)
  deltas <- vector(length=t)
  bhats<-vector(length=t)
  
  var.bhats.ratio <- vector(length=t)
  d <- vector(length=t)
  
  stats <- vector(length=t)
  ps <- vector(length=t)
  
  for (i in 1:t){
    
    vids <- !is.na(xs[,i])
    xv <- xs[vids,i]
    
    yv <- ys[vids]
    x0v<-X0[vids,]
    
    MLE1 <- emma.REMLE.c (yv, x0v, K=1, xv, ngrids=100, llim=-10, ulim=10,esp=1e-10, eig.L = NULL, eig.R = NULL)
    if(length(MLE1$vg)!=0){
      ML1s[i]<-MLE1$REML
      ML0s[i]<-MLE0$REML
      vgs[i]<-MLE1$vg
      ves[i]<-MLE1$ve
      deltas[i] <- MLE1$delta
      
      nv<-length(MLE1$etas)
      Lam<-diag(c(1/(MLE1$delta*MLE1$lambda+1),rep(1,nv-1)))
      Lam1 <- diag(c(1/(MLE1$delta*MLE1$lambda+1)^2,rep(1,nv-1)))
      temp <- crossprod(xv,MLE1$U_R)
      bhats[i] <- MLE1$delta*temp%*%Lam%*%MLE1$etas
      var.bhats.ratio[i] <- MLE1$delta^2*temp%*%Lam%*%t(temp)%*%temp%*%Lam%*%t(temp)+MLE1$delta*temp%*%Lam1%*%t(temp)
      
      d[i] <- (1-var.bhats.ratio[i])*((1-var.bhats.ratio[i])>=0) #to record me=sum(d)
      stats[i]<- 2*(MLE1$REML-MLE0$REML)
      ps[i]<-if(stats[i]<=1e-100) 1 else pchisq(stats[i],1,lower.tail=F)/2
      
    }else{
      ps[i]<-1
    }
  }
   return(list(ID=1:t,ps=ps,bhats=bhats,deltas=deltas,d=d,ML1s=ML1s,ML0s=ML0s,stats=stats,vbs=vgs,ves=ves))
  
}
######################################

fixed.REMLE0.c <- function(Y_c,W_c){
  
  n <- length(Y_c)
  stopifnot(nrow(W_c)==n)
  M_c <-diag(1,n)-W_c%*%solve(crossprod(W_c,W_c))%*%t(W_c)
  t <-qr(W_c)$rank
  v <-n-t
  etas<-crossprod(M_c,Y_c)

  LL <- 0.5*v*(log(v/(2*pi))-1-log(sum(etas*etas)))
  return(list(REML=LL))
  
}

fixed.REML.LRT.c.sim<-function(ys, xs, Z, X0){ 
  
  
  ys <- Z%*%ys          
  xs <- Z%*%xs
  X0 <- Z%*%X0
  
  ys<-as.matrix(ys)
  xs<-as.matrix(xs)
  X0<-as.matrix(X0)
  
  
  n <- nrow(ys)
  m <- nrow(xs)
  t <- ncol(xs)
  stats <- vector(length=t)
  ps <- vector(length=t)
  
  
  ML1s <- vector(length=t)
  ML0s <- vector(length=t)
  REML0.c<-fixed.REMLE0.c(ys,X0)
  
  for(i in 1:t){
    vids <- !is.na(xs[,i])
    xv <- xs[vids,i]
    
    yv <- ys[vids]
    x0v<-X0[vids,]
    
    xv.new<-cbind(x0v,xv)
    REML1.c<-fixed.REMLE0.c(yv,xv.new)
    
    stats[i]<- 2*(REML1.c$REML-REML0.c$REML)
    ps[i]<-pchisq(stats[i],1,lower.tail=F)
    
    ML1s[i]<-REML1.c$REML
    ML0s[i]<-REML0.c$REML
    
    
  }
  return(list(ID=1:t,ps=ps,ML1s=ML1s,ML0s=ML0s,stats=stats))

} 
##########################

  
fixed.ML.LRT.c.sim<-function(ys, xs, Z, X0){
 
  ys <- Z%*%ys              
  xs <- Z%*%xs
  X0 <- Z%*%X0
  
  ys<-as.matrix(ys)
  xs<-as.matrix(xs)
  X0<-as.matrix(X0)
  
  
  n <- nrow(ys)
  m <- nrow(xs)
  t <- ncol(xs)

  stats <- vector(length=t)
  ps <- vector(length=t)
  
  
  ML1s <- vector(length=t)
  ML0s <- vector(length=t)
  
  REML0.c<-emma.MLE0.c(ys,X0)#ys=Y_c,X0=W_c
  
  for(i in 1:t){
    vids <- !is.na(xs[,i])
    xv <- xs[vids,i]
    
    yv <- ys[vids]
    x0v<-X0[vids,]
    
    xv.new<-cbind(x0v,xv)
    
    REML1.c<-emma.MLE0.c(yv,xv.new)
    
    stats[i]<- 2*(REML1.c$ML-REML0.c$ML)
    ps[i]<-pchisq(stats[i],1,lower.tail=F)
    
    ML1s[i]<-REML1.c$ML
    ML0s[i]<-REML0.c$ML
    
    
  }
  return(list(ID=1:t,ps=ps,ML1s=ML1s,ML0s=ML0s,stats=stats))
  
} 


######################################
replaceNaN<-  function(LL) {
  #handler of grids with NaN log 
  index=(LL=="NaN")
  if(length(index)>0) theMin=min(LL[!index])
  if(length(index)<1) theMin="NaN"
  LL[index]=theMin
  return(LL)    
}

#########################################################################
peak.id<-function(Lod.temp){
  m<-length(Lod.temp)
  optids<-vector(length=0)
  if(Lod.temp[1]>Lod.temp[2])   optids<-append(optids,1)
  
  for(j in 2:(m-1)){
    if ((Lod.temp[j-1]<Lod.temp[j]) & (Lod.temp[j]>Lod.temp[j+1])) {
      optids<-append(optids,j)
    }
  }
  if(Lod.temp[m]>Lod.temp[m-1])  optids<-append(optids,m)
  return(optids)
}
##############################
#Multinormal distribution density function
multinormal<-function(y,mean,sigma)
{
  pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
  return (pdf_value)
}
#LOD value test #library(MASS)
######################################3
likelihood.a.d.F2<-function(xxn,xxx,yn,bbo,intercept)
  #xxn:fix matrix;xxx:gene matrix;yn:pheno matrix;bbo:gene effect from adalasso
{
  nq<-ncol(xxx)
  ns<-nrow(yn)
  at1<-0
  ww1<-as.matrix(which(abs(bbo)>1e-5))

  ww.a<-ww1[ww1%%2==1]
  ww.d<-ww1[ww1%%2==0]
  
  ww.a.new<-c(ww.a,ww.a+1)
  ww.d.new<-c(ww.d,ww.d-1)
  
  ww1.new<-union(ww.a.new,ww.d.new)
  
  ww1.new<-ww1.new[order(ww1.new)]
  ww1.new<-as.matrix(ww1.new)

  at1<-dim(ww1.new)[1]
  lod<-matrix(rep(0,nq),nq,1)
  ps<-matrix(rep(1,nq),nq,1)
  
  ad<-if(at1>0.5) cbind(xxn,xxx[,ww1.new]) else xxn
  
  bb<-if (is.null(intercept)) ginv(crossprod(ad,ad))%*%crossprod(ad,yn) else c(intercept,bbo[ww1.new])

  vv1<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns)
  ll1<-sum(log(abs(multinormal(yn,ad%*%bb,vv1))))
  
  sub<-1:ncol(ad)

  at2<-if(at1>1) seq(1,at1,by=2) else 1
  if(at1>0.5)
  {
    for(i in at2)
    {
      ij<-which((sub!=sub[i+ncol(xxn)])&(sub!=sub[i+ncol(xxn)+1]))

      ad1<-ad[,ij,drop=F]

      bb1<-if (is.null(intercept)) ginv(crossprod(ad1,ad1))%*%crossprod(ad1,yn) else c(intercept,bbo[ww1.new])[ij]
 
      vv0<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns)
      ll0<-sum(log(abs(multinormal(yn,ad1%*%bb1,vv0))))
      lod[ww1.new[i]]<--2.0*(ll0-ll1)/(2.0*log(10))
      
      ps[ww1.new[i]]<-pchisq(-2.0*(ll0-ll1),2,lower.tail = F)
    }
  }
  return (list(lod=lod,ps=ps))
}

