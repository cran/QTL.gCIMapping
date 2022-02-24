#' To perform QTL mapping with Wen method
#'
#' @param pheRaw phenotype matrix.
#' @param genRaw genotype matrix.
#' @param mapRaw1 linkage map matrix.
#' @param WalkSpeed Walk speed for Genome-wide Scanning.
#' @param CriLOD Critical LOD scores for significant QTL.
#' @param dir file path in your computer.
#'
#' @return a list
#' @export
#'
#' @examples
#' data(F2data)
#' readraw<-Readdata(file=F2data,fileFormat="GCIM",
#' method="GCIM-QEI",filecov=NULL,
#' MCIMmap=NULL,MultiEnv=TRUE)
#' DoResult<-Dodata(fileFormat="GCIM",
#' Population="F2",method="GCIM-QEI",
#' Model="Random",readraw,MultiEnv=TRUE)
#' ZhouMatrices<-ZhouF(pheRaw=DoResult$pheRaw,
#' genRaw=DoResult$genRaw,
#' mapRaw1=DoResult$mapRaw1,
#' WalkSpeed=1,CriLOD=3,
#' dir=tempdir())
ZhouF<-function(pheRaw=NULL,genRaw=NULL,mapRaw1=NULL,WalkSpeed=NULL,CriLOD=NULL,dir=NULL){
  cl<-WalkSpeed
  sLOD<-CriLOD
  # yygg<-NULL
  # mx=NULL;phe=NULL;chr_name=NULL;v.map=NULL
  if(is.null(genRaw)==TRUE){
    warning("Please input correct genotype dataset!")
  }
  if(is.null(pheRaw)==TRUE){
    warning("Please input correct phenotype dataset!")
  }
  if(is.null(mapRaw1)==TRUE){
    warning("Please input correct linkage map dataset!")
  }
  if((is.null(genRaw)==FALSE)&&(is.null(pheRaw)==FALSE)&&(is.null(mapRaw1)==FALSE)&&(cl<0)){
    warning("Please input Walk Speed: >0!")
  }
  if((is.null(genRaw)==FALSE)&&(is.null(pheRaw)==FALSE)&&(is.null(mapRaw1)==FALSE)&&(cl>0)&&(sLOD<0)){
    warning("Please input critical LOD score: >0!")
  }

  mapRaw<-as.matrix(mapRaw1)
  chr_name<-unique(mapRaw[,2])
  chr_secon<-as.matrix(mapRaw[,2])

  mm<-numeric()
  map_chr<-numeric()
  for(i in 1:length(chr_name)){
    chr_i<-length(which(chr_secon[]==chr_name[i]))
    mm<-c(mm,chr_i)
    chr_name[i]<-i
    map_chr<-c(map_chr,rep(i,chr_i))
  }
  mm<-matrix(mm,ncol=1)
  map_chr<-matrix(map_chr,ncol=1)

  mapRaw[,2]<-map_chr

  chr<-length(chr_name)

  for(i in 1:chr){
    pos1<-as.matrix(mapRaw[which(mapRaw[,2]==i),3])
    delerow<-which(duplicated(pos1))
    if(length(delerow)!=0){
      break
    }
  }
  if(length(delerow)!=0){
    warning("Please check linkage maps (linkage groups) to make sure whether all the marker positions are different!")
  }else{

    blank<-matrix("",nrow=3,ncol=dim(pheRaw)[2])
    blank[1,]<-colnames(pheRaw)
    # p1<-as.matrix(pheRaw)
    p1<-pheRaw
    colnames(p1)<-NULL
    p1<-t(rbind(blank,p1))
    g1<-cbind(mapRaw,genRaw)
    colnames(p1)<-NULL
    colnames(g1)<-NULL

    pgcombine<-rbind(p1,g1)
    write.table(pgcombine,file=paste(dir,"/listeria_rotY",".csv",sep=""),sep=",",row.names = F,col.names = F)

    ##########  calculate conditional probability for K matrix
    f2<-read.cross("csvr",dir,"listeria_rotY.csv",genotypes=c("A","H","B","D","C"),na.strings = "-",crosstype="f2")
    # f2<-jittermap(f2)# Jitter the marker positions in a genetic map so that no two markers are on top of each other   # jittermap(object, amount=1e-6)
    simf2<-calc.genoprob(f2, step=0,error.prob = 0.0001)

    ##########  Access to chromosome information
    genoname<-apply(mapRaw[,2:3],2,as.numeric)
    genoname<-data.frame(marker=mapRaw[,1],chr=genoname[,1],pos=genoname[,2])
    chr_n<-as.numeric(genoname[,2])
    chr_n<-chr_n[!duplicated(chr_n)]
    maxdistance<-0
    for(i in 1:length(chr_n)){
      maxdistance<-max(maxdistance,   max(diff(as.matrix(genoname[which(genoname[,2]==i),3]))))
    }

    ##########  calculate conditional probability for K matrix
    Ax0<-NULL;Hx0<-NULL;Bx0<-NULL
    for(i in 1:length(chr_n)){
      map_gen<-simf2$geno[[i]]$prob
      A_gen<-round(map_gen[,,1],digits=15)
      H_gen<-round(map_gen[,,2],digits=15)
      B_gen<-round(map_gen[,,3],digits=15)
      Ax0<-cbind(Ax0,A_gen)
      Hx0<-cbind(Hx0,H_gen)
      Bx0<-cbind(Bx0,B_gen)
    }# dim(Ax0) # mn<-dim(Ax)[2]

    ##########  Whether need to insert markers
    if(maxdistance>cl){#user's options
      simf2<-calc.genoprob(f2, step=cl,error.prob = 0.0001)
      Ax<-NULL;Hx<-NULL;Bx<-NULL;regenoname<-NULL
      for(i in 1:length(chr_n)){
        map_gen<-simf2$geno[[i]]$prob
        A_gen<-round(map_gen[,,1],digits=15)
        H_gen<-round(map_gen[,,2],digits=15)
        B_gen<-round(map_gen[,,3],digits=15)
        Ax<-cbind(Ax,A_gen)
        Hx<-cbind(Hx,H_gen)
        Bx<-cbind(Bx,B_gen)
        nowpos<-attr(map_gen,"map")
        nowbin<-names(nowpos)
        nowchr<-rep(as.numeric(chr_n[i]),length(nowpos))
        nowdata<-data.frame(marker=nowbin,chr=nowchr,pos=nowpos)
        regenoname<-rbind(regenoname,nowdata)
      }# dim(Ax)
      rownames(regenoname)<-NULL
      regenoname<-cbind(regenoname,seq(1,dim(regenoname)[1],1))
      colnames(regenoname)<-c("marker","chr","pos","id.all")
      # mn<-dim(regenoname)[1]
      genoname<-regenoname
    }else{
      Ax<-Ax0;Hx<-Hx0;Bx<-Bx0
      genoname<-cbind(genoname,seq(1,nrow(genoname),by=1))
      colnames(genoname)<-c("marker","chr","pos","id.all")
    }
  }
  output<-list(genoname=genoname,mapRaw=mapRaw,
               Ax0=Ax0,Hx0=Hx0,Bx0=Bx0,Ax=Ax,Hx=Hx,Bx=Bx)# yygg=yygg,pheRaw=pheRaw,chr_n=chr_n,
  return(output)
}

#' The second step of Zhou method for single environment
#'
#' @param Model Random or fixed model.
#' @param pheRaw phenotype matrix.
#' @param genRaw genotype matrix.
#' @param mapRaw linkage map matrix.
#' @param CriLOD Critical LOD scores for significant QTL.
#' @param NUM The serial number of the trait to be analyzed.
#' @param yygg covariate matrix.
#' @param genoname linkage map matrix with pseudo markers inserted.
#' @param Ax0 AA genotype matrix.
#' @param Hx0 Aa genotype matrix.
#' @param Bx0 aa genotype matrix.
#' @param Ax AA genotype matrix with pseudo markers inserted.
#' @param Hx Aa genotype matrix with pseudo markers inserted.
#' @param Bx aa genotype matrix with pseudo markers inserted.
#' @param dir file storage path.
#' @param CriDis The distance of optimization.
#' @param CLO Number of CPUs.
#'
#' @return a list
#' @export
#'
#' @examples
#' data(F2data)
#' readraw<-Readdata(file=F2data,fileFormat="GCIM",
#' method="GCIM-QEI",filecov=NULL,
#' MCIMmap=NULL,MultiEnv=FALSE)
#' DoResult<-Dodata(fileFormat="GCIM",Population="F2",
#' method="GCIM-QEI",Model="Random",
#' readraw,MultiEnv=FALSE)
#' ZhouMatrices<-ZhouF(pheRaw=DoResult$pheRaw,
#' genRaw=DoResult$genRaw,mapRaw1=DoResult$mapRaw1,
#' WalkSpeed=1,CriLOD=3,dir=tempdir())
#' OutputZhou<-ZhouMethod_single_env(Model="Random",
#' pheRaw=DoResult$pheRaw,genRaw=DoResult$genRaw,
#' mapRaw=ZhouMatrices$mapRaw,CriLOD=3,NUM=1,
#' yygg=DoResult$yygg1,genoname=ZhouMatrices$genoname,
#' Ax0=ZhouMatrices$Ax0,Hx0=ZhouMatrices$Hx0,
#' Bx0=ZhouMatrices$Bx0,Ax=ZhouMatrices$Ax,
#' Hx=ZhouMatrices$Hx,Bx=ZhouMatrices$Bx,
#' dir=tempdir(),CriDis=5,CLO=2)
ZhouMethod_single_env<-function(Model=NULL,pheRaw=NULL,genRaw=NULL,mapRaw=NULL,CriLOD=NULL,NUM=NULL,yygg=NULL,genoname=NULL,
                                Ax0=NULL,Hx0=NULL,Bx0=NULL,Ax=NULL,Hx=NULL,Bx=NULL,dir=NULL,CriDis=NULL,CLO=NULL){# chr_n=NULL,

  #########################################  function  #########################################
  kinship_every<-function(coded_gen){
    kk<-coded_gen%*%t(coded_gen)
    kk<-kk/mn
    return(kk)
  }
  p3d_method<-function(x,y,kinship){
    # estimate the value of λ & λk        kinship=K;
    P3D<-function(x,y,kinship){
      # iteration function(H=K*λ+I);estimate λ
      iter_p3d<-function(ga){
        lambda<-exp(ga)
        diag_element<-lambda*K_value+1
        logH<-sum(log(diag_element))
        RH_value<-1/(diag_element)
        yuRHyu<-sum(yu*RH_value*yu)
        yuRHxu<-matrix(0,nrow = 1,ncol = q)
        xuRHxu<-matrix(0,nrow = q,ncol = q)
        for(i in 1:q){
          yuRHxu[,i]<-sum(yu*RH_value*xu[,i])
          for(j in 1:q){
            xuRHxu[i,j]<-sum(xu[,i]*RH_value*xu[,j])
          }
        }
        logxuRHxu<-log(det(xuRHxu))
        logyuPyu<-log(yuRHyu-yuRHxu%*%tcrossprod(solve(xuRHxu),yuRHxu))
        output<- -0.5*(logH+logxuRHxu+(n-q)*logyuPyu)
        return(-output)
      }

      q<-ncol(x)
      n<-nrow(y)
      eigenK<-eigen(kinship)
      K_vector<-eigenK$vectors
      K_value<-eigenK$values
      # rm(eigenK);gc()
      xu<-crossprod(K_vector,x)
      yu<-crossprod(K_vector,y)

      ga0<-0
      optimp3d<-optim(par=ga0,fn=iter_p3d,hessian = TRUE,method="L-BFGS-B",lower=-50,upper=10)
      lambda<-exp(optimp3d$par)

      return(list(lambda,K_vector,K_value))
    }

    q<-ncol(x)
    value1<-P3D(x=x,y=y,kinship=kinship)
    lambda<-value1[[1]]
    uu<-as.matrix(value1[[2]])#  The eigenvector of K matrix
    vv<-value1[[3]] # The eigenvalue of K matrix
    # RH_value<-1/(vv*lambda+1) # rm(value1);gc()
    return(list(lambda=lambda,k_vector=uu,k_value=vv))
  }
  single_locus_model<-function(x,y,zz,lambda,uu,vv,CLO){
    value2<-matrix(c(0.5,-0.5,0,1,-0.5,-0.5),2,3)
    # iteration function(R=zu%*%t(zu)*λk+D*λ+I);estimate λk
    rqtl<-function(ga){
      lambdak<-exp(ga)
      Hk_term<-zuRHzu*lambdak+diag(1,3)
      logHk<-sum(log(lambda*vv+1))+log(det(Hk_term))
      RHk_term<-solve(Hk_term)*lambdak
      yuRHkyu<-yuRHyu-yuRHzu%*%tcrossprod(RHk_term,yuRHzu)
      yuRHkxu<-yuRHxu-yuRHzu%*%tcrossprod(RHk_term,xuRHzu)
      xuRHkxu<-xuRHxu-xuRHzu%*%tcrossprod(RHk_term,xuRHzu)
      yuPkyu<-yuRHkyu-yuRHkxu%*%tcrossprod(solve(xuRHkxu),yuRHkxu)
      rqtl<- -0.5*( logHk + log(det(xuRHkxu)) + (n-q)*log(yuPkyu) )
      return(-rqtl)
    }

    #estimate the value of γ
    gamma_estimate<-function(lambdak){
      Hk_term<-zuRHzu*lambdak+diag(1,3)
      RHk_term<-solve(Hk_term)*lambdak
      yuRHkyu<-yuRHyu-yuRHzu%*%tcrossprod(RHk_term,yuRHzu)
      yuRHkxu<-yuRHxu-yuRHzu%*%tcrossprod(RHk_term,xuRHzu)
      xuRHkxu<-xuRHxu-xuRHzu%*%tcrossprod(RHk_term,xuRHzu)
      zuRHkxu<-t(xuRHzu)-zuRHzu%*%tcrossprod(RHk_term,xuRHzu)
      zuRHkyu<-t(yuRHzu)-zuRHzu%*%tcrossprod(RHk_term,yuRHzu)
      zuRHkzu<-zuRHzu-zuRHzu%*%RHk_term%*%zuRHzu
      beta<-solve(xuRHkxu,t(yuRHkxu))
      beta<-matrix(beta,ncol = 1)
      yuPkyu<-yuRHkyu-yuRHkxu%*%tcrossprod(solve(xuRHkxu),yuRHkxu)
      sigma<-yuPkyu/(n-q)
      sigma<-as.numeric(sigma)
      gamma<-lambdak*zuRHkyu-lambdak*zuRHkxu%*%beta
      var<-abs((lambdak*diag(1,3)-lambdak*zuRHkzu*lambdak)*sigma)
      stderr<-sqrt(diag(var))
      phi.k<-sigma*lambdak#Φk
      return(list(gamma,beta,var,phi.k,sigma,stderr))
    }

    # estimate the value of logp
    logp_estimate<-function(L, g_k, pn){
      var.1<-L%*%tcrossprod(var_ga,L)
      Wk.1<-crossprod(g_k,ginv(var.1))%*%g_k
      rank1<-qr(L)$rank

      tr1<-sum(diag(ginv(tcrossprod(L,L))%*%var.1))
      dk1<-rank1-tr1/phi_k

      p1<-pgamma(Wk.1,shape=pn/2,scale=2*dk1,lower.tail = FALSE,log.p = FALSE)
      log1<-pgamma(Wk.1,shape=pn/2,scale=2*dk1,lower.tail = FALSE,log.p = TRUE)
      log1<--log1*log10(exp(1))
      return(list(p=p1,log=log1))
    }

    RH_value<-1/(vv*lambda+1)
    mn<-ncol(zz)/3
    n<-nrow(y)
    q<-ncol(x)

    xu<-crossprod(uu,x)
    yu<-crossprod(uu,y)

    yuRHyu<-sum(yu*RH_value*yu)
    yuRHxu<-matrix(0,nrow=1,ncol=q)
    xuRHxu<-matrix(0,nrow=q,ncol=q)
    for(i in 1:q){
      yuRHxu[,i]<-sum(yu*RH_value*xu[,i])
      for(j in 1:q){
        xuRHxu[j,i]<-sum(xu[,j]*RH_value*xu[,i])
      }
    }
    # yuRHyu<-crossprod(yu,diag(RH_value))%*%yu
    # yuRHxu<-crossprod(yu,diag(RH_value))%*%xu
    # xuRHxu<-crossprod(xu,diag(RH_value))%*%xu
    if(is.null(CLO)==TRUE){
      cl.cores <- detectCores()
    if(cl.cores<=2){
      cl.cores<-1
    }else{
      if(cl.cores>10){
        cl.cores <-10
      }else{
        cl.cores <- detectCores()-1
      }
    }
    # cl.cores<-2
    }else{
      cl.cores <-CLO
    }

    cl <- makeCluster(cl.cores)
    registerDoParallel(cl)

    SR_i<-numeric()
    result<-foreach(SR_i=1:mn,.combine=rbind)%dopar%{
      # library(MASS)
      z<-zz[,((SR_i-1)*3+1):(SR_i*3),drop=F]
      zu<-crossprod(uu,z)

      xuRHzu<-matrix(0,nrow=q,ncol=3)
      yuRHzu<-matrix(0,nrow=1,ncol=3)
      zuRHzu<-matrix(0,nrow=3,ncol=3)
      for(i in 1:3){
        yuRHzu[,i]<-sum(yu*RH_value*zu[,i])
        for(j in 1:q){
          xuRHzu[j,i]<-sum(xu[,j]*RH_value*zu[,i])
        }
        for(j in 1:3){
          zuRHzu[j,i]<-sum(zu[,j]*RH_value*zu[,i])
        }
      }
      # xuRHzu<-crossprod(xu,diag(RH_value))%*%zu
      # yuRHzu<-crossprod(yu,diag(RH_value))%*%zu
      # zuRHzu<-crossprod(zu,diag(RH_value))%*%zu

      ga<-0
      par<-optim(par=ga,fn=rqtl,hessian = TRUE,method="L-BFGS-B",lower=-10,upper=10)
      lambdak<-exp(par$par)

      value3<-gamma_estimate(lambdak)
      gamma_k<-value3[[1]]
      beta<-value3[[2]]#estimate β   (y=Xβ+Zγ+ξ+ε)
      var_ga<-value3[[3]]
      phi_k<-value3[[4]]
      sigma_2<-value3[[5]]

      main_effect<-value2%*%gamma_k

      logvalue1<-logp_estimate(L=value2, g_k=main_effect, pn=2)#the -log(10)p of qtl effect
      p1<-logvalue1$p
      log1<-logvalue1$log

      result<-cbind(lambdak,beta,matrix(gamma_k,1,3),p1,log1,phi_k,sigma_2)
    }
    stopCluster(cl)

    lambda_k<-result[,1]
    mu_beta<-result[,2]
    gamma_all<-result[,3:5]
    p1<-result[,6]
    log_p1<-result[,7]
    phi_k<-result[,8]
    sigma_2<-result[,9]

    return(list(lambda_k=lambda_k, fixed=mu_beta, gamma=gamma_all,
                p1=p1, log1=log_p1, phi_k=phi_k, sigma_2=sigma_2))

  }
  single_locus_model_Fixed<-function(x,y,zz,lambda,uu,vv,CLO){
    value2<-matrix(c(0.5,-0.5,0,1,-0.5,-0.5),2,3)
    # estimate the value of logp
    logp_estimate<-function(L, g_k, pn){
      var.1<-L%*%tcrossprod(var_ga,L)
      Wk.1<-crossprod(g_k,ginv(var.1))%*%g_k# Test statistic
      p1<-pchisq(Wk.1,df=pn,lower.tail = FALSE,log.p = FALSE)
      log1<-pchisq(Wk.1,df=pn,lower.tail = FALSE,log.p = TRUE)
      log1<--log1*log10(exp(1))
      return(list(p=p1,log=log1))
    }

    Hsolve<-1/(vv*lambda+1)

    mn<-ncol(zz)/3
    n<-nrow(y)
    q<-ncol(x)

    # xu<-crossprod(uu,x)
    yu<-crossprod(uu,y)# Equation deformation

    if(is.null(CLO)==TRUE){
      cl.cores <- detectCores()
      if(cl.cores<=2){
        cl.cores<-1
      }else{
        if(cl.cores>10){
          cl.cores <-10
        }else{
          cl.cores <- detectCores()-1
        }
      }
    }else{
      cl.cores <-CLO
    }
    cl <- makeCluster(cl.cores)
    registerDoParallel(cl)

    SF_i<-numeric()
    result<-foreach(SF_i=1:mn,.combine=rbind)%dopar%{
      # library(MASS)
      z<-zz[,((SF_i-1)*3+1):(SF_i*3),drop=F]
      uxz<-crossprod(uu,cbind(x,z))
      x_gamma<-ginv(t(uxz)%*%diag(Hsolve)%*%uxz)%*%t(uxz)%*%diag(Hsolve)%*%yu
      q<-qr(uxz)$rank
      sig_e2<-as.numeric(t(yu-uxz%*%x_gamma)%*%diag(Hsolve)%*%(yu-uxz%*%x_gamma)/(dim(uxz)[1]-q))
      x_gamma_covmatr<-sig_e2*ginv(t(uxz)%*%diag(Hsolve)%*%uxz)
      gamma<-x_gamma[-c(1:dim(x)[2])]
      var_ga<-x_gamma_covmatr[-c(1:dim(x)[2]),-c(1:dim(x)[2]),drop=F]

      main_effect<-value2%*%gamma

      logvalue1<-logp_estimate(L=value2, g_k=main_effect, pn=2)#the -log(10)p of qtl effect
      p1<-logvalue1$p
      log1<-logvalue1$log

      result<-cbind(x_gamma[c(1:dim(x)[2])],matrix(gamma,1,3),p1,log1,sig_e2)
    }
    stopCluster(cl)

    mu_beta<-result[,1:dim(x)[2]]
    gamma_all<-result[,(dim(x)[2]+1):(dim(x)[2]+3)]
    p1<-result[,(dim(x)[2]+4)]
    log_p1<-result[,(dim(x)[2]+5)]
    sigma_2<-result[,(dim(x)[2]+6)]

    return(list(fixed=mu_beta, gamma=gamma_all, p1=p1, log1=log_p1, sigma_2=sigma_2))
  }
  peak_selection<-function(log_value,genoname){
    peak_pos<-function(Lod.temp){
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
    chr_all<-as.matrix(genoname[,2])
    chr_kind<-chr_all[!duplicated(chr_all)]
    id_pos<-NULL
    for(jjj in chr_kind){
      now_id<-which(chr_all%in%jjj)
      id_pos<-c(id_pos,now_id[peak_pos(log_value[now_id])])
    }
    return(sort(id_pos))
  }
  multi_peak_new<-function(gencoded,peak_id){
    enk<-3
    mut_peak_id<-NULL
    term1<-seq(1,3,1)
    for(i in 1:length(peak_id)){
      mut_peak_id<-c(mut_peak_id,(rep(peak_id[i],enk)-1)*enk+term1)
    }
    return(list(z=gencoded[,sort(mut_peak_id)],order=sort(rep(peak_id,enk))))
  }
  Zhou_lars<-function(peak,CodeMatrix,n){
    multi_value<-multi_peak_new(CodeMatrix,peak)
    DesignMatrix<-multi_value$z
    order0<-multi_value$order# length(order0); length(peak_id)
    if(length(peak)>=n){
      larstep<-length(order0)%/%(3*5)
      lar_result<-lars(x=DesignMatrix, y=y, type = "lar",trace = FALSE, normalize = TRUE, intercept = TRUE, eps = .Machine$double.eps, use.Gram=FALSE,max.steps = larstep)
      lar_result.0<-lar_result$beta[nrow(lar_result$beta),]
      lar_pos0<-order0[which(lar_result.0!=0)]
      lar_pos<-lar_pos0[!duplicated(lar_pos0)]# length(lar_pos)
      multi_value1<-multi_peak_new(CodeMatrix,lar_pos)
      DesignMatrix1<-multi_value1$z # coefficient matrix of selected peak loci
      order1<-multi_value1$order
    }else{
      lar_pos<-peak
      DesignMatrix1<-DesignMatrix
      order1<-order0
    }# length(lar_pos)
    return(list(lar_pos=lar_pos,Matrix=DesignMatrix1,order=order1))
  }
  sblgwas<-function(x,y,z,t,max.iter=200,min.err=1e-6){
    x<-as.matrix(x)
    y<-as.matrix(y)
    z<-as.matrix(z)
    n<-length(y)
    q<-ncol(x)
    m<-ncol(z)
    b0<-solve(t(x)%*%x,tol=1e-50)%*%(t(x)%*%y)
    s2<-sum((y-x%*%b0)^2)/(n-q)
    b0<-matrix(0,q,1)
    b<-b0
    g0<-matrix(0,m,1)
    g<-g0
    lambda<-matrix(0,m,1)
    tau<-g0
    v<-g0
    xx<-NULL
    xy<-NULL
    for(i in 1:q){
      xx<-c(xx,sum(x[,i]^2))
      xy<-c(xy,sum(x[,i]*y))
    }
    zz<-NULL
    zy<-NULL
    for(k in 1:m){
      zz<-c(zz,sum(z[,k]^2))
      zy<-c(zy,sum(z[,k]*y))
    }
    d<-numeric(m)
    a<-matrix(0,n,1)
    iter<-0
    err<-1e8
    my.iter<-NULL
    while(iter < max.iter & err > min.err){
      for(i in 1:q){
        a<-a-x[,i]*b0[i]
        ai<-sum(x[,i]*a)
        b[i]<-(xy[i]-ai)/xx[i]
        a<-a+x[,i]*b[i]
      }
      df<-0
      for(k in 1:m){
        a<-a-z[,k]*g0[k]
        ak<-sum(z[,k]*a)
        c1<- -(t+3)*zz[k]^2
        c2<- -(2*t+5)*zz[k]+(zy[k]-ak)^2
        c3<- -(t+2)
        if( ((c2^2-4*c1*c3) < 0) | (c2 < 0) ){
          tau[k]<-0
        } else {
          tau[k]<-(-c2-sqrt(c2^2-4*c1*c3))/(2*c1)
        }
        lambda[k]<-tau[k]/s2
        g[k]<-lambda[k]*(zy[k]-ak)-lambda[k]^2*zz[k]*(zy[k]-ak)/(lambda[k]*zz[k]+1)
        d[k]<-lambda[k]*(zz[k]-lambda[k]*zz[k]^2/(lambda[k]*zz[k]+1))
        v[k]<-tau[k]-tau[k]*d[k]
        df<-df+d[k]
        a<-a+z[,k]*g[k]
      }

      if((n-q-df) > 0){s2<-sum((y-a)^2)/(n-q-df)
      }else{
        s2<-sum((y-a)^2)/(n-q)
      }

      iter<-iter+1
      err<-sum((g-g0)^2)/m
      g0<-g
      b0<-b
      my.iter<-rbind(my.iter,cbind(iter,err,s2,t(b),t(g)))
    }
    my.parm<-data.frame(iter,err,s2,b,df)
    names(my.parm)<-c("iter","error","s2","beta","df")

    posv<-which(v!=0)
    m<-length(g)
    wald<-c(rep(0,m))
    gg<-g[posv]
    vv<-v[posv]
    wald[posv]<-gg^2/vv
    p<-pchisq(wald,1,lower.tail=FALSE)

    my.blup<-data.frame(g,v,wald,p)
    names(my.blup)<-c("gamma","vg","wald","p_wald")

    var.beta<-NULL
    for(i in 1:q){
      var.beta<-c(var.beta,paste("beta",i,sep=""))
    }
    var.gamma<-NULL
    for(k in 1:m){
      var.gamma<-c(var.gamma,paste("gamma",k,sep=""))
    }
    var.names<-c(c("iter","error","s2"),var.beta,var.gamma)
    my.iter<-data.frame(my.iter)
    names(my.iter)<-var.names

    out<-list(my.iter,my.parm,my.blup)
    names(out)<-c("iteration","parm","blup")
    return(out)
  }
  selection<-function(posx,genoname,svrad){
    chose_peak<-c(posx[1])
    order_now<-1
    while(order_now<length(posx)){
      order_now<-order_now+1
      repeat_pos<-which( abs(chose_peak-as.numeric(posx[order_now]))<=(svrad)  )
      if(length(repeat_pos)>0){
        if_condition<-length(which(  genoname[chose_peak[repeat_pos],2]==as.numeric(genoname[posx[order_now],2])  ))==0
        if(if_condition){
          chose_peak<-c(chose_peak,posx[order_now])
        }
      }else{
        chose_peak<-c(chose_peak,posx[order_now])
      }
    }
    return(chose_peak)
  }
  selection2<-function(posx1,posx2,genoname,svrad){
    chose_peak<-NULL
    order_now<-0
    while(order_now<length(posx1)){
      order_now<-order_now+1
      repeat_pos<-which( abs(posx2-as.numeric(posx1[order_now]))<=(svrad)  )
      if(length(repeat_pos)>0){
        if_condition<-length(which(  genoname[posx2[repeat_pos],2]==as.numeric(genoname[posx1[order_now],2])  ))==0
        if(if_condition){
          chose_peak<-c(chose_peak,posx1[order_now])
        }
      }else{
        chose_peak<-c(chose_peak,posx1[order_now])
      }
    }
    return(chose_peak)
  }
  ebayes_EM<-function(x,z,y,v0,v,tau,err_max){
    n<-nrow(z);k<-ncol(z)
    mk<-3; kn<-k/mk
    v0<-as.numeric(v0)
    v<-matrix(v,ncol=1)
    if(abs(min(eigen(crossprod(x,x))$values))<1e-6){
      try_b<-try({  b<-chol2inv(chol(crossprod(x,x)+diag(ncol(x))*1e-8))%*%crossprod(x,y)  },silent=TRUE)
      if('try-error' %in% class(try_b)){
        try_c<-try({  b<-solve(crossprod(x,x))%*%crossprod(x,y)   },silent=TRUE)
        if('try-error' %in% class(try_c)){   b<-ginv(crossprod(x,x))%*%crossprod(x,y)  }
      }
    }else{
      try_b<-try({  b<-chol2inv(chol(crossprod(x,x)))%*%(crossprod(x,y))  },silent=TRUE)
      if('try-error' %in% class(try_b)){
        try_c<-try({  b<-solve(crossprod(x,x))%*%(crossprod(x,y))   },silent=TRUE)
        if('try-error' %in% class(try_c)){  b<-ginv(crossprod(x,x))%*%(crossprod(x,y))   }
      }
    }
    u<-matrix(0,nrow=mk,ncol=kn)# E(γk)
    w<-matrix(0,nrow=mk,ncol=k)# var(γk)
    s<-matrix(0,nrow=kn,ncol=1)# tr(var(γk))
    vv<-matrix(0,n,n)
    for(i in 1:kn){
      nc<-( (i-1)*mk+1 ):(i*mk)
      zz<-z[,nc]# Zk
      vv=vv+tcrossprod(zz,zz)*v[i,]
    }
    vv<-vv+diag(n)*v0 # V

    L<-matrix(c(0.5,-0.5,0,1,-0.5,-0.5),2,3)
    rank_1<-qr(L)$rank

    iter<-0;err<-1000;iter_max<-500;
    omega<-0
    while( (iter<iter_max)&&(err>err_max) ){
      iter<-iter+1
      v01<-v0# v01 is the initial σ^2
      v1<-v# v1 is the initial σk^2
      b1<-b# b1 is the initial β
      #s1<-s
      try_a<-try({ vi<-chol2inv(chol(vv)) },silent=TRUE)# solve(V)
      if('try-error' %in% class(try_a)){
        try_aa<-try({ vi<-solve(vv) },silent=TRUE)
        if('try-error' %in% class(try_aa)){ vi<-ginv(vv) }
      }
      xtv<-crossprod(x,vi)# t(X)%*%solve(V)

      if(ncol(x)==1){
        b<-((xtv%*%x)^(-1))*(xtv%*%y)
      }else{
        if(abs(min(Mod(eigen(xtv%*%x)$values)))<1e-6){
          try_b<-try({  b<-chol2inv(chol((xtv%*%x)+diag(ncol(x))*1e-8))%*%(xtv%*%y)  },silent=TRUE)
          if('try-error' %in% class(try_b)){
            try_c<-try({  b<-solve((xtv%*%x))%*%(xtv%*%y)  },silent=TRUE)
            if('try-error' %in% class(try_c)){  b<-ginv((xtv%*%x))%*%(xtv%*%y)  }
          }
        }else{
          try_b<-try({ b<-chol2inv(chol(xtv%*%x))%*%(xtv%*%y) },silent=TRUE)
          if('try-error' %in% class(try_b)){
            try_c<-try({  b<-solve((xtv%*%x))%*%(xtv%*%y)  },silent=TRUE)
            if('try-error' %in% class(try_c)){  b<-ginv((xtv%*%x))%*%(xtv%*%y)  }
          }
        }
      }
      r<-y-x%*%b# y-Xβ
      ss<-matrix(0,nrow=n,ncol=1)
      vv<-matrix(0,n,n)# new V

      for(i in 1:kn){
        nc<-( (i-1)*mk+1 ):(i*mk)
        zz<-z[,nc]# Zk
        zztvi<-crossprod(zz,vi)# t(Zk)%*%solve(V)
        u[,i]<-v[i,]*zztvi%*%r# E(γk)
        w[,nc]<-v[i,]*( diag(1,mk)-zztvi%*%zz*v[i,] )# var(γk)
        s[i,]<-sum(diag(w[,nc]))# tr(var(γk))
        v[i,]<-(crossprod(u[,i,drop=F],u[,i,drop=F])+s[i,]+omega)/(tau+2+mk)
        ss<-ss+zz%*%u[,i,drop=F]
        vv<-vv+tcrossprod(zz,zz)*v[i,]# ∑( Zk%*%t(Zk)*(σk^2) )
      }
      v0<-as.numeric(crossprod(r,(r-ss))/n)# new σ^2
      vv<-vv+diag(n)*v0# new V

      err<-(crossprod((b1-b),(b1-b))+(v01-v0)^2+crossprod((v1-v),(v1-v)))/(1+ncol(x)+kn)
      beta<-t(b)
      sigma2<-v0
    }

    u1<-matrix(0,nrow=2,ncol=kn)# main-E(γk)
    p1<-matrix(1,kn,1)
    # pvalue<-matrix(1,kn,1)

    for(i in 1:kn){
      nc<-( (i-1)*mk+1 ):(i*mk)
      gammak<-u[,i,drop=F]

      u1[,i]<-L%*%gammak

      var_1<-L%*%w[,nc,drop=F]%*%t(L)
      tr_1<-sum(diag(ginv(tcrossprod(L))%*%var_1))##tr[...]
      dk1<-abs(rank_1-tr_1/v[i,])

      p1[i,]<-1-pchisq(  t(u1[,i,drop=F])%*%ginv(L%*%w[,nc]%*%t(L))%*%u1[,i,drop=F],     2)
    }
    return(list(b=b,u=u,u1=u1,sigma2=sigma2,p1=p1,iter=iter))
  }
  Zhou_sbl<-function(peak,Order,DesignMatrix,CodeMatrix,genoname,sbl_t,sbl_p,tau,err_max,fix_p,Sigma,SigmaK){

    chr_n<-as.numeric(genoname[,2])
    chr_n<-chr_n[!duplicated(chr_n)]
    maxdistance<-0
    for(i in 1:length(chr_n)){
      maxdistance<-max(maxdistance,   max(diff(as.matrix(genoname[which(genoname[,2]==i),3]))))
    }
    result_sblgwas<-sblgwas(x,y,DesignMatrix,sbl_t)
    sbl_p_wald<-result_sblgwas$blup[,4]#sbl_par<-result_sblgwas$blup[,1]# sbl_p_wald<-p.adjust(sbl_p_wald, method = "bonferroni")
    sbl_pos_order<-Order[order(sbl_p_wald)]
    p_order<-sort(sbl_p_wald)
    id1<-which( p_order< (1-pchisq(sbl_p*2*log(10),1)))
    id2<-which((p_order>=(1-pchisq(sbl_p*2*log(10),1)))&(p_order<1))

    if(length(id1)>0){

      sbl_pos_order1<-sbl_pos_order[id1]
      sbl_pos_order2<-sbl_pos_order[id2]
      sbl_pos_order1<-sbl_pos_order1[!duplicated(sbl_pos_order1)]
      sbl_pos_order2<-sbl_pos_order2[!duplicated(sbl_pos_order2)]

      sort_order1<-sort(sbl_pos_order1)
      result_emba<-ebayes_EM(x,multi_peak_new(CodeMatrix,sort_order1)$z,y,Sigma,SigmaK[sort_order1],tau,err_max)
      emID1<-which(result_emba$p1<(1-pchisq(fix_p*2*log(10),2)))
      emID2<-order(result_emba$p1)[seq(1,5,1)]

      if(length(emID1)>5){
        emID<-sort(emID1)
        fix_pos<-sort_order1[emID]

        if(maxdistance<=1){
          sbl_pos_order1<-selection2(sbl_pos_order1,fix_pos,genoname,1)
          sbl_pos_order1<-selection(sbl_pos_order1,genoname,1)
          sbl_pos_order2<-selection2(sbl_pos_order2,fix_pos,genoname,1)
          sbl_pos_order2<-selection2(sbl_pos_order2,sbl_pos_order1,genoname,1)
          sbl_pos_order2<-selection(sbl_pos_order2,genoname,1)
        }
      }else{
        emID<-sort(union(emID1,emID2))
        fix_pos<-sort_order1[emID]

        if(maxdistance<=1){
          sbl_pos_order1<-selection(sbl_pos_order1,genoname,1)
          sbl_pos_order2<-selection2(sbl_pos_order2,sbl_pos_order1,genoname,1)
          sbl_pos_order2<-selection(sbl_pos_order2,genoname,1)
        }
      }

      sbl_pos<-sort(c(sbl_pos_order1,sbl_pos_order2))# length(union(sbl_pos,fix_pos))

    }else{

      result_emba<-ebayes_EM(x,multi_peak_new(CodeMatrix,peak)$z,y,Sigma,SigmaK[peak],tau,err_max)
      emID1<-which(result_emba$p1<(1-pchisq(fix_p*2*log(10),2)))
      emID2<-order(result_emba$p1)[seq(1,5,1)]

      if(length(emID1)>5){
        emID<-sort(emID1)
        fix_pos<-peak[emID]

        if(maxdistance<=1){
          sbl_pos_order1<-selection2(peak,fix_pos,genoname,1)
          sbl_pos_order1<-selection(sbl_pos_order1,genoname,1)
        }else{
          sbl_pos_order1<-peak[-emID]
        }
      }else{
        emID<-sort(union(emID1,emID2))
        fix_pos<-peak[emID]

        if(maxdistance<=1){
          sbl_pos_order1<-selection2(peak,fix_pos,genoname,1)
          sbl_pos_order1<-selection(sbl_pos_order1,genoname,1)
        }else{
          sbl_pos_order1<-peak[-emID]
        }
      }

      sbl_pos<-sort(sbl_pos_order1)
    }

    sbl_fix_pos<-sort(fix_pos)
    xin<-cbind(x,multi_peak_new(CodeMatrix,sbl_fix_pos)$z)
    return(list(fix=sbl_fix_pos,pos=sbl_pos,xin=xin))
  }
  into_vector<-function(xmatrix){
    xvector<-NULL
    for(i in 1:dim(xmatrix)[2]){
      xvector<-c(xvector,xmatrix[,i])
    }
    return(xvector)
  }
  multinormal<-function(y,mean,sigma){
    pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
    return (pdf_value)
  }
  LRT_F2<-function(xxn,xxx,yn,par,mk){
    # mk<-2+2*(en-1)# the number of genotypes at per locus
    xn<-ncol(as.matrix(xxn))
    nq<-ncol(xxx)
    ns<-nrow(yn)
    kn<-nq/mk
    at1<-nq

    ad<-if(at1>0.5) cbind(xxn,xxx) else xxn
    if(length(par)==0){
      bb<-if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6) solve(crossprod(ad,ad)+diag(ncol(ad))*1e-8)%*%crossprod(ad,yn) else solve(crossprod(ad,ad))%*%crossprod(ad,yn)
    }else{
      bb<-par
    }
    vv<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns)##y-(X Z)t(β γ)
    ll<-sum(log(abs(multinormal(yn,ad%*%bb,vv))))
    lod<-matrix(0,kn,1)

    if(at1>0.5){
      for(m in 1:kn){
        i1<-(((m-1)*mk+1):(m*mk));# i2<-((m-1)*mk+1):((m-1)*mk+2);  i3<-((m-1)*mk+3):(m*mk)
        m1<-seq(1,ncol(ad),1)[-c(i1+xn)];# m2<-sub[-c(i2+xn)];  m3<-sub[-c(i3+xn)]

        ad1<-ad[,m1,drop=F]
        if(length(par)==0){
          bb1<-if(abs(min(eigen(crossprod(ad1,ad1))$values))<1e-6) solve(crossprod(ad1,ad1)+diag(ncol(ad1))*1e-8)%*%crossprod(ad1,yn) else solve(crossprod(ad1,ad1))%*%crossprod(ad1,yn)
        }else{
          bb1<-par[m1]
        }
        vv1<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns)
        ll1<-sum(log(abs(multinormal(yn,ad1%*%bb1,vv1))))
        lod[m,]<--2.0*(ll1-ll)/(2.0*log(10))
      }
    }

    return(lod)
  }
  optimize_every_posx<-function(xpos,z,yn,genoname,rr,tau,err_max){
    chr_n<-as.numeric(genoname[,2])
    chr_n<-chr_n[!duplicated(chr_n)]
    maxdistance<-0
    for(i in 1:length(chr_n)){
      maxdistance<-max(maxdistance,   max(diff(as.matrix(genoname[which(genoname[,2]==i),3]))))
    }
    if(maxdistance<rr){
      rr<-rr/maxdistance
      ad<-cbind(x,multi_peak_new(z,xpos)$z)
      if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6){
        try_a<-try({   bb<- chol2inv(chol(crossprod(ad,ad)+diag(ncol(ad))*1e-8))%*%crossprod(ad,yn)  },silent=TRUE)
        if('try-error' %in% class(try_a)){
          try_aa<-try({bb<- solve(crossprod(ad,ad))%*%crossprod(ad,yn)},silent=TRUE)
          if('try-error' %in% class(try_aa)){  bb<- ginv(crossprod(ad,ad))%*%crossprod(ad,yn)}
        }
      }else{
        try_a<-try({  bb<-chol2inv(chol(crossprod(ad,ad)))%*%crossprod(ad,yn) },silent=TRUE)
        if('try-error' %in% class(try_a)){
          try_aa<-try({bb<- solve(crossprod(ad,ad))%*%crossprod(ad,yn)},silent=TRUE)
          if('try-error' %in% class(try_aa)){  bb<- ginv(crossprod(ad,ad))%*%crossprod(ad,yn)}
        }
      }

      par<-bb[-c(1:dim(x)[2])]
      result_pos<-xpos
      chr_sum<-NULL
      for(i in 1:length(chr_n)){
        chr_sum<-c(chr_sum,length(which(genoname[,2]==i)))
      }
      chr_sum<-c(0,chr_sum)
      for(i in 1:length(xpos)){
        yy<-y-multi_peak_new(z,xpos[-i])$z%*%par[-seq((i-1)*3+1,i*3,1)]
        chr_now<-apply(genoname[,2,drop=F],2,as.numeric)[xpos[i]]
        if(i==1){
          left_rr<-min(xpos[i]-1-sum(chr_sum[seq(1,chr_now,1)]),rr)
        }else{
          if(genoname[xpos[i-1],2]==genoname[xpos[i],2]){
            left_rr<-min(0.5*(xpos[i]-xpos[i-1]),rr)
          }else{
            left_rr<-min(xpos[i]-1-sum(chr_sum[seq(1,chr_now,1)]),rr)
          }
        }
        if(i==length(xpos)){
          right_rr<-min(sum(chr_sum[seq(1,chr_now+1,1)])-xpos[i],rr)
        }else{
          if(genoname[xpos[i+1],2]==genoname[xpos[i],2]){
            right_rr<-min(0.5*(xpos[i+1]-xpos[i]),rr)
          }else{
            right_rr<-min(sum(chr_sum[seq(1,chr_now+1,1)])-xpos[i],rr)
          }
        }
        left_rr<-floor(left_rr)
        right_rr<-floor(right_rr)
        least_pos<-xpos[-i]
        now_pos<-c((xpos[i]-left_rr):(xpos[i]+right_rr))
        try_x<-try({
          result_embax<-ebayes_EM(x,multi_peak_new(z,now_pos)$z,yy,initial_sigma,initial_sigmak[now_pos],tau,err_max)
        },silent=TRUE)
        if('try-error' %in% class(try_x)){
          max_pos<-now_pos[which.min(result_embax$p1)]
          result_pos[i]<-max_pos# rm(result_embax)
        }
      }
    }else{
      result_pos<-xpos
    }

    return(result_pos)
  }
  multi_code_classic<-function(n_id,peak_id){
    mk<-2# the number of genotypes at per locus
    lengthpeak<-length(peak_id)
    gen_A<-(Ax-Bx)[n_id,peak_id,drop=F]
    gen_D<-Hx[n_id,peak_id,drop=F]
    adgen3<-matrix(0,nrow=n,ncol=lengthpeak*mk)
    adgen3[,seq(1,lengthpeak*mk,mk)]<-gen_A
    adgen3[,seq(2,lengthpeak*mk,mk)]<-gen_D
    return(adgen3)
  }
  effect_estimation<-function(n_id,xpos){
    xmatrix<-multi_code_classic(n_id,xpos)

    ad<-cbind(x,xmatrix)
    bb<-if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6) solve(crossprod(ad,ad)+diag(ncol(ad))*1e-8)%*%crossprod(ad,y) else solve(crossprod(ad,ad))%*%crossprod(ad,y)
    sig_e2<-as.numeric(crossprod(y-ad%*%bb)/length(y))
    bb<-bb[-1]

    Into_matrix<-function(vector_x,row_n){
      col_n<-length(vector_x)/row_n
      result_x<-matrix(0,nrow=row_n,ncol=col_n)
      for(i in 1:col_n){
        result_x[,i]<-vector_x[((i-1)*row_n+1):(i*row_n)]
      }
      return(result_x)
    }
    effect_all<-t(Into_matrix(bb,2))

    ef_Q<-effect_all
    sig_Q<-0.5*(ef_Q[,1])^2+0.25*(ef_Q[,2])^2
    sig_y<-max(var(y),(sum(sig_Q)+sig_e2))
    pve<-(sig_Q/sig_y)*100

    return(list(effect_all,pve,sig_Q,sig_e2,sig_y))
  }
  LeftRight_marker<-function(map,ChrPos){
    LR_result<-NULL
    for(i in 1:dim(ChrPos)[1]){
      now_id<-which(as.numeric(map[,2])==as.numeric(ChrPos[i,1]))
      now_pos<-as.numeric(ChrPos[i,2])
      all_pos<-as.numeric(map[now_id,3])
      if(now_pos<min(all_pos)){
        left_mar<-""
      }else{
        left_id<-max(which(all_pos<=now_pos))
        left_mar<-map[now_id,1][left_id]
      }
      if(now_pos>max(all_pos)){
        right_mar<-""
      }else{
        right_id<-min(which(all_pos>=now_pos))
        right_mar<-map[now_id,1][right_id]
      }
      LR_result<-rbind(LR_result,c(left_mar,right_mar))
    }
    return(LR_result)
  }
  ######################################### input and basic setup  #########################################
  #*#########  environment and phenotype  #
  pheno<-pheRaw[,NUM,drop=F]
  yes_id<-which(pheno!="-")
  y<-as.numeric(pheno[yes_id])
  n<-length(y)
  y<-as.matrix(y)
  #*#########  genotype  #
  # genRaw<-as.matrix(genRaw)

  #*#########  calculate Z matrix for K matrix  #
  mn<-dim(Ax0)[2]
  Z<-matrix(0,nrow=n,ncol=mn*3)
  Z[,seq(1,mn*3,3) ]<-Ax0[yes_id,]
  Z[,seq(2,mn*3,3) ]<-Hx0[yes_id,]
  Z[,seq(3,mn*3,3) ]<-Bx0[yes_id,]# dim(Z)
  #*#########  calculate K matrix  #
  K<-kinship_every(Z)

  #*#########  calculate Z matrix  for the subsequent algorithm #
  mn<-dim(Ax)[2]
  Z<-matrix(0,nrow=n,ncol=mn*3)
  Z[,seq(1,mn*3,3) ]<-Ax[yes_id,]
  Z[,seq(2,mn*3,3) ]<-Hx[yes_id,]
  Z[,seq(3,mn*3,3) ]<-Bx[yes_id,]# dim(Z)

  #*#########  X matrix;     y=Xβ+Zγ+ξ+ε  #
  x<-matrix(1,nrow=n,ncol=1)#
  if(is.null(yygg)==FALSE){
    x<-cbind(x,yygg[yes_id,,drop=F])
  }# dim(x)
  if(det(crossprod(x,x))==0){
    warning("X is singular")
  }
  ReduceDim_x<-TRUE
  if(ReduceDim_x){
    x_effect<-if(abs(min(eigen(crossprod(x,x))$values))<1e-6) solve(crossprod(x,x)+diag(ncol(x))*1e-8)%*%crossprod(x,y) else solve(crossprod(x,x))%*%crossprod(x,y)
    yygg_effect<-x_effect[-1,1,drop=F]
    y<-y-x[,-1,drop=F]%*%yygg_effect
    x<-matrix(1,nrow=n,ncol=1)
  }
  #*#########  name  #

  ######################################### single_locus_scanning  #########################################
  #*#########  single locus scanning  #
  p3d_result<-p3d_method(x,y,K)

  if(Model=="Random"){
    single_locus_model_result<-single_locus_model(x=x,y=y,zz=Z,lambda=p3d_result$lambda,uu=p3d_result$k_vector,vv=p3d_result$k_value,CLO=CLO)
    initial_sigma<-mean(single_locus_model_result$sigma_2)
    initial_sigmak<-single_locus_model_result$phi_k
  }else if(Model=="Fixed"){
    single_locus_model_result<-single_locus_model_Fixed(x=x,y=y,zz=Z,lambda=p3d_result$lambda,uu=p3d_result$k_vector,vv=p3d_result$k_value,CLO=CLO)
    initial_sigma<-mean(single_locus_model_result$sigma_2)
    initial_sigmak<-rep(1,mn)
  }else{
    warning("Please enter Model!")
  }

  #*#########  pick the peaks  #
  peak_id<-peak_selection(single_locus_model_result$log1,genoname)# length(peak_id)
  ######################################### multi_locus_scanning  #########################################
  multi_locus_result1<-Zhou_lars(peak_id,Z,n) # length(multi_locus_result1$lar_pos)

  multi_locus_result2<-Zhou_sbl(peak=multi_locus_result1$lar_pos,Order=multi_locus_result1$order,
                                DesignMatrix=multi_locus_result1$Matrix,CodeMatrix=Z,
                                genoname=genoname,
                                sbl_t=-1,sbl_p=3,tau=0,err_max=1e-6,fix_p=1.5,
                                Sigma=initial_sigma,SigmaK=initial_sigmak)# larpos=multi_locus_result1$lar_pos0,larbeta=multi_locus_result1$beta
  emba_p<-3
  result_emba<-ebayes_EM(multi_locus_result2$xin,multi_peak_new(Z,multi_locus_result2$pos)$z,
                         y,initial_sigma,initial_sigmak[multi_locus_result2$pos],
                         tau=-2,err_max=1e-8)
  emba_pos0<-which(result_emba$p1<(1-pchisq(emba_p*2*log(10),2)))# cbind(result_emba$p1,result_emba$p2)
  emba_all_pos<-sort(c(multi_locus_result2$fix,multi_locus_result2$pos[emba_pos0]))


  result_emba1<-ebayes_EM(x,multi_peak_new(Z,emba_all_pos)$z,
                          y,initial_sigma,initial_sigmak[emba_all_pos],
                          tau=-2,err_max=1e-8)
  emba1_pos<-emba_all_pos
  emba1_par_E<-c(result_emba1$b)
  emba1_par<-result_emba1$u
  emba1_par<-into_vector(emba1_par)

  multi_value4<-multi_peak_new(Z,emba1_pos)
  z_M4<-multi_value4$z
  order4<-multi_value4$order
  LRT_lod<-LRT_F2(xxn=x,xxx=z_M4, yn=y,par=c(emba1_par_E,emba1_par),mk=1)# cbind(order4,LRT_lod)
  lrt_pos<-order4[which(LRT_lod>2.5)]
  lrt_pos<-lrt_pos[!duplicated(lrt_pos)]# length(lrt_pos)
  ######################################### Optimization and output #########################################
  if(length(lrt_pos)>0){
    if(CriDis<=4){
      optimize_pos<-optimize_every_posx(xpos=lrt_pos,z=Z,yn=y,genoname,rr=CriDis,tau=0,err_max=1e-8)
    }else{
      optimize_pos<-lrt_pos
    }
    emba3_pos<-optimize_pos

    # CriLOD<-3
    lod_Q<-LRT_F2(xxn=x,xxx=multi_peak_new(Z,emba3_pos)$z, yn=y,par=NULL,mk=3)# cbind(emba3_pos,lod_Q)
    lrt2_pos<-emba3_pos[which(lod_Q>=CriLOD)]# length(lrt2_pos)
    last_lod<-lod_Q[which(lod_Q>=CriLOD)]

    if(length(lrt2_pos)>0){
      IC_data<-cbind(x,multi_peak_new(Z,lrt2_pos)$z)
      lm_IC<-lm(y~IC_data-1)
      AIC(lm_IC)
      BIC(lm_IC)

      LR_marker<-LeftRight_marker(map=mapRaw,ChrPos=genoname[lrt2_pos,2:3,drop=F])
      result_all<-effect_estimation(yes_id,lrt2_pos)
      var_e<-matrix("",nrow=length(lrt2_pos),ncol=1)
      var_y<-matrix("",nrow=length(lrt2_pos),ncol=1)
      var_e[1]<-round(result_all[[4]],4)
      var_y[1]<-round(result_all[[5]],4)

      data.all<-data.frame(genoname[lrt2_pos,2:3,drop=F],
                           round(result_all[[1]],4),round(last_lod,4),
                           LR_marker,
                           round(result_all[[3]],4),
                           round(result_all[[2]],4),
                           var_e,var_y)
      # rep(AIC(lm_IC),length(lrt2_pos)),
      # rep(BIC(lm_IC),length(lrt2_pos)))
      rownames(data.all)<-NULL
      colnames(data.all)<-c("Chr","Position(cM)","Effect.a","Effect.d","LOD",
                            "Left_marker","right_marker",
                            "Var_Genet","r2(%)",
                            "Var_Error","Var_Phen(total)")
      reslt_list<-list(result=data.all,p_Q=single_locus_model_result$log1)
    }else{
      reslt_list<-NULL
      warning("No QTL were detected!")
    }
  }else{
    reslt_list<-NULL
    warning("No QTL were detected!")
  }

  return(reslt_list)
}
#' The second step of Zhou method for multiple environments
#'
#' @param Model Random or fixed model.
#' @param pheRaw phenotype matrix.
#' @param genRaw genotype matrix.
#' @param mapRaw linkage map matrix.
#' @param CriLOD Critical LOD scores for significant QTL.
#' @param NUM The serial number of the trait to be analyzed.
#' @param EnvNum The number of environments for each trait is a vector.
#' @param yygg covariate matrix.
#' @param genoname linkage map matrix with pseudo markers inserted.
#' @param Ax0 AA genotype matrix.
#' @param Hx0 Aa genotype matrix.
#' @param Bx0 aa genotype matrix.
#' @param Ax AA genotype matrix with pseudo markers inserted.
#' @param Hx Aa genotype matrix with pseudo markers inserted.
#' @param Bx aa genotype matrix with pseudo markers inserted.
#' @param dir file storage path.
#' @param CriDis The distance of optimization.
#' @param CLO Number of CPUs.
#'
#' @return a list
#' @export
#'
#' @examples
#' data(F2data)
#' readraw<-Readdata(file=F2data,fileFormat="GCIM",
#' method="GCIM-QEI",filecov=NULL,
#' MCIMmap=NULL,MultiEnv=TRUE)
#' DoResult<-Dodata(fileFormat="GCIM",
#' Population="F2",method="GCIM-QEI",
#' Model="Random",readraw,MultiEnv=TRUE)
#' ZhouMatrices<-ZhouF(pheRaw=DoResult$pheRaw,
#' genRaw=DoResult$genRaw,mapRaw1=DoResult$mapRaw1,
#' WalkSpeed=1,CriLOD=3,dir=tempdir())
#' OutputZhou<-ZhouMethod(Model="Random",
#' pheRaw=DoResult$pheRaw,genRaw=DoResult$genRaw,
#' mapRaw=ZhouMatrices$mapRaw,CriLOD=3,NUM=1,
#' EnvNum=DoResult$EnvNum,yygg=DoResult$yygg1,
#' genoname=ZhouMatrices$genoname,
#' Ax0=ZhouMatrices$Ax0,Hx0=ZhouMatrices$Hx0,
#' Bx0=ZhouMatrices$Bx0,Ax=ZhouMatrices$Ax,
#' Hx=ZhouMatrices$Hx,Bx=ZhouMatrices$Bx,
#' dir=tempdir(),CriDis=5,CLO=2)
ZhouMethod<-function(Model=NULL,pheRaw=NULL,genRaw=NULL,mapRaw=NULL,CriLOD=NULL,NUM=NULL,EnvNum=NULL,yygg=NULL,genoname=NULL,
                     Ax0=NULL,Hx0=NULL,Bx0=NULL,Ax=NULL,Hx=NULL,Bx=NULL,dir=NULL,CriDis=NULL,CLO=NULL){# chr_n=NULL,

  #########################################  function  #########################################
  kinship_all<-function(coded_gen,n_id,en){
    kinship_every<-function(coded_gen){
      kk<-coded_gen%*%t(coded_gen)
      kk<-kk/(dim(coded_gen)[2]/en/3)
      return(kk)
    }
    k_all<-matrix(0,n,n)
    sum_n<-0
    for (i in 1:en){
      row_col<-(sum_n+1):(sum_n+length(n_id[[i]]))
      k_all[row_col,row_col]<-kinship_every(coded_gen[(sum_n+1):(sum_n+length(n_id[[i]])),])
      sum_n<-sum_n+length(n_id[[i]])
    }
    return(k_all)
  }
  fixed_x<-function(n_id,en){
    x0<-matrix(1,n,1)
    col.E<-matrix(0,nrow = n,ncol = en-1)
    col.E[(n-length(n_id[[en]])+1):n,]<--1
    sum_n<-0
    for(i in 1:(en-1)){
      col.E[(sum_n+1):(sum_n+length(n_id[[i]])),i]<-1
      sum_n<-sum_n+length(n_id[[i]])
    }
    x<-cbind(x0,col.E)
    return(x)
  }
  name_function<-function(en){
    effect_name<-NULL
    for(i in 1:en){
      effect_name<-c(effect_name,paste("Effect.aE",i,sep = ""))
      effect_name<-c(effect_name,paste("Effect.dE",i,sep = ""))
    }
    effect_name<-c("Effect.a","Effect.d",effect_name)
    return(effect_name)
  }
  p3d_method<-function(x,y,kinship){
    # estimate the value of λ & λk        kinship=K;
    P3D<-function(x,y,kinship){
      iter_p3d<-function(ga){
        lambda<-exp(ga)
        diag_element<-lambda*K_value+1
        logH<-sum(log(diag_element))
        RH_value<-1/(diag_element)
        yuRHyu<-sum(yu*RH_value*yu)
        yuRHxu<-matrix(0,nrow = 1,ncol = q)
        xuRHxu<-matrix(0,nrow = q,ncol = q)
        for(i in 1:q){
          yuRHxu[,i]<-sum(yu*RH_value*xu[,i])
          for(j in 1:q){
            xuRHxu[i,j]<-sum(xu[,i]*RH_value*xu[,j])
          }
        }
        logxuRHxu<-log(det(xuRHxu))
        logyuPyu<-log(yuRHyu-yuRHxu%*%tcrossprod(solve(xuRHxu),yuRHxu))
        output<- -0.5*(logH+logxuRHxu+(n-q)*logyuPyu)
        return(-output)
      }

      q<-ncol(x)
      n<-nrow(y)
      eigenK<-eigen(kinship)
      K_vector<-eigenK$vectors
      K_value<-eigenK$values# rm(eigenK);gc()
      xu<-crossprod(K_vector,x)
      yu<-crossprod(K_vector,y)

      ga0<-0
      optimp3d<-optim(par=ga0,fn=iter_p3d,hessian = TRUE,method="L-BFGS-B",lower=-50,upper=10)
      lambda<-exp(optimp3d$par)

      return(list(lambda,K_vector,K_value))
    }

    q<-ncol(x)
    value1<-P3D(x=x,y=y,kinship=kinship)
    lambda<-value1[[1]]
    uu<-as.matrix(value1[[2]])
    vv<-value1[[3]]
    # RH_value<-1/(vv*lambda+1)
    rm(value1);gc()
    return(list(lambda=lambda,k_vector=uu,k_value=vv))
  }
  single_locus_model<-function(x,y,zz,lambda,uu,vv,en,CLO){
    # genotype effect transform to additive dominance(transformation matrix)
    L_coefficient<-function(en){
      e.seq<-rep(1,3*en)
      a11<-matrix(0,1,3*en)
      a12<-matrix(0,1,3*en)
      a13<-matrix(0,1,3*en)
      a.seq<-seq(1,3*en,by=3)
      a11[a.seq]<-1
      a12[a.seq+1]<-1
      a13[a.seq+2]<-1
      a123<-rbind(a11,a12,a13)

      a1<-1/en*a11-1/(3*en)*e.seq
      a2<-1/en*a12-1/(3*en)*e.seq
      a3<-1/en*a13-1/(3*en)*e.seq

      L1<-rbind(a1,a2,a3)
      a4<-matrix(c(0.5,-0.5,0,1,-0.5,-0.5),2,3)
      LL1<-a4%*%L1

      L2<-matrix(0,en,3*en)
      L3<-matrix(0,3*en,3*en)
      c4<-matrix(0,2*en,3*en)
      for(i in 1:en){

        b11<-matrix(0,1,3*en)
        b11[((i-1)*3+1):(i*3)]<-1
        L2[i,]<-1/3*b11-1/(3*en)*e.seq

        c4[((i-1)*2+1):(i*2),((i-1)*3+1):(i*3)]<-a4

        for(i0 in 1:3){
          seq.c<-(i-1)*3+i0
          c0<-matrix(0,1,3*en)
          c0[seq.c]<-1
          L3[seq.c,]<--1/en*a123[i0,]-1/3*b11+1/(3*en)*e.seq+c0
        }
      }
      LL3<-c4%*%L3
      return(list(matrix_C1=L1, matrix_C2=L2, matrix_C3=L3 , LL1=LL1, LL3=LL3, L1=a4, L3=c4))
    }
    value2<-L_coefficient(en)# L coefficient matrix # rm(L_coefficient);gc()

    # iteration function(R=zu%*%t(zu)*λk+D*λ+I);estimate λk
    rqtl<-function(ga){
      lambdak<-exp(ga)
      Hk_term<-zuRHzu*lambdak+diag(1,en*3)
      logHk<-sum(log(lambda*vv+1))+log(det(Hk_term))
      RHk_term<-solve(Hk_term)*lambdak
      yuRHkyu<-yuRHyu-yuRHzu%*%tcrossprod(RHk_term,yuRHzu)
      yuRHkxu<-yuRHxu-yuRHzu%*%tcrossprod(RHk_term,xuRHzu)
      xuRHkxu<-xuRHxu-xuRHzu%*%tcrossprod(RHk_term,xuRHzu)
      yuPkyu<-yuRHkyu-yuRHkxu%*%tcrossprod(solve(xuRHkxu),yuRHkxu)
      rqtl<- -0.5*( logHk + log(det(xuRHkxu)) + (n-q)*log(yuPkyu) )
      return(-rqtl)
    }

    # estimate the value of γ
    gamma_estimate<-function(lambdak){
      Hk_term<-zuRHzu*lambdak+diag(1,en*3)
      RHk_term<-solve(Hk_term)*lambdak
      yuRHkyu<-yuRHyu-yuRHzu%*%tcrossprod(RHk_term,yuRHzu)
      yuRHkxu<-yuRHxu-yuRHzu%*%tcrossprod(RHk_term,xuRHzu)
      xuRHkxu<-xuRHxu-xuRHzu%*%tcrossprod(RHk_term,xuRHzu)
      zuRHkxu<-t(xuRHzu)-zuRHzu%*%tcrossprod(RHk_term,xuRHzu)
      zuRHkyu<-t(yuRHzu)-zuRHzu%*%tcrossprod(RHk_term,yuRHzu)
      zuRHkzu<-zuRHzu-zuRHzu%*%RHk_term%*%zuRHzu
      beta<-solve(xuRHkxu,t(yuRHkxu))
      beta<-matrix(beta,ncol = 1)
      yuPkyu<-yuRHkyu-yuRHkxu%*%tcrossprod(solve(xuRHkxu),yuRHkxu)
      sigma<-yuPkyu/(n-q)
      sigma<-as.numeric(sigma)
      gamma<-lambdak*zuRHkyu-lambdak*zuRHkxu%*%beta
      var<-abs((lambdak*diag(1,en*3)-lambdak*zuRHkzu*lambdak)*sigma)
      stderr<-sqrt(diag(var))
      phi.k<-sigma*lambdak # Φk
      return(list(gamma,beta,var,phi.k,sigma,stderr))
    }

    # estimate the value of logp
    logp_estimate<-function(L, g_k, pn){
      var.1<-L%*%tcrossprod(var_ga,L)
      Wk.1<-crossprod(g_k,ginv(var.1))%*%g_k
      rank1<-qr(L)$rank

      tr1<-sum(diag(ginv(tcrossprod(L,L))%*%var.1))
      dk1<-rank1-tr1/phi_k

      p1<-pgamma(Wk.1,shape=pn/2,scale=2*dk1,lower.tail = FALSE,log.p = FALSE)
      log1<-pgamma(Wk.1,shape=pn/2,scale=2*dk1,lower.tail = FALSE,log.p = TRUE)
      log1<--log1*log10(exp(1))
      return(list(p=p1,log=log1))
    }

    RH_value<-1/(vv*lambda+1)
    mn<-ncol(zz)/(en*3)
    n<-nrow(y)
    q<-ncol(x)

    xu<-crossprod(uu,x)
    yu<-crossprod(uu,y)

    yuRHyu<-sum(yu*RH_value*yu)
    yuRHxu<-matrix(0,nrow=1,ncol=q)
    xuRHxu<-matrix(0,nrow=q,ncol=q)
    for(i in 1:q){
      yuRHxu[,i]<-sum(yu*RH_value*xu[,i])
      for(j in 1:q){
        xuRHxu[j,i]<-sum(xu[,j]*RH_value*xu[,i])
      }
    }
    # yuRHyu<-crossprod(yu,diag(RH_value))%*%yu
    # yuRHxu<-crossprod(yu,diag(RH_value))%*%xu
    # xuRHxu<-crossprod(xu,diag(RH_value))%*%xu

    if(is.null(CLO)==TRUE){
      cl.cores <- detectCores()
      if(cl.cores<=2){
        cl.cores<-1
      }else{
        if(cl.cores>10){
          cl.cores <-10
        }else{
          cl.cores <- detectCores()-1
        }
      }
    }else{
      cl.cores <-CLO
    }
    cl <- makeCluster(cl.cores)
    registerDoParallel(cl)

    MR_i<-numeric()
    result<-foreach(MR_i=1:mn,.combine=rbind)%dopar%{
      # library(MASS)
      z<-zz[,((MR_i-1)*3*en+1):(MR_i*3*en),drop=F]
      zu<-crossprod(uu,z)

      xuRHzu<-matrix(0,nrow=q,ncol=3*en)
      yuRHzu<-matrix(0,nrow=1,ncol=3*en)
      zuRHzu<-matrix(0,nrow=3*en,ncol=3*en)
      for(i in 1:(3*en)){
        yuRHzu[,i]<-sum(yu*RH_value*zu[,i])
        for(j in 1:q){
          xuRHzu[j,i]<-sum(xu[,j]*RH_value*zu[,i])
        }
        for(j in 1:(3*en)){
          zuRHzu[j,i]<-sum(zu[,j]*RH_value*zu[,i])
        }
      }
      # xuRHzu<-crossprod(xu,diag(RH_value))%*%zu
      # yuRHzu<-crossprod(yu,diag(RH_value))%*%zu
      # zuRHzu<-crossprod(zu,diag(RH_value))%*%zu
      ga<-0
      par<-optim(par=ga,fn=rqtl,hessian = TRUE,method="L-BFGS-B",lower=-10,upper=10)
      lambdak<-exp(par$par)

      value3<-gamma_estimate(lambdak)
      gamma_k<-value3[[1]]
      beta<-value3[[2]]# estimate β
      var_ga<-value3[[3]]
      phi_k<-value3[[4]]
      sigma_2<-value3[[5]]
      # stderr<-value3[[6]]

      gamma_k2<-gamma_k
      gamma_main_k <-value2$matrix_C1%*%gamma_k2
      gamma_env_k  <-value2$matrix_C2%*%gamma_k2
      gamma_inter_k<-value2$matrix_C3%*%gamma_k2

      main_effect<-value2$L1%*%gamma_main_k
      interact_effect<-value2$L3%*%gamma_inter_k

      logvalue1<-logp_estimate(L=value2$LL1, g_k=main_effect, pn=2)# the -log(10)p of qtl effect
      p1<-logvalue1$p
      log1<-logvalue1$log

      logvalue2<-logp_estimate(L=value2$matrix_C2, g_k=gamma_env_k, pn=en-1)# the -log(10)p of environment effect
      p2<-logvalue2$p
      log2<-logvalue2$log

      logvalue3<-logp_estimate(L=value2$LL3, g_k=interact_effect, pn=2*(en-1))# the -log(10)p of interaction effect
      p3<-logvalue3$p
      log3<-logvalue3$log

      result<-cbind(lambdak,matrix(beta,1,en),matrix(gamma_k,1,3*en),p1,p2,p3,log1,log2,log3,phi_k,sigma_2)

    }
    stopCluster(cl)
    # rm(RH_value);gc()

    lambda_k<-result[,1]
    mu_beta<-result[,2:(1+en)]
    gamma_all<-result[,(2+en):(1+4*en)]
    p1<-result[,(2+4*en)]
    p2<-result[,(3+4*en)]
    p3<-result[,(4+4*en)]
    log_p1<-result[,(5+4*en)]
    log_p2<-result[,(6+4*en)]
    log_p3<-result[,(7+4*en)]
    phi_k<-result[,(8+4*en)]
    sigma_2<-result[,(9+4*en)]

    return(list(lambda_k=lambda_k, fixed=mu_beta, gamma=gamma_all,
                p1=p1, p2=p2, p3=p3,
                log1=log_p1, log2=log_p2, log3=log_p3,
                phi_k=phi_k, sigma_2=sigma_2))

  }
  single_locus_model_Fixed<-function(x,y,zz,lambda,uu,vv,en,CLO){
    # genotype effect transform to additive dominance(transformation matrix)
    L_coefficient<-function(en){
      e.seq<-rep(1,3*en)
      a11<-matrix(0,1,3*en)
      a12<-matrix(0,1,3*en)
      a13<-matrix(0,1,3*en)
      a.seq<-seq(1,3*en,by=3)
      a11[a.seq]<-1
      a12[a.seq+1]<-1
      a13[a.seq+2]<-1
      a123<-rbind(a11,a12,a13)

      a1<-1/en*a11-1/(3*en)*e.seq
      a2<-1/en*a12-1/(3*en)*e.seq
      a3<-1/en*a13-1/(3*en)*e.seq

      L1<-rbind(a1,a2,a3)
      a4<-matrix(c(0.5,-0.5,0,1,-0.5,-0.5),2,3)
      LL1<-a4%*%L1

      L2<-matrix(0,en,3*en)
      L3<-matrix(0,3*en,3*en)
      c4<-matrix(0,2*en,3*en)
      for(i in 1:en){

        b11<-matrix(0,1,3*en)
        b11[((i-1)*3+1):(i*3)]<-1
        L2[i,]<-1/3*b11-1/(3*en)*e.seq

        c4[((i-1)*2+1):(i*2),((i-1)*3+1):(i*3)]<-a4

        for(i0 in 1:3){
          seq.c<-(i-1)*3+i0
          c0<-matrix(0,1,3*en)
          c0[seq.c]<-1
          L3[seq.c,]<--1/en*a123[i0,]-1/3*b11+1/(3*en)*e.seq+c0
        }
      }
      LL3<-c4%*%L3
      return(list(matrix_C1=L1, matrix_C2=L2, matrix_C3=L3 , LL1=LL1, LL3=LL3, L1=a4, L3=c4))
    }
    value2<-L_coefficient(en)# L coefficient matrix # rm(L_coefficient);gc()

    logp_estimate<-function(L, g_k, pn){
      var.1<-L%*%tcrossprod(var_ga,L)
      Wk.1<-crossprod(g_k,ginv(var.1))%*%g_k

      p1<-pchisq(Wk.1,df=pn,lower.tail = FALSE,log.p = FALSE)
      log1<-pchisq(Wk.1,df=pn,lower.tail = FALSE,log.p = TRUE)
      log1<--log1*log10(exp(1))
      return(list(p=p1,log=log1))
    }

    yu<-crossprod(uu,y)
    Hsolve<-1/(vv*lambda+1)

    if(is.null(CLO)==TRUE){
      cl.cores <- detectCores()
      if(cl.cores<=2){
        cl.cores<-1
      }else{
        if(cl.cores>10){
          cl.cores <-10
        }else{
          cl.cores <- detectCores()-1
        }
      }
      # cl.cores<-2
    }else{
      cl.cores <-CLO
    }
    cl <- makeCluster(cl.cores)
    registerDoParallel(cl)

    MF_i<-numeric()
    result<-foreach(MF_i=1:mn,.combine=rbind)%dopar%{
      # library(MASS)
      z<-zz[,((MF_i-1)*3*en+1):(MF_i*3*en),drop=F]
      uxz<-crossprod(uu,cbind(x,z))
      x_gamma<-ginv(t(uxz)%*%diag(Hsolve)%*%uxz)%*%t(uxz)%*%diag(Hsolve)%*%yu
      q<-qr(uxz)$rank
      sig_e2<-as.numeric(t(yu-uxz%*%x_gamma)%*%diag(Hsolve)%*%(yu-uxz%*%x_gamma)/(dim(uxz)[1]-q))
      x_gamma_covmatr<-sig_e2*ginv(t(uxz)%*%diag(Hsolve)%*%uxz)
      gamma<-x_gamma[-c(1:dim(x)[2])]
      var_ga<-x_gamma_covmatr[-c(1:dim(x)[2]),-c(1:dim(x)[2]),drop=F]

      gamma_main_k <-value2$matrix_C1%*%gamma
      gamma_env_k  <-value2$matrix_C2%*%gamma
      gamma_inter_k<-value2$matrix_C3%*%gamma

      main_effect<-value2$L1%*%gamma_main_k
      interact_effect<-value2$L3%*%gamma_inter_k

      logvalue1<-logp_estimate(L=value2$LL1, g_k=main_effect, pn=2)#the -log(10)p of qtl effect
      p1<-logvalue1$p
      log1<-logvalue1$log

      logvalue2<-logp_estimate(L=value2$matrix_C2, g_k=gamma_env_k, pn=en-1)#the -log(10)p of environment effect
      p2<-logvalue2$p
      log2<-logvalue2$log

      logvalue3<-logp_estimate(L=value2$LL3, g_k=interact_effect, pn=2*(en-1))#the -log(10)p of interaction effect
      p3<-logvalue3$p
      log3<-logvalue3$log

      result<-cbind(matrix(x_gamma[c(1:dim(x)[2])],1,en),matrix(gamma,1,3*en),p1,p2,p3,log1,log2,log3,sig_e2)
    }
    stopCluster(cl)
    # rm(RH_value);gc()

    mu_beta<-result[,1:dim(x)[2]]
    gamma_all<-result[,(dim(x)[2]+1):(dim(x)[2]+3*en)]
    p1<-result[,(dim(x)[2]+3*en+1)]
    p2<-result[,(dim(x)[2]+3*en+2)]
    p3<-result[,(dim(x)[2]+3*en+3)]
    log_p1<-result[,(dim(x)[2]+3*en+4)]
    log_p2<-result[,(dim(x)[2]+3*en+5)]
    log_p3<-result[,(dim(x)[2]+3*en+6)]
    sigma_2<-result[,(dim(x)[2]+3*en+7)]

    return(list(fixed=mu_beta, gamma=gamma_all,
                p1=p1, p2=p2, p3=p3,
                log1=log_p1, log2=log_p2, log3=log_p3,
                sigma_2=sigma_2))

  }
  peak_selection<-function(log_value,genoname){
    peak_pos<-function(Lod.temp){
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
    chr_all<-as.matrix(genoname[,2])
    chr_kind<-chr_all[!duplicated(chr_all)]
    id_pos<-NULL
    for(jjj in chr_kind){
      now_id<-which(chr_all%in%jjj)
      id_pos<-c(id_pos,now_id[peak_pos(log_value[now_id])])
    }
    return(sort(id_pos))
  }
  multi_peak_new<-function(gencoded,peak_id,en){
    enk<-3*en
    mut_peak_id<-NULL
    term1<-seq(1,3*en,1)
    for(i in 1:length(peak_id)){
      mut_peak_id<-c(mut_peak_id,(rep(peak_id[i],enk)-1)*enk+term1)
    }
    return(list(z=gencoded[,sort(mut_peak_id)],order=sort(rep(peak_id,enk))))
  }
  Zhou_lars<-function(peak,CodeMatrix,n,en){
    multi_value<-multi_peak_new(CodeMatrix,peak,en)
    DesignMatrix<-multi_value$z
    order0<-multi_value$order# length(order0); length(peak_id)
    if(length(peak)>=n){
      lar_result<-lars(x=cbind(x[,-1,drop=F],DesignMatrix), y=y, type = "lar",trace = FALSE, normalize = TRUE, intercept = TRUE, eps = .Machine$double.eps, use.Gram=FALSE)

      lar_result.0<-lar_result$beta[nrow(lar_result$beta),][-c(1:(dim(x)[2]-1))]
      lar_pos0<-order0[which(lar_result.0!=0)]
      lar_pos<-lar_pos0[!duplicated(lar_pos0)]# length(lar_pos)
      multi_value1<-multi_peak_new(CodeMatrix,lar_pos,en)
      DesignMatrix1<-multi_value1$z
      order1<-multi_value1$order
    }else{
      lar_pos<-peak
      DesignMatrix1<-DesignMatrix
      order1<-order0
    }# length(lar_pos)
    return(list(lar_pos=lar_pos,Matrix=DesignMatrix1,order=order1))
  }
  sblgwas<-function(x,y,z,t,max.iter=200,min.err=1e-6){
    x<-as.matrix(x)
    y<-as.matrix(y)
    z<-as.matrix(z)
    n<-length(y)
    q<-ncol(x)
    m<-ncol(z)
    b0<-solve(t(x)%*%x,tol=1e-50)%*%(t(x)%*%y)
    s2<-sum((y-x%*%b0)^2)/(n-q)
    b0<-matrix(0,q,1)
    b<-b0
    g0<-matrix(0,m,1)
    g<-g0
    lambda<-matrix(0,m,1)
    tau<-g0
    v<-g0
    xx<-NULL
    xy<-NULL
    for(i in 1:q){
      xx<-c(xx,sum(x[,i]^2))
      xy<-c(xy,sum(x[,i]*y))
    }
    zz<-NULL
    zy<-NULL
    for(k in 1:m){
      zz<-c(zz,sum(z[,k]^2))
      zy<-c(zy,sum(z[,k]*y))
    }
    d<-numeric(m)
    a<-matrix(0,n,1)
    iter<-0
    err<-1e8
    my.iter<-NULL
    while(iter < max.iter & err > min.err){
      for(i in 1:q){
        a<-a-x[,i]*b0[i]
        ai<-sum(x[,i]*a)
        b[i]<-(xy[i]-ai)/xx[i]
        a<-a+x[,i]*b[i]
      }
      df<-0
      for(k in 1:m){
        a<-a-z[,k]*g0[k]
        ak<-sum(z[,k]*a)
        c1<- -(t+3)*zz[k]^2
        c2<- -(2*t+5)*zz[k]+(zy[k]-ak)^2
        c3<- -(t+2)
        if( ((c2^2-4*c1*c3) < 0) | (c2 < 0) ){
          tau[k]<-0
        } else {
          tau[k]<-(-c2-sqrt(c2^2-4*c1*c3))/(2*c1)
        }
        lambda[k]<-tau[k]/s2
        g[k]<-lambda[k]*(zy[k]-ak)-lambda[k]^2*zz[k]*(zy[k]-ak)/(lambda[k]*zz[k]+1)
        d[k]<-lambda[k]*(zz[k]-lambda[k]*zz[k]^2/(lambda[k]*zz[k]+1))
        v[k]<-tau[k]-tau[k]*d[k]
        df<-df+d[k]
        a<-a+z[,k]*g[k]
      }

      if((n-q-df) > 0){s2<-sum((y-a)^2)/(n-q-df)
      }else{
        s2<-sum((y-a)^2)/(n-q)
      }

      iter<-iter+1
      err<-sum((g-g0)^2)/m
      g0<-g
      b0<-b
      my.iter<-rbind(my.iter,cbind(iter,err,s2,t(b),t(g)))
    }
    my.parm<-data.frame(iter,err,s2,b,df)
    names(my.parm)<-c("iter","error","s2","beta","df")

    posv<-which(v!=0)
    m<-length(g)
    wald<-c(rep(0,m))
    gg<-g[posv]
    vv<-v[posv]
    wald[posv]<-gg^2/vv
    p<-pchisq(wald,1,lower.tail=FALSE)

    my.blup<-data.frame(g,v,wald,p)
    names(my.blup)<-c("gamma","vg","wald","p_wald")

    var.beta<-NULL
    for(i in 1:q){
      var.beta<-c(var.beta,paste("beta",i,sep=""))
    }
    var.gamma<-NULL
    for(k in 1:m){
      var.gamma<-c(var.gamma,paste("gamma",k,sep=""))
    }
    var.names<-c(c("iter","error","s2"),var.beta,var.gamma)
    my.iter<-data.frame(my.iter)
    names(my.iter)<-var.names

    out<-list(my.iter,my.parm,my.blup)
    names(out)<-c("iteration","parm","blup")
    return(out)
  }
  selection<-function(posx,genoname,svrad){
    chose_peak<-c(posx[1])
    order_now<-1
    while(order_now<length(posx)){
      order_now<-order_now+1
      repeat_pos<-which( abs(chose_peak-as.numeric(posx[order_now]))<=(svrad)  )
      if(length(repeat_pos)>0){
        if_condition<-length(which(  genoname[chose_peak[repeat_pos],2]==as.numeric(genoname[posx[order_now],2])  ))==0
        if(if_condition){
          chose_peak<-c(chose_peak,posx[order_now])
        }
      }else{
        chose_peak<-c(chose_peak,posx[order_now])
      }
    }
    return(chose_peak)
  }
  selection2<-function(posx1,posx2,genoname,svrad){
    chose_peak<-NULL
    order_now<-0
    while(order_now<length(posx1)){
      order_now<-order_now+1
      repeat_pos<-which( abs(posx2-as.numeric(posx1[order_now]))<=(svrad)  )
      if(length(repeat_pos)>0){
        if_condition<-length(which(  genoname[posx2[repeat_pos],2]==as.numeric(genoname[posx1[order_now],2])  ))==0
        if(if_condition){
          chose_peak<-c(chose_peak,posx1[order_now])
        }
      }else{
        chose_peak<-c(chose_peak,posx1[order_now])
      }
    }
    return(chose_peak)
  }
  ebayes_EM<-function(x,z,y,en,v0,v,tau,err_max){
    n<-nrow(z);k<-ncol(z)
    mk<-3*en; kn<-k/mk

    v0<-as.numeric(v0)
    v<-matrix(v,ncol=1)
    if(abs(min(eigen(crossprod(x,x))$values))<1e-6){
      try_b<-try({  b<-chol2inv(chol(crossprod(x,x)+diag(ncol(x))*1e-8))%*%crossprod(x,y)  },silent=TRUE)
      if('try-error' %in% class(try_b)){
        try_c<-try({  b<-solve(crossprod(x,x))%*%crossprod(x,y)   },silent=TRUE)
        if('try-error' %in% class(try_c)){   b<-ginv(crossprod(x,x))%*%crossprod(x,y)  }
      }
    }else{
      try_b<-try({  b<-chol2inv(chol(crossprod(x,x)))%*%(crossprod(x,y))  },silent=TRUE)
      if('try-error' %in% class(try_b)){
        try_c<-try({  b<-solve(crossprod(x,x))%*%(crossprod(x,y))   },silent=TRUE)
        if('try-error' %in% class(try_c)){  b<-ginv(crossprod(x,x))%*%(crossprod(x,y))   }
      }
    }# β: fixed effect-rough estimate
    u<-matrix(0,nrow=mk,ncol=kn)# E(γk)
    w<-matrix(0,nrow=mk,ncol=k)# var(γk)
    s<-matrix(0,nrow=kn,ncol=1)# tr(var(γk))
    vv<-matrix(0,n,n)# V
    for(i in 1:kn){
      nc<-( (i-1)*mk+1 ):(i*mk)
      zz<-z[,nc]# Zk

      vv=vv+tcrossprod(zz,zz)*v[i,]# ∑( Zk%*%t(Zk)*(σk^2) )
    }
    vv<-vv+diag(n)*v0# V : the covariance matrix for y

    # genotype effect transform to additive dominance(transformation matrix)
    L_coefficient<-function(en){
      e.seq<-rep(1,3*en)
      a11<-matrix(0,1,3*en)
      a12<-matrix(0,1,3*en)
      a13<-matrix(0,1,3*en)
      a.seq<-seq(1,3*en,by=3)
      a11[a.seq]<-1
      a12[a.seq+1]<-1
      a13[a.seq+2]<-1
      a123<-rbind(a11,a12,a13)

      a1<-1/en*a11-1/(3*en)*e.seq
      a2<-1/en*a12-1/(3*en)*e.seq
      a3<-1/en*a13-1/(3*en)*e.seq

      L1<-rbind(a1,a2,a3)
      a4<-matrix(c(0.5,-0.5,0,1,-0.5,-0.5),2,3)
      LL1<-a4%*%L1

      L2<-matrix(0,en,3*en)
      L3<-matrix(0,3*en,3*en)
      c4<-matrix(0,2*en,3*en)
      for(i in 1:en){

        b11<-matrix(0,1,3*en)
        b11[((i-1)*3+1):(i*3)]<-1
        L2[i,]<-1/3*b11-1/(3*en)*e.seq

        c4[((i-1)*2+1):(i*2),((i-1)*3+1):(i*3)]<-a4

        for(i0 in 1:3){
          seq.c<-(i-1)*3+i0
          c0<-matrix(0,1,3*en)
          c0[seq.c]<-1
          L3[seq.c,]<--1/en*a123[i0,]-1/3*b11+1/(3*en)*e.seq+c0
        }
      }
      LL3<-c4%*%L3
      return(list(matrix_C1=L1, matrix_C2=L2, matrix_C3=L3 , LL1=LL1, LL3=LL3, L1=a4, L3=c4))
    }
    L<-L_coefficient(en)#L coefficient matrix # rm(L_coefficient);gc()
    rank_1<-qr(L$LL1)$rank
    rank_2<-qr(L$LL3)$rank

    iter<-0;err<-1000;iter_max<-500;
    omega<-0
    while( (iter<iter_max)&&(err>err_max) ){
      iter<-iter+1
      v01<-v0# v01 is the initial σ^2
      v1<-v# v1 is the initial σk^2
      b1<-b# b1 is the initial β
      #s1<-s

      try_a<-try({ vi<-chol2inv(chol(vv)) },silent=TRUE)# solve(V)
      if('try-error' %in% class(try_a)){
        try_aa<-try({ vi<-solve(vv) },silent=TRUE)
        if('try-error' %in% class(try_aa)){ vi<-ginv(vv) }
      }

      xtv<-crossprod(x,vi)# t(X)%*%solve(V)

      if(ncol(x)==1){
        b<-((xtv%*%x)^(-1))*(xtv%*%y)
      }else{
        if(abs(min(Mod(eigen(xtv%*%x)$values)))<1e-6){
          try_b<-try({  b<-chol2inv(chol((xtv%*%x)+diag(ncol(x))*1e-8))%*%(xtv%*%y) },silent=TRUE)
          if('try-error' %in% class(try_b)){
            try_c<-try({  b<-solve((xtv%*%x))%*%(xtv%*%y)  },silent=TRUE)
            if('try-error' %in% class(try_c)){  b<-ginv((xtv%*%x))%*%(xtv%*%y)  }
          }
        }else{
          try_b<-try({ b<-chol2inv(chol(xtv%*%x))%*%(xtv%*%y) },silent=TRUE)
          if('try-error' %in% class(try_b)){
            try_c<-try({  b<-solve((xtv%*%x))%*%(xtv%*%y)  },silent=TRUE)
            if('try-error' %in% class(try_c)){  b<-ginv((xtv%*%x))%*%(xtv%*%y)  }
          }
        }
      }
      r<-y-x%*%b# y-Xβ
      ss<-matrix(0,nrow=n,ncol=1)
      vv<-matrix(0,n,n)# new V

      for(i in 1:kn){
        nc<-( (i-1)*mk+1 ):(i*mk)
        zz<-z[,nc]# Zk
        zztvi<-crossprod(zz,vi)# t(Zk)%*%solve(V)
        u[,i]<-v[i,]*zztvi%*%r# E(γk)
        w[,nc]<-v[i,]*( diag(1,mk)-zztvi%*%zz*v[i,] )# var(γk)
        s[i,]<-sum(diag(w[,nc]))# tr(var(γk))
        v[i,]<-(crossprod(u[,i,drop=F],u[,i,drop=F])+s[i,]+omega)/(tau+2+mk)# new (σk^2)
        ss<-ss+zz%*%u[,i,drop=F]
        vv<-vv+tcrossprod(zz,zz)*v[i,]# ∑( Zk%*%t(Zk)*(σk^2) )
      }
      v0<-as.numeric(crossprod(r,(r-ss))/n)# new σ^2
      vv<-vv+diag(n)*v0# new V

      err<-(crossprod((b1-b),(b1-b))+(v01-v0)^2+crossprod((v1-v),(v1-v)))/(1+ncol(x)+kn)
      beta<-t(b)
      sigma2<-v0
    }

    u1<-matrix(0,nrow=2,ncol=kn)# main-E(γk)
    u2<-matrix(0,nrow=2*en,ncol=kn)# interaction-E(γk)
    p1<-matrix(1,kn,1)
    p2<-matrix(1,kn,1)
    # pvalue<-matrix(1,kn,1)

    for(i in 1:kn){
      nc<-( (i-1)*mk+1 ):(i*mk)
      gammak<-u[,i,drop=F]

      u1[,i]<-L$LL1%*%gammak
      u2[,i]<-L$LL3%*%gammak

      var_1<-L$LL1%*%w[,nc]%*%t(L$LL1)
      tr_1<-sum(diag(ginv(tcrossprod(L$LL1))%*%var_1))
      dk1<-abs(rank_1-tr_1/v[i,])

      var_2<-L$LL3%*%w[,nc]%*%t(L$LL3)
      tr_2<-sum(diag(ginv(tcrossprod(L$LL3))%*%var_2))
      dk2<-abs(rank_2-tr_2/v[i,])

      p1[i,]<-1-pchisq(  t(u1[,i,drop=F])%*%ginv(L$LL1%*%w[,nc]%*%t(L$LL1))%*%u1[,i,drop=F],     2)
      p2[i,]<-1-pchisq(  t(u2[,i,drop=F])%*%ginv(L$LL3%*%w[,nc]%*%t(L$LL3))%*%u2[,i,drop=F],   2*(en-1))
    }
    return(list(b=b,u=u,u1=u1,u2=u2,sigma2=sigma2,p1=p1,p2=p2,iter=iter))
  }
  Zhou_sbl<-function(peak,Order,DesignMatrix,CodeMatrix,genoname,en,sbl_t,sbl_p,tau,err_max,fix_p,Sigma,SigmaK){

    chr_n<-as.numeric(genoname[,2])
    chr_n<-chr_n[!duplicated(chr_n)]
    maxdistance<-0
    for(i in 1:length(chr_n)){
      maxdistance<-max(maxdistance,   max(diff(as.matrix(genoname[which(genoname[,2]==i),3]))))
    }

    result_sblgwas<-sblgwas(x=x,y=y,z=DesignMatrix,t=sbl_t)
    sbl_p_wald<-result_sblgwas$blup[,4]# sbl_par<-result_sblgwas$blup[,1]# sbl_p_wald<-p.adjust(sbl_p_wald, method = "bonferroni")
    sbl_pos_order<-Order[order(sbl_p_wald)]
    p_order<-sort(sbl_p_wald)
    id1<-which( p_order< (1-pchisq(sbl_p*2*log(10),1)))
    id2<-which((p_order>=(1-pchisq(sbl_p*2*log(10),1)))&(p_order<1))

    if(length(id1)>0){

      sbl_pos_order1<-sbl_pos_order[id1]
      sbl_pos_order2<-sbl_pos_order[id2]
      sbl_pos_order1<-sbl_pos_order1[!duplicated(sbl_pos_order1)]
      sbl_pos_order2<-sbl_pos_order2[!duplicated(sbl_pos_order2)]

      sort_order1<-sort(sbl_pos_order1)
      result_emba<-ebayes_EM(x,multi_peak_new(CodeMatrix,sort_order1,en)$z,y,en,Sigma,SigmaK[sort_order1],tau,err_max)
      emID1<-which((result_emba$p1<(1-pchisq(fix_p*2*log(10),2)))|(result_emba$p2<(1-pchisq(fix_p*2*log(10),2*(en-1)))))
      emID2<-union(order(result_emba$p1)[seq(1,5,1)],order(result_emba$p2)[seq(1,5,1)])

      if(length(emID1)>5){
        emID<-sort(emID1)
        fix_pos<-sort_order1[emID]

        if(maxdistance<=1){
          sbl_pos_order1<-selection2(sbl_pos_order1,fix_pos,genoname,1)
          sbl_pos_order1<-selection(sbl_pos_order1,genoname,1)
          sbl_pos_order2<-selection2(sbl_pos_order2,fix_pos,genoname,1)
          sbl_pos_order2<-selection2(sbl_pos_order2,sbl_pos_order1,genoname,1)
          sbl_pos_order2<-selection(sbl_pos_order2,genoname,1)
        }
      }else{
        emID<-sort(union(emID1,emID2))
        fix_pos<-sort_order1[emID]

        if(maxdistance<=1){
          sbl_pos_order1<-selection(sbl_pos_order1,genoname,1)
          sbl_pos_order2<-selection2(sbl_pos_order2,sbl_pos_order1,genoname,1)
          sbl_pos_order2<-selection(sbl_pos_order2,genoname,1)
        }
      }

      sbl_pos<-sort(c(sbl_pos_order1,sbl_pos_order2))# length(union(sbl_pos,fix_pos))

    }else{

      result_emba<-ebayes_EM(x,multi_peak_new(CodeMatrix,peak,en)$z,y,en,Sigma,SigmaK[peak],tau,err_max)
      emID1<-which((result_emba$p1<(1-pchisq(fix_p*2*log(10),2)))|(result_emba$p2<(1-pchisq(fix_p*2*log(10),2*(en-1)))))
      emID2<-union(order(result_emba$p1)[seq(1,5,1)],order(result_emba$p2)[seq(1,5,1)])

      if(length(emID1)>5){
        emID<-sort(emID1)
        fix_pos<-peak[emID]

        if(maxdistance<=1){
          sbl_pos_order1<-selection2(peak,fix_pos,genoname,1)
          sbl_pos_order1<-selection(sbl_pos_order1,genoname,1)
        }else{
          sbl_pos_order1<-peak[-emID]
        }

      }else{
        emID<-sort(union(emID1,emID2))
        fix_pos<-peak[emID]

        if(maxdistance<=1){
          sbl_pos_order1<-selection2(peak,fix_pos,genoname,1)
          sbl_pos_order1<-selection(sbl_pos_order1,genoname,1)
        }else{
          sbl_pos_order1<-peak[-emID]
        }
      }

      sbl_pos<-sort(sbl_pos_order1)
    }

    sbl_fix_pos<-sort(fix_pos)
    xin<-cbind(x,multi_peak_new(CodeMatrix,sbl_fix_pos,en)$z)
    return(list(fix=sbl_fix_pos,pos=sbl_pos,xin=xin))
  }
  multi_code_classic<-function(peak,ee,n_id){
    mk<-2+2*(en-1)# the number of genotypes at per locus
    lengthpeak<-length(peak)
    gen_A<-NULL
    gen_D<-NULL
    for(i in 1:en){
      gen_A<-rbind(gen_A,(Ax-Bx)[n_id[[i]],peak,drop=F])
      gen_D<-rbind(gen_D,     Hx[n_id[[i]],peak,drop=F])
    }
    adgen3<-matrix(0,nrow=n,ncol=lengthpeak*mk)

    adgen3[,seq(1,lengthpeak*mk,mk)]<-gen_A
    adgen3[,seq(2,lengthpeak*mk,mk)]<-gen_D

    for(i in 1:lengthpeak){
      col.1<-seq( ((i-1)*mk+3),(i*mk),by=2 )
      col.2<-seq( ((i-1)*mk+4),(i*mk),by=2 )
      adgen3[,col.1]<-gen_A[,i]*ee
      adgen3[,col.2]<-gen_D[,i]*ee
    }
    return(adgen3)
  }
  multinormal<-function(y,mean,sigma){
    pdf_value<-(1/sqrt(2*3.14159265358979323846*sigma))*exp(-(y-mean)*(y-mean)/(2*sigma));
    return (pdf_value)
  }
  LRT_F2<-function(xxn,xxx,yn,par,mk){
    # mk<-2+2*(en-1)# the number of genotypes at per locus
    xn<-ncol(as.matrix(xxn))
    nq<-ncol(xxx)
    ns<-nrow(yn)
    kn<-nq/mk
    at1<-nq

    ad<-if(at1>0.5) cbind(xxn,xxx) else xxn
    if(length(par)==0){
      bb<-if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6) solve(crossprod(ad,ad)+diag(ncol(ad))*1e-8)%*%crossprod(ad,yn) else solve(crossprod(ad,ad))%*%crossprod(ad,yn)
    }else{
      bb<-par
    }
    vv<-as.numeric(crossprod((yn-ad%*%bb),(yn-ad%*%bb))/ns)# y-(X Z)t(β γ)
    ll<-sum(log(abs(multinormal(yn,ad%*%bb,vv))))
    lod<-matrix(0,kn,1)

    if(at1>0.5){
      for(m in 1:kn){
        i1<-(((m-1)*mk+1):(m*mk));# i2<-((m-1)*mk+1):((m-1)*mk+2);  i3<-((m-1)*mk+3):(m*mk)
        m1<-seq(1,ncol(ad),1)[-c(i1+xn)];# m2<-sub[-c(i2+xn)];  m3<-sub[-c(i3+xn)]

        ad1<-ad[,m1,drop=F]
        if(length(par)==0){
          bb1<-if(abs(min(eigen(crossprod(ad1,ad1))$values))<1e-6) solve(crossprod(ad1,ad1)+diag(ncol(ad1))*1e-8)%*%crossprod(ad1,yn) else solve(crossprod(ad1,ad1))%*%crossprod(ad1,yn)
        }else{
          bb1<-par[m1]
        }
        vv1<-as.numeric(crossprod((yn-ad1%*%bb1),(yn-ad1%*%bb1))/ns)
        ll1<-sum(log(abs(multinormal(yn,ad1%*%bb1,vv1))))
        lod[m,]<--2.0*(ll1-ll)/(2.0*log(10))
      }
    }

    return(lod)
  }
  optimize_every_posx<-function(xpos,z,yn,genoname,en,rr,tau,err_max){
    chr_n<-as.numeric(genoname[,2])
    chr_n<-chr_n[!duplicated(chr_n)]
    maxdistance<-0
    for(i in 1:length(chr_n)){
      maxdistance<-max(maxdistance,   max(diff(as.matrix(genoname[which(genoname[,2]==i),3]))))
    }
    if(maxdistance<rr){
      rr<-rr/maxdistance
      ad<-cbind(x,multi_peak_new(z,xpos,en)$z)

      if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6){
        try_a<-try({   bb<- chol2inv(chol(crossprod(ad,ad)+diag(ncol(ad))*1e-8))%*%crossprod(ad,yn)  },silent=TRUE)
        if('try-error' %in% class(try_a)){
          try_aa<-try({bb<- solve(crossprod(ad,ad))%*%crossprod(ad,yn)},silent=TRUE)
          if('try-error' %in% class(try_aa)){  bb<- ginv(crossprod(ad,ad))%*%crossprod(ad,yn)}
        }
      }else{
        try_a<-try({  bb<-chol2inv(chol(crossprod(ad,ad)))%*%crossprod(ad,yn) },silent=TRUE)
        if('try-error' %in% class(try_a)){
          try_aa<-try({bb<- solve(crossprod(ad,ad))%*%crossprod(ad,yn)},silent=TRUE)
          if('try-error' %in% class(try_aa)){  bb<- ginv(crossprod(ad,ad))%*%crossprod(ad,yn)}
        }
      }
      par<-bb[-c(1:dim(x)[2])]
      result_pos<-xpos
      chr_sum<-NULL

      for(i in 1:length(chr_n)){
        chr_sum<-c(chr_sum,length(which(genoname[,2]==i)))
      }
      chr_sum<-c(0,chr_sum)
      for(i in 1:length(xpos)){
        yy<-y-multi_peak_new(z,xpos[-i],en)$z%*%par[-seq((i-1)*3*en+1,i*3*en,1)]
        chr_now<-apply(genoname[,2,drop=F],2,as.numeric)[xpos[i]]
        if(i==1){
          left_rr<-min(xpos[i]-1-sum(chr_sum[seq(1,chr_now,1)]),rr)
        }else{
          if(genoname[xpos[i-1],2]==genoname[xpos[i],2]){
            left_rr<-min(0.5*(xpos[i]-xpos[i-1]),rr)
          }else{
            left_rr<-min(xpos[i]-1-sum(chr_sum[seq(1,chr_now,1)]),rr)
          }
        }
        if(i==length(xpos)){
          right_rr<-min(sum(chr_sum[seq(1,chr_now+1,1)])-xpos[i],rr)
        }else{
          if(genoname[xpos[i+1],2]==genoname[xpos[i],2]){
            right_rr<-min(0.5*(xpos[i+1]-xpos[i]),rr)
          }else{
            right_rr<-min(sum(chr_sum[seq(1,chr_now+1,1)])-xpos[i],rr)
          }
        }
        left_rr<-floor(left_rr)
        right_rr<-floor(right_rr)
        least_pos<-xpos[-i]
        now_pos<-c((xpos[i]-left_rr):(xpos[i]+right_rr))
        try({
          result_emba2<-ebayes_EM(x,multi_peak_new(z,now_pos,en)$z,yy,en,initial_sigma,initial_sigmak[now_pos],tau,err_max)
          maxp1<-min(result_emba2$p1)
          maxp2<-min(result_emba2$p2)
          max_pos1<-now_pos[which.min(result_emba2$p1)]
          max_pos2<-now_pos[which.min(result_emba2$p2)]
          max_pos1
          max_pos2
          if((maxp1!=1)|(maxp2!=1)){
            if(max_pos1==max_pos2){
              result_pos[i]<-max_pos1
            }else{
              result_pos[i]<-c(max_pos1,max_pos2)[which.min(c(maxp1,maxp2))]
            }
          }
        })
      }
    }else{
      result_pos<-xpos
    }

    return(result_pos)
  }
  effect_estimation<-function(n_id,xpos,lod,en,ee){
    xmatrix<-multi_code_classic(xpos,ee,n_id)
    na_id1<-which(lod[,2]==0)
    na_id2<-which(lod[,3]==0)
    mk<-2+2*(en-1)
    na_id_Q<-sort(c((na_id1-1)*mk+1,(na_id1-1)*mk+2))
    na_id_QE<-NULL
    for(jj in 1:length(na_id2)){
      na_id_QE<-c(na_id_QE,sort(rep((na_id2[jj]-1)*mk,2*(en-1))+seq(3,mk,1)))
    }

    xmatrix0<-xmatrix
    xmatrix0<-cbind(x,xmatrix0)
    xmatrix[,c(na_id_Q,na_id_QE)]<-0
    ad<-cbind(x,xmatrix)
    bb0<-if(abs(min(eigen(crossprod(ad,ad))$values))<1e-6) solve(crossprod(ad,ad)+diag(ncol(ad))*1e-8)%*%crossprod(ad,y) else solve(crossprod(ad,ad))%*%crossprod(ad,y)
    bb<-bb0[-seq(1,dim(x)[2],1)]
    sig_e2<-as.numeric(crossprod(y-xmatrix0%*%bb0)/length(y))

    Into_matrix<-function(vector_x,row_n){
      col_n<-length(vector_x)/row_n
      result_x<-matrix(0,nrow=row_n,ncol=col_n)
      for(i in 1:col_n){
        result_x[,i]<-vector_x[((i-1)*row_n+1):(i*row_n)]
      }
      return(result_x)
    }
    effect_all_0<-t(Into_matrix(bb,2+2*(en-1)))

    last_effect<-function(effect){
      a<-effect[,seq(3,2+2*(en-1),2),drop=F]
      d<-effect[,seq(4,2+2*(en-1),2),drop=F]
      last_matrix<-matrix(0,nrow=dim(a)[1],ncol=2)
      for(i in 1:(en-1)){
        last_matrix[,1]<-last_matrix[,1]-a[,i]
        last_matrix[,2]<-last_matrix[,2]-d[,i]
      }
      return(last_matrix)
    }
    effect_all<-cbind(effect_all_0,last_effect(effect_all_0))# dim(effect_all)

    ef_Q<-effect_all[,c(1:2),drop=F]
    ef_QE<-effect_all[,-c(1:2),drop=F]
    sig_Q<-0.5*(ef_Q[,1])^2+0.25*(ef_Q[,2])^2
    sig_QE<-matrix(0,nrow=dim(effect_all)[1],ncol=1)
    for(i in 1:en){
      sig_QE<-sig_QE+(1/en)*0.5*(ef_QE[,1+(i-1)*2])^2+(1/en)*0.5*(ef_QE[,2+(i-1)*2])^2
    }
    sig_y<-max(var(y),(sum(sig_Q)+sum(sig_QE)+sig_e2))
    pve<-cbind((sig_Q/sig_y)*100,(sig_QE/sig_y)*100)
    pve<-cbind(as.matrix(pve[,1]+pve[,2]),pve)

    return(list(effect_all,pve,sig_Q,sig_QE,sig_e2,sig_y))
  }
  LeftRight_marker<-function(map,ChrPos){
    LR_result<-NULL
    for(i in 1:dim(ChrPos)[1]){
      now_id<-which(as.numeric(map[,2])==as.numeric(ChrPos[i,1]))
      now_pos<-as.numeric(ChrPos[i,2])
      all_pos<-as.numeric(map[now_id,3])
      if(now_pos<min(all_pos)){
        left_mar<-""
      }else{
        left_id<-max(which(all_pos<=now_pos))
        left_mar<-map[now_id,1][left_id]
      }
      if(now_pos>max(all_pos)){
        right_mar<-""
      }else{
        right_id<-min(which(all_pos>=now_pos))
        right_mar<-map[now_id,1][right_id]
      }

      LR_result<-rbind(LR_result,c(left_mar,right_mar))
    }
    return(LR_result)
  }
  ######################################### input and basic setup  #########################################
  #*#########  environment and phenotype  #
  en<-EnvNum[NUM]
  sum_en<-sum(EnvNum[0:(NUM-1)])
  pheno<-t(pheRaw[,(sum_en+1):(sum_en+en),drop=F])
  # rownames(pheno)<-NULL
  yes_id<-NULL
  for(i in 1:dim(pheno)[1]){
    yes_id[[i]]<-which(pheno[i,]!="-")
  }
  # pheno<-as.matrix(pheno)
  y<-NULL;yall<-NULL
  for(i in 1:dim(pheno)[1]){
    y<-c(y,as.numeric(pheno[i,yes_id[[i]]]))
    yall<-c(yall,pheno[i,])
  }
  n0<-dim(pheno)[2]# The number of individuals in each environment
  n<-length(y)# The number of individuals in all environments after deleting the missing values
  nn<-length(yall)# The number of individuals in all environments before deleting the missing values
  y<-as.matrix(y) # rm(pheno);gc()

  #*#########  genotype  #
  # genRaw<-as.matrix(genRaw)

  #*#########  calculate Z matrix for K matrix  #
  mn<-dim(Ax0)[2]
  Z<-matrix(0,nrow=n,ncol=mn*en*3)
  sum_n<-0
  for(j in 1:en){
    Z[  (sum_n+1):(sum_n+length(yes_id[[j]])),  seq((j-1)*3+1,mn*en*3,3*en)  ]<-Ax0[yes_id[[j]],]
    Z[  (sum_n+1):(sum_n+length(yes_id[[j]])),  seq((j-1)*3+2,mn*en*3,3*en)  ]<-Hx0[yes_id[[j]],]
    Z[  (sum_n+1):(sum_n+length(yes_id[[j]])),  seq((j-1)*3+3,mn*en*3,3*en)  ]<-Bx0[yes_id[[j]],]
    sum_n<-sum_n+length(yes_id[[j]])
  }# dim(Z)

  #*#########  calculate K matrix  #
  K<-kinship_all(Z,yes_id,en)# rm(kinship_every,kinship_all,K0);gc()

  #*#########  calculate Z matrix  for the subsequent algorithm #
  mn<-dim(Ax)[2]
  Z<-matrix(0,nrow=n,ncol=mn*en*3)
  sum_n<-0
  for(j in 1:en){
    Z[  (sum_n+1):(sum_n+length(yes_id[[j]])),  seq((j-1)*3+1,mn*en*3,3*en)  ]<-Ax[yes_id[[j]],]
    Z[  (sum_n+1):(sum_n+length(yes_id[[j]])),  seq((j-1)*3+2,mn*en*3,3*en)  ]<-Hx[yes_id[[j]],]
    Z[  (sum_n+1):(sum_n+length(yes_id[[j]])),  seq((j-1)*3+3,mn*en*3,3*en)  ]<-Bx[yes_id[[j]],]
    sum_n<-sum_n+length(yes_id[[j]])
  }# dim(Z)

  #*#########  X matrix;     y=Xβ+Zγ+ξ+ε  #
  x<-fixed_x(yes_id,en)
  if(is.null(yygg)==FALSE){
    yygg_x<-NULL
    for(i in 1:en){
      yygg_x<-rbind(yygg_x,yygg[yes_id[[i]],])
    }# dim(yygg_x)
    x<-cbind(x,yygg_x)
  }# dim(x)
  if(det(crossprod(x,x))==0){
    warning("X is singular")
  }
  ReduceDim_x<-TRUE
  if(ReduceDim_x){
    x_effect<-if(abs(min(eigen(crossprod(x,x))$values))<1e-6) solve(crossprod(x,x)+diag(ncol(x))*1e-8)%*%crossprod(x,y) else solve(crossprod(x,x))%*%crossprod(x,y)
    # fix_ef<-bb[seq(1,dim(x)[2],1)]
    yygg_effect<-x_effect[-seq(1,en,1),1,drop=F]
    y<-y-x[,-seq(1,en,1),drop=F]%*%yygg_effect
    x<-fixed_x(yes_id,en)
  }
  #*#########  name  #
  effect_name<-name_function(en)
  ######################################### single_locus_scanning  #########################################
  #*#########  single locus scanning  #
  p3d_result<-p3d_method(x,y,K)

  if(Model=="Random"){
    single_locus_model_result<-single_locus_model(x=x,y=y,zz=Z,lambda=p3d_result$lambda,uu=p3d_result$k_vector,vv=p3d_result$k_value,en=en,CLO=CLO)
    initial_sigma<-mean(single_locus_model_result$sigma_2)
    initial_sigmak<-single_locus_model_result$phi_k# write.table(single_locus_model_result,file=paste("single_locus_model_result.csv",sep = ""),sep=",",row.names = F,col.names = T)
  }else if(Model=="Fixed"){
    single_locus_model_result<-single_locus_model_Fixed(x=x,y=y,zz=Z,lambda=p3d_result$lambda,uu=p3d_result$k_vector,vv=p3d_result$k_value,en=en,CLO=CLO)
    initial_sigma<-mean(single_locus_model_result$sigma_2)
    initial_sigmak<-rep(1,mn)
  }else{
    warning("Please enter Model!")
  }
  #*#########  pick the peaks  #
  peak_id1<-peak_selection(single_locus_model_result$log1,genoname)
  peak_id3<-peak_selection(single_locus_model_result$log3,genoname)
  peak_id<-sort(union(peak_id1,peak_id3))# length(peak_id)
  ######################################### multi_locus_scanning  #########################################
  multi_locus_result1<-Zhou_lars(peak_id,Z,n,en)
  # length(multi_locus_result1$lar_pos)
  multi_locus_result2<-Zhou_sbl(peak=multi_locus_result1$lar_pos,Order=multi_locus_result1$order,
                                DesignMatrix=multi_locus_result1$Matrix,CodeMatrix=Z,
                                genoname=genoname,en=en,
                                sbl_t=-1,sbl_p=3,tau=0,err_max=1e-6,fix_p=1.5,
                                Sigma=initial_sigma,SigmaK=initial_sigmak)# larpos=multi_locus_result1$lar_pos0,larbeta=multi_locus_result1$beta

  emba_p<-1.5
  t1<-proc.time()
  result_emba<-ebayes_EM(multi_locus_result2$xin,multi_peak_new(Z,multi_locus_result2$pos,en)$z,
                         y,en,initial_sigma,initial_sigmak[multi_locus_result2$pos],
                         tau=-2,err_max=1e-6)
  t2<-proc.time()
  (t2-t1)[3]
  emba_pos0<-which((result_emba$p1<(1-pchisq(emba_p*2*log(10),2)))|(result_emba$p2<(1-pchisq(emba_p*2*log(10),2*(en-1)))))# cbind(result_emba$p1,result_emba$p2)
  emba_all_pos<-sort(c(multi_locus_result2$fix,multi_locus_result2$pos[emba_pos0]))

  if(length(multi_locus_result2$pos[emba_pos0])>0){
    emba_pos_Q<-multi_locus_result2$pos[emba_pos0]
    emba_pos_QE<-multi_locus_result2$pos[emba_pos0]
    emba_pos<-multi_locus_result2$pos[emba_pos0]
    if(length(emba_pos_Q)>0){
      z_M4_Q <-multi_code_classic(peak=emba_pos_Q, ee=x[,seq(2,en),drop=F],n_id=yes_id)
      order_Q<-sort(c(seq(1,dim(z_M4_Q)[2],2+2*(en-1)),seq(2,dim(z_M4_Q)[2],2+2*(en-1))))
      z_M4_Q<-z_M4_Q[,order_Q,drop=F]
      lod_Q<-LRT_F2(xxn=multi_locus_result2$xin,xxx=z_M4_Q, yn=y,par=NULL,mk=2)# cbind(emba_pos_Q,lod_Q)
      lrt_pos_Q<-emba_pos_Q[which(lod_Q>2.5)]
      emba_pos_Q[which(lod_Q>2.5)]# cbind(emba2_pos_Q[which(lod_Q>2.5)],lod_Q[which(lod_Q>2.5)])
    }else{
      lrt_pos_Q<-NULL
    }
    if(length(emba_pos_QE)>0){
      z_M4_QE<-multi_code_classic(peak=emba_pos_QE,ee=x[,seq(2,en),drop=F],n_id=yes_id)
      order_QE<-sort(c(seq(1,dim(z_M4_QE)[2],2+2*(en-1)),seq(2,dim(z_M4_QE)[2],2+2*(en-1))))
      z_M4_QE<-z_M4_QE[,-order_QE,drop=F]
      lod_QE<-LRT_F2(xxn=multi_locus_result2$xin,xxx=z_M4_QE, yn=y,par=NULL,mk=2*en-2)# cbind(emba_pos_QE,lod_QE)
      lrt_pos_QE<-emba_pos_QE[which(lod_QE>2.5)]
      emba_pos_QE[which(lod_QE>2.5)]# cbind(emba2_pos_QE[which(lod_QE>2.5)],lod_QE[which(lod_QE>2.5)])
    }else{
      lrt_pos_QE<-NULL
    }
    lrt_pos<-sort(union(lrt_pos_Q,lrt_pos_QE))
    lrt_pos<-sort(union(multi_locus_result2$fix,lrt_pos))# length(lrt_pos)

  }else{
    lrt_pos<-multi_locus_result2$fix# length(lrt_pos)
  }
  ######################################### Optimization and output #########################################
  if(length(lrt_pos)>0){
    optimize_pos<-optimize_every_posx(xpos=lrt_pos,z=Z,yn=y,genoname,en,rr=CriDis,tau=0,err_max=1e-6)

    emba3_p<-3
    result_emba3<-ebayes_EM(x,multi_peak_new(Z,optimize_pos,en)$z,y,en,initial_sigma,initial_sigmak[optimize_pos],tau=0,err_max = 1e-8)
    emba3_pos1<-optimize_pos[which(result_emba3$p1<(1-pchisq(emba3_p*2*log(10),2)))]
    emba3_pos2<-optimize_pos[which(result_emba3$p2<(1-pchisq(emba3_p*2*log(10),2*(en-1))))]
    emba3_pos3<-optimize_pos[which((result_emba3$p1>=(1-pchisq(emba3_p*2*log(10),2)))&(result_emba3$p2>=(1-pchisq(emba3_p*2*log(10),2*(en-1)))))]
    emba3_pos_Q<-sort(union(emba3_pos1,emba3_pos3))
    emba3_pos_QE<-sort(union(emba3_pos2,emba3_pos3))
    emba3_pos<-sort(union(emba3_pos_Q,emba3_pos_QE))

    # CriLOD<-3
    if(length(emba3_pos_Q)>0){
      z_M5_Q<-multi_code_classic(emba3_pos_Q,x[,seq(2,en),drop=F],yes_id)
      order_Q<-sort(c(seq(1,dim(z_M5_Q)[2],2+2*(en-1)),seq(2,dim(z_M5_Q)[2],2+2*(en-1))))
      z_M5_Q<-z_M5_Q[,order_Q]
    }else{
      z_M5_Q<-NULL
    }
    if(length(emba3_pos_QE)>0){
      z_M5_QE<-multi_code_classic(emba3_pos_QE,x[,seq(2,en),drop=F],yes_id)
      order_QE<-sort(c(seq(1,dim(z_M5_QE)[2],2+2*(en-1)),seq(2,dim(z_M5_QE)[2],2+2*(en-1))))
      z_M5_QE<-z_M5_QE[,-order_QE,drop=F]
    }else{
      z_M5_QE<-NULL
    }
    if(length(emba3_pos_Q)>0){
      lod_Q<-LRT_F2(xxn=cbind(x,z_M5_QE),xxx=z_M5_Q, yn=y,par=NULL,mk=2)# cbind(emba3_pos_Q,lod_Q)
      lrt_pos_Q<-emba3_pos_Q[which(lod_Q>=CriLOD)]
    }else{
      lrt_pos_Q<-NULL
    }
    if(length(emba3_pos_QE)>0){
      lod_QE<-LRT_F2(xxn=cbind(x,z_M5_Q),xxx=z_M5_QE, yn=y,par=NULL,mk=2*en-2)# cbind(emba3_pos_QE,lod_QE)
      lrt_pos_QE<-emba3_pos_QE[which(lod_QE>=CriLOD)]
    }else{
      lrt_pos_QE<-NULL
    }

    lrt2_pos<-sort(union(lrt_pos_Q,lrt_pos_QE))
    if(length(lrt2_pos)>0){
      last_lod<-matrix(0,nrow=length(lrt2_pos),ncol=3)
      last_lod[which(lrt2_pos%in%lrt_pos_Q),2]<-lod_Q[which(lod_Q>=CriLOD)]
      last_lod[which(lrt2_pos%in%lrt_pos_QE),3]<-lod_QE[which(lod_QE>=CriLOD)]
      last_lod[,1]<-last_lod[,2]+last_lod[,3]

      IC_data<-cbind(x,multi_peak_new(Z,lrt2_pos,en)$z)
      lm_IC<-lm(y~IC_data-1)
      AIC(lm_IC)
      BIC(lm_IC)

      if(length(lrt_pos_Q)>0){
        zM_Q<-multi_code_classic(lrt_pos_Q,x[,seq(2,en),drop=F],yes_id)
        order_Q<-sort(c(seq(1,dim(zM_Q)[2],2+2*(en-1)),seq(2,dim(zM_Q)[2],2+2*(en-1))))
        zM_Q<-zM_Q[,order_Q,drop=F]
      }else{
        zM_Q<-NULL
      }
      if(length(lrt_pos_QE)>0){
        zM_QE<-multi_code_classic(lrt_pos_QE,x[,seq(2,en),drop=F],yes_id)
        order_QE<-sort(c(seq(1,dim(zM_QE)[2],2+2*(en-1)),seq(2,dim(zM_QE)[2],2+2*(en-1))))
        zM_QE<-zM_QE[,-order_QE,drop=F]
      }else{
        zM_QE<-NULL
      }
      IC_data<-cbind(cbind(x,zM_Q),zM_QE)
      lm_IC<-lm(y~IC_data-1)
      AIC(lm_IC)
      BIC(lm_IC)

      LR_marker<-LeftRight_marker(map=mapRaw,ChrPos=genoname[lrt2_pos,2:3,drop=F])
      result_all<-effect_estimation(n_id=yes_id,xpos=lrt2_pos,lod=last_lod,en,ee=x[,seq(2,en),drop=F])
      var_e<-matrix("",nrow=length(lrt2_pos),ncol=1)
      var_y<-matrix("",nrow=length(lrt2_pos),ncol=1)
      var_e[1]<-round(result_all[[5]],4)
      var_y[1]<-round(result_all[[6]],4)

      data.all<-data.frame(genoname[lrt2_pos,2:3,drop=F],
                           round(result_all[[1]],4),round(last_lod,4),
                           LR_marker,
                           round(result_all[[3]]+result_all[[4]],4),
                           round(result_all[[3]],4),round(result_all[[4]],4),
                           round(result_all[[2]],4),
                           var_e,var_y)
      # rep(AIC(lm_IC),length(lrt2_pos)),
      # rep(BIC(lm_IC),length(lrt2_pos)))
      rownames(data.all)<-NULL
      colnames(data.all)<-c("Chr","Position(cM)",effect_name,
                            "LOD","LOD_QTL","LOD_QEI",
                            "Left_marker","right_marker",
                            "Var_Genet","Var_Genet_QTL","Var_Genet_QEI",
                            "r2(%)","r2_QTL(%)","r2_QEI(%)",
                            "Var_Error","Var_Phen(total)")
      reslt_list<-list(result=data.all,p_Q=single_locus_model_result$log1,p_QE=single_locus_model_result$log3)
    }else{
      reslt_list<-NULL
      warning("No QTL or QEI were detected!")
    }
  }else{
    reslt_list<-NULL
    warning("No QTL or QEI were detected!")
  }

  return(reslt_list)
}
