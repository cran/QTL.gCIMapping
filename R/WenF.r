WenF<-function(pheRaw=NULL,genRaw=NULL,mapRaw1=NULL,yygg1=NULL,cov_en=NULL,
              WalkSpeed=NULL,CriLOD=NULL,dir=NULL){
  
  cl<-WalkSpeed;sLOD<-CriLOD;yygg<-NULL
  mx=NULL;phe=NULL;chr_name=NULL;v.map=NULL
  a.gen.orig=NULL;d.gen.orig=NULL;g1=NULL;X.ad.t4=NULL
  
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
  
  if(is.null(yygg1)==FALSE){
    cov_en<-as.matrix(cov_en)
    yygg1<-as.matrix(yygg1)
    covname<-cov_en[2:nrow(cov_en),1]
    yygg1<-cbind(covname,yygg1)
    mapRaw10<-as.matrix(mapRaw1[-1,])
    chr_name<-unique(mapRaw10[,2])
    chr_secon<-as.matrix(mapRaw10[,2])
    mm<-numeric()
    map_chr<-numeric()
    for(i in 1:length(chr_name)){
      chr_i<-which(chr_secon[]==chr_name[i])
      len<-matrix(length(chr_i),,1)
      mm<-rbind(mm,len)
      chr_name[i]<-i
    }
    map_chr<-numeric()
    for(i in 1:length(chr_name)){
      chr_name<-as.numeric(chr_name)
      chr_pos1<-matrix(chr_name[i],mm[i],1)
      map_chr<-rbind(map_chr,chr_pos1)
    }
    map_marker<-matrix(mapRaw10[,1],,1)
    map_pos<-matrix(mapRaw10[,3],,1)
    mapRaw<-cbind(map_marker,map_chr,map_pos)
    nameMap<-matrix(mapRaw[,1],,1)
    nameGenrow<-matrix(genRaw[1,],1,)#individual's name in genotype 
    nameGencol<-matrix(genRaw[,1],,1)#marker's name in genotype
    namePhe<-as.matrix(pheRaw[,1],,1)
    nameCov<-matrix(yygg1[,1],,1)
    if(nameGenrow[2]=="1"){
      phee<-as.matrix(pheRaw[-1,-1])
      phe<-matrix(as.numeric(phee),nrow(phee),ncol(phee))
      mapname<-mapRaw
      genn<-genRaw[-1,-1]
      genn<-as.matrix(genn)
      mapnametwo<-matrix(as.numeric(mapname[,2:3]),nrow(mapname),2)
      mx<-cbind(mapnametwo,genn)
      chrRaw_name<-unique(mapRaw10[,2])
      newyygg<-as.matrix(yygg1)
      yyggChar<-as.matrix(newyygg[,-1])
      yygg<-matrix(as.numeric(yyggChar),nrow(yyggChar),ncol(yyggChar))
      mapname<-mapname
      chrRaw_name<-chrRaw_name
      chr_name<-chr_name
    }else{
      sameName_MG<-intersect(nameMap,nameGencol)
      sameName_PG<-intersect(namePhe,nameGenrow)
      locPhe<-match(sameName_PG,namePhe)
      locMap<-match(sameName_MG,nameMap)
      locGen_PG<-match(sameName_PG,nameGenrow)
      locGen_MG<-match(sameName_MG,nameGencol)
      locCov<-match(sameName_PG,nameCov)
      newyygg<-as.matrix(yygg1[locCov,])
      yyggChar<-as.matrix(newyygg[,-1])
      yygg<-matrix(as.numeric(yyggChar),nrow(yyggChar),ncol(yyggChar))
      newPhe<-as.matrix(pheRaw[locPhe,])
      newMap<-as.matrix(mapRaw[locMap,])
      newGenrow<-as.matrix(genRaw[,locGen_PG])
      newGen<-as.matrix(newGenrow[locGen_MG,])
      gen_two<-newMap[,2:3]
      genChar<-cbind(gen_two,newGen)
      if(ncol(newPhe)==2)
      {
        pheChar<-as.matrix(newPhe[,2:ncol(newPhe)])
      }else if(ncol(newPhe)>2){
        pheChar<-newPhe[,2:ncol(newPhe)]
      }
      chrRaw_name<-unique(mapRaw10[,2])
      chr_name<-chr_name
      mx<-as.matrix(genChar)
      phe<-matrix(as.numeric(pheChar),nrow(pheChar),ncol(pheChar))
      mapname<-newMap
    }
  } else{
    mapRaw10<-as.matrix(mapRaw1[-1,])
    chr_name<-unique(mapRaw10[,2])
    chr_secon<-as.matrix(mapRaw10[,2])
    mm<-numeric()
    map_chr<-numeric()
    for(i in 1:length(chr_name)){
      chr_i<-which(chr_secon[]==chr_name[i])
      len<-matrix(length(chr_i),,1)
      mm<-rbind(mm,len)
      chr_name[i]<-i
    }
    map_chr<-numeric()
    for(i in 1:length(chr_name)){
      chr_name<-as.numeric(chr_name)
      chr_pos1<-matrix(chr_name[i],mm[i],1)
      map_chr<-rbind(map_chr,chr_pos1)
    }
    map_marker<-matrix(mapRaw10[,1],,1)
    map_pos<-matrix(mapRaw10[,3],,1)
    mapRaw<-cbind(map_marker,map_chr,map_pos)
    nameMap<-matrix(mapRaw[,1],,1)
    nameGenrow<-matrix(genRaw[1,],1,)
    nameGencol<-matrix(genRaw[,1],,1)
    namePhe<-as.matrix(pheRaw[,1],,1)
    
    if(nameGenrow[2]=="1"){
      phee<-as.matrix(pheRaw[-1,-1])
      phe<-matrix(as.numeric(phee),nrow(phee),ncol(phee))
      mapname<-mapRaw
      genn<-genRaw[-1,-1]
      genn<-as.matrix(genn)
      mapnametwo<-matrix(as.numeric(mapname[,2:3]),nrow(mapname),2)
      mx<-cbind(mapnametwo,genn)
      chrRaw_name<-unique(mapRaw10[,2])
      chr_name<-chr_name
      mapname<-mapname
      chrRaw_name<-chrRaw_name
    }else{
      sameName_MG<-intersect(nameMap,nameGencol)
      sameName_PG<-intersect(namePhe,nameGenrow)
      locPhe<-match(sameName_PG,namePhe)
      locMap<-match(sameName_MG,nameMap)
      locGen_PG<-match(sameName_PG,nameGenrow)
      locGen_MG<-match(sameName_MG,nameGencol)
      newPhe<-as.matrix(pheRaw[locPhe,])
      newMap<-as.matrix(mapRaw[locMap,])
      newGenrow<-as.matrix(genRaw[,locGen_PG])
      newGen<-as.matrix(newGenrow[locGen_MG,])
      gen_two<-newMap[,2:3]
      genChar<-cbind(gen_two,newGen)
      if(ncol(newPhe)==2)
      {
        pheChar<-as.matrix(newPhe[,2:ncol(newPhe)])
      }else if(ncol(newPhe)>2){
        pheChar<-newPhe[,2:ncol(newPhe)]
      }
      chrRaw_name<-unique(mapRaw10[,2])
      chr_name<-chr_name
      mx<-as.matrix(genChar)
      phe<-matrix(as.numeric(pheChar),nrow(pheChar),ncol(pheChar))
      mapname<-newMap
    }
  }
  mapp1<-as.numeric(mapname[,2:3])
  mapp1<-matrix(mapp1,,2)
  chr<-length(unique(mapp1[,1]))
  for(i in 1:chr){
    pos1<-as.matrix(mapp1[which(mapp1[,1]==i),])
    delerow<-which(duplicated(pos1[,2]))
    if(length(delerow)!=0){
      break
    }else{
      mapname<-mapname
    }
  }
  if(length(delerow)!=0){
    warning("Please correct the map file to make sure whether all the marker positions are different!")
  }else{
    
    phe1<-as.matrix(phe[,1])
    mx<-as.matrix(mx)
    
    g1<-cbind(as.matrix(mapname[,1]),mx)
    blank<-as.matrix(c("phe","",""))
    p1<-t(rbind(blank,phe1))
    pgcombine<-rbind(p1,g1)
    
    # suppressWarnings(dir.create(path="temp",recursive = T))
    # dir1<-"temp/"
    
    write.table(pgcombine,file=paste(dir,"/listeria_rotY",".csv",sep=""),sep=",",row.names = F,col.names = F)
    
    Genotypes=c("A","H","B","D","C")
    Crosstype="f2"
    data.qtl<-read.cross("csvr",dir,"listeria_rotY.csv",genotypes=Genotypes,na.strings = "-",crosstype=Crosstype)
    data.qtl<-jittermap(data.qtl)
    
    bin<-cl  #user's options 1cM
    data.qtl.2<-calc.genoprob(data.qtl,step=bin,error.prob = 0.001)
    rownames(g1)<-NULL
    map.raw<-as.data.frame(g1[,1:3],stringsAsFactors = F)
    map.raw[,2:3]<-sapply(map.raw[,2:3],as.numeric)
    names(map.raw)<-c("marker","chr","pos")
    m0<-nrow(g1)
    
    gen.raw1<-cbind(mapRaw1,genRaw[,2:ncol(genRaw)])
    gen.raw<-gen.raw1[-1,]
    colnames(gen.raw)<-c("id","","",gen.raw1[1,4:ncol(gen.raw)])
    nchr<-max(as.numeric(mx[,1]))
    
    gen3<-calc.genoprob(data.qtl,step=0,error.prob = 0.001)
    marker.aa<-NULL
    marker.dd<-NULL
    for(rr in 1:nchr){
      mapgen<-gen3$geno[[rr]]$prob
      aam.gen<-mapgen[,,1]-mapgen[,,3]
      marker.aa<-cbind(marker.aa,aam.gen)
      ddm.gen<-mapgen[,,2]
      marker.dd<-cbind(marker.dd,ddm.gen)
    }
    a.gen.orig<-t(marker.aa)
    d.gen.orig<-t(marker.dd)
    
    aa.0<-NULL
    dd.0<-NULL
    v.map<-NULL
    nchr<-max(as.numeric(mx[,1]))
    for(ii in 1:nchr){
      map.gen<-data.qtl.2$geno[[ii]]$prob
      aa.gen<-map.gen[,,1]-map.gen[,,3]
      dd.gen<-map.gen[,,2]
      aa.0<-cbind(aa.0,aa.gen)
      dd.0<-cbind(dd.0,dd.gen)
      pos<-attr(map.gen,"map")
      pos.gen.1<-data.frame(pos)
      pos.gen.2<-data.frame(marker=row.names(pos.gen.1),chr=ii,pos=pos.gen.1,locus=1:dim(pos.gen.1)[1])
      id.insert<-grep(pattern = "loc",x=c(as.character(pos.gen.2$marker)))
      insertflag<-rep(1,dim(pos.gen.2)[1])
      insertflag[id.insert]<-0
      pos.gen.3<-data.frame(pos.gen.2,insertflag=insertflag,leftmarker=pos.gen.2$marker,rightmarker=pos.gen.2$marker)
      id.left<-findInterval(pos.gen.3[id.insert,]$pos,map.raw[map.raw$chr==ii,]$pos)
      namesleft<-map.raw[map.raw$chr==ii,]$marker[id.left]
      namesright<-map.raw[map.raw$chr==ii,]$marker[id.left+1]
      pos.gen.3[id.insert,]$leftmarker<-namesleft
      pos.gen.3[id.insert,]$rightmarker<-namesright
      pos.gen.3[id.insert,]$marker<-namesleft
      v.map<-rbind(v.map,pos.gen.3)
    }
    
    ##########################
    names.insert<-v.map
    ######################
    n<-ncol(g1)-3
    m.a<-dim(v.map)[1]
    #############################
    ad<-matrix(NA,nrow = (2*m.a),ncol=n)
    ad[seq(1,(2*m.a),by=2),]<-t(aa.0)
    ad[seq(2,(2*m.a),by=2),]<-t(dd.0)
    ad.gen.insert<-ad
    names.insert1<-names.insert[rep(1:nrow(names.insert),each=2),]
    names.insert2<-cbind(names.insert1,id.all=1:(2*m.a))
    X.ad.tran.data<-cbind(names.insert2,ad.gen.insert)
    X.ad.t4<-t(ad)
    
  }
   output<-list(yygg=yygg,mx=mx,phe=phe,chr_name=chr_name,v.map=v.map,gen.raw=gen.raw,a.gen.orig=a.gen.orig,d.gen.orig=d.gen.orig,n=n,
                names.insert2=names.insert2,X.ad.tran.data=X.ad.tran.data,X.ad.t4=X.ad.t4)
   return(output)
}



WenS<-function(flag=NULL,CriLOD=NULL,NUM=NULL,pheRaw=NULL,Likelihood=NULL,
               flagrqtl=NULL,yygg=NULL,mx=NULL,phe=NULL,chr_name=NULL,v.map=NULL,
               gen.raw=NULL,a.gen.orig=NULL,d.gen.orig=NULL,n=NULL,
               names.insert2=NULL,X.ad.tran.data=NULL,X.ad.t4=NULL,dir=NULL){
  
  sLOD<-CriLOD;result=NULL;mxmp=NULL;galaxyy1=NULL;res1a=NULL;res1d=NULL
  
  mxmp<-mx[,1:2]
  rownames(mxmp)<-NULL
  mxmp<-as.data.frame(mxmp,stringsAsFactors = F)
  mxmp[,1:2]<-sapply(mxmp[,1:2],as.numeric)

  m.a<-dim(v.map)[1]
  Y.orig<-as.matrix(phe[,NUM])
  Y.name<-pheRaw[1,NUM+1]
  ######################################
  
  W.orig<-matrix(1,n,1) 
  id.ind<-which(!is.na(Y.orig))
  Y.part.1<-Y.orig[id.ind,,drop=F]
  
  if(is.null(yygg)==FALSE){
    xenvir<-cbind(matrix(1,n,1),yygg) 
    xenvir<-xenvir[id.ind,]
    beta<-solve(t(xenvir)%*%xenvir)%*%t(xenvir)%*%Y.part.1
    Y.part.1<-Y.part.1-xenvir%*%beta+W.orig[id.ind,]
  }
  bb.adalasso<-NULL
  #kinship matrix#################
  a.K<-kinship.F2(t(a.gen.orig[,id.ind]))
  d.K<-kinship.F2(t(d.gen.orig[,id.ind]))
  kinship<-list(a.K,d.K)
  data.wy<-cbind(Y.part.1,rep(1,times=length(Y.part.1)))
  vars<-suppressWarnings(mixed.vars(data.wy,kinship,optim.speed=TRUE))
  delta.aK.init<-vars$tau.kk[1]/vars$tau.kk[3]
  delta.dK.init<-vars$tau.kk[2]/vars$tau.kk[3]
  
  remle.B<-IRMMA.aK.dK.effects.B(Z=NULL,a.K,d.K,delta.aK.init,delta.dK.init,complete=TRUE)
  C2<-remle.B$C.g
  if(flag==1){
    if (Likelihood=="REML"){
      REML.LRT.c2<-emma.REML.LRT.c.noalpha(ys=Y.orig[id.ind,,drop=F], xs=X.ad.t4[id.ind,], K=1, Z=C2, X0=W.orig[id.ind,,drop=F], ngrids=100, llim=-10, ulim=10, esp=1e-10)
    }else{
      REML.LRT.c2<-emma.ML.LRT.c.noalpha(ys=Y.orig[id.ind,,drop=F], xs=X.ad.t4[id.ind,], K=1, Z=C2, X0=W.orig[id.ind,,drop=F], ngrids=100, llim=-10, ulim=10, esp=1e-10)
    }
  }else{
    if(Likelihood=="REML"){
      REML.LRT.c2<-fixed.REML.LRT.c.sim(ys=Y.orig[id.ind,,drop=F], xs=X.ad.t4[id.ind,], Z=C2, X0=W.orig[id.ind,,drop=F])
    }else{
      REML.LRT.c2<-fixed.ML.LRT.c.sim(ys=Y.orig[id.ind,,drop=F], xs=X.ad.t4[id.ind,], Z=C2, X0=W.orig[id.ind,,drop=F])
    }
  }
  REML.LRT.c3<-data.frame(REML.LRT.c2)
  id.odd<-which(REML.LRT.c3$ID%%2==1)
  id.even<-which(REML.LRT.c3$ID%%2==0)
  LOD.a<--log10(REML.LRT.c3[id.odd,]$ps)
  LOD.d<--log10(REML.LRT.c3[id.even,]$ps)
  #########################3
  optids.a<-peak.id(LOD.a)
  optids.d<-peak.id(LOD.d)
  
  res1a<-cbind(v.map[,2:3],as.matrix(LOD.a))
  res1d<-cbind(v.map[,2:3],as.matrix(LOD.d))
  ###############################################
  ################################33
  a.eb.loc.1<-optids.a*2-1
  a.eb.loc.2<-c(a.eb.loc.1,a.eb.loc.1+1)
  d.eb.loc.1<-optids.d*2
  d.eb.loc.2<-c(d.eb.loc.1,d.eb.loc.1-1)
  ad.eb.loc.0<-union(a.eb.loc.2,d.eb.loc.2)
  ad.eb.loc.0<-ad.eb.loc.0[order(ad.eb.loc.0)]#xdata,a+d
  #########################################3
  xdata.ad<-X.ad.t4[id.ind,ad.eb.loc.0]
  X.ad.eb.names.0<-X.ad.tran.data[ad.eb.loc.0,1:8]
  
  set.seed(11001)  #user's options
  model1<-adalasso(X=xdata.ad,y=Y.part.1,k=10,use.Gram = T,both=T)
  adalasso.beta<-model1$coefficients.adalasso
  
  bb<-matrix(0,nrow=2*m.a,ncol=1)
  bb[ad.eb.loc.0,]<-adalasso.beta
  
  lrt.0<-likelihood.a.d.F2(xxn=W.orig[id.ind,,drop=F],xxx=xdata.ad,yn=Y.part.1,bbo=bb[ad.eb.loc.0,],intercept=model1$intercept.adalasso)
  id.lrt.0.a<-which(lrt.0$lod>=sLOD) #
  
  if(length(id.lrt.0.a)==0){
    warning("No QTL were detected!")
  }else if(length(id.lrt.0.a)>0){
    
  
  id.lrt.0.d<-id.lrt.0.a+1
  id.lrt.0.ad<-c(id.lrt.0.a,id.lrt.0.d)
  id.lrt.0.ad<-id.lrt.0.ad[order(id.lrt.0.ad)]
  
  opt.adalasso<-data.frame(X.ad.eb.names.0[id.lrt.0.a,],effect.a=bb[ad.eb.loc.0,][id.lrt.0.a],effect.d=bb[ad.eb.loc.0,][id.lrt.0.d],LOD=lrt.0$lod[id.lrt.0.a],ps=lrt.0$ps[id.lrt.0.a])
  
  ##########################
  if(flagrqtl=="TRUE"){
    
    # suppressWarnings(dir.create(path="temp",recursive = T))
    # dir1<-"temp/"
    
    b.qtl.1<-bb[ad.eb.loc.0,][id.lrt.0.ad]
    b.id.qtl<-which(b.qtl.1!=0)
    X.qtl<-X.ad.t4[id.ind,ad.eb.loc.0][,id.lrt.0.ad][,b.id.qtl]
    XY.qtl<-data.frame(Y.part.1,X.qtl)
    model.qtl.lm<-lm(Y.part.1~X.qtl,data = XY.qtl)
    
    b.qtl.2<-matrix(0,nrow = length(b.qtl.1),ncol=1)
    b.qtl.2[b.id.qtl]<-model.qtl.lm$coefficients[-1]
    
    ######################
    
    dif.1<-diff(opt.adalasso$locus)
    dif.2<-diff(opt.adalasso$pos)
    
    pos.options<-5 
    locus.options<-6 
    
    
    dif.id<-which(dif.1>0&dif.1<=locus.options&dif.2>0&dif.2<=pos.options)
    dif.id.ad.0<-c(dif.id,dif.id+1)
    
    dif.id.ad.0<-dif.id.ad.0[order(dif.id.ad.0)]
    qtl.id<-opt.adalasso[dif.id.ad.0,]
    dif.id.n<-length(dif.id.ad.0)
    dif.id.ad.1<-c()
    if(dif.id.n!=0){
      opt.qtl.all<-data.frame()
      for(ii in seq(1,dif.id.n,by=2)){
        if(((qtl.id[ii,]$effect.a==0&qtl.id[ii,]$effect.d!=0)&(qtl.id[ii+1,]$effect.a!=0&qtl.id[ii+1,]$effect.d==0))
           |((qtl.id[ii,]$effect.a!=0&qtl.id[ii,]$effect.d==0)&(qtl.id[ii+1,]$effect.a==0&qtl.id[ii+1,]$effect.d!=0))
        ){
          dif.id.ad.1<-c(dif.id.ad.1,dif.id.ad.0[c(ii,ii+1)])
          
          b.qtl.3<-b.qtl.2
          id.omit<-c(2*c(dif.id.ad.0[c(ii,ii+1)])-1,2*c(dif.id.ad.0[c(ii,ii+1)]))
          b.qtl.3[id.omit]<-0
          Y.part.2<-Y.part.1-X.ad.t4[id.ind,ad.eb.loc.0][,id.lrt.0.ad]%*%b.qtl.3
          Y.part.3<-matrix(NA,nrow=n,ncol = 1)
          Y.part.3[id.ind,]<-Y.part.2
          
          im.phe.rot<-t(cbind(Y.part.3,1:n))
          row.names(im.phe.rot)<-c(Y.name,"id")
          
          write.table(im.phe.rot,file=paste(dir,"/",NUM,"_phe_rot_",(ii+1)/2,".csv",sep=""),sep=",",row.names = T,col.names = F)
          ################
          marker.options<-5
          gen.qtl.id.left<-which(gen.raw[,1]==qtl.id[ii,]$leftmarker)
          gen.qtl.id.right<-which(gen.raw[,1]==qtl.id[ii+1,]$rightmarker)
          gen.qtl<-gen.raw[c((gen.qtl.id.left-marker.options):(gen.qtl.id.right+marker.options)),]
          
          im.gen.rot<-rbind(c("id","","",1:n),gen.qtl)
          write.table(im.gen.rot,file=paste(dir,"/",NUM,"_gen_rot_",(ii+1)/2,".csv",sep=""),sep=",",row.names = F,col.names = F)
          
          #####################################
          #rqtl
          data.qtl<-read.cross("csvsr",dir,paste(NUM,"_gen_rot_",(ii+1)/2,".csv",sep=""),paste(NUM,"_phe_rot_",(ii+1)/2,".csv",sep=""))
          data.qtl.1<-jittermap(data.qtl)
          data.qtl.2<-calc.genoprob(data.qtl.1,step=1,error.prob = 0.001)
          
          out.em<-suppressWarnings(scanone(data.qtl.2,method="em"))
          
          out1<-summary(out.em,threshold = 2.5)
          
          data.qtl.3<-sim.geno(data.qtl.1,step=1,n.draws = 1,error.prob =0.001 )
          qtl<-makeqtl(data.qtl.3,chr =out1$chr,pos=out1$pos )
          out2<-suppressWarnings(fitqtl(data.qtl.3,qtl=qtl,formula = y~Q1,dropone=F,get.ests = T))
          
          out3<-summary(out2)
          row.names(out3$ests)[2]
          
          opt.qtl<-data.frame(id.1=row.names(out1),id.2=row.names(out3$ests)[2],out1[,1:2],effect.a=out3$ests[2],effect.d=out3$ests[3],LOD=out1[,3],ps=out3$result.full[1,6],stringsAsFactors = F)
          
          if(opt.qtl$id.1%in%gen.raw[,1]){
            
            names.qtl.1<-names.insert2[which(names.insert2$chr==as.integer(as.character(opt.qtl$chr))),]
            
            names.left<-findInterval(opt.qtl$pos,names.qtl.1$pos)
            names.qtl.2<-names.qtl.1[names.left-1,]
            opt.qtl.lr<-data.frame(names.qtl.2,effect.a=opt.qtl$effect.a,effect.d=opt.qtl$effect.d,LOD=opt.qtl$LOD,ps=opt.qtl$ps,stringsAsFactors=F)
            
          }else{
            
            gen.qtl.1<-gen.qtl[which(gen.qtl[,2]==as.integer(as.character(opt.qtl$chr))),]
            left1<-findInterval(opt.qtl$pos,gen.qtl.1[,3])
            
            #################
            names.qtl.1<-names.insert2[which(names.insert2$chr==as.integer(as.character(opt.qtl$chr))),]
            names.left<-findInterval(opt.qtl$pos,names.qtl.1$pos)
            names.qtl.2<-names.qtl.1[names.left,]
            opt.qtl.lr<-data.frame(names.qtl.2[,1:2],locus=NA,pos=opt.qtl$pos,insertflag=0,leftmarker=names.qtl.2$leftmarker,rightmarker=gen.qtl.1[left1+1,1],id.all=names.qtl.2$id.all-1,effect.a=opt.qtl$effect.a,effect.d=opt.qtl$effect.d,LOD=opt.qtl$LOD,ps=opt.qtl$ps,stringsAsFactors=F)
            
          }
          opt.qtl.all<-rbind(opt.qtl.all,opt.qtl.lr)
        }
      }
      opt.im.all<-rbind(opt.adalasso[-dif.id.ad.1,],opt.qtl.all)
      opt.im.all.0<-opt.im.all[order(opt.im.all$chr,opt.im.all$pos),]
      
    }else{
      opt.im.all.0<-opt.adalasso
    }#if (length(dif.id.ad.0)!=0) end
    
  }else{
    opt.im.all.0<-opt.adalasso
  }
  ###################################
  X.r.id<-c(opt.im.all.0$id.all,opt.im.all.0$id.all+1)
  X.r.id<-X.r.id[order(X.r.id)]
  
  X.r<-X.ad.t4[id.ind,X.r.id]
  num.b<-length(X.r.id)
  b.r<-matrix(0,nrow =num.b)
  b.r[seq(1,num.b,by=2),]<-opt.im.all.0$effect.a
  b.r[seq(2,num.b,by=2),]<-opt.im.all.0$effect.d
  
  X.intc<-matrix(1,nrow=length(id.ind),ncol=1)
  intc<-model1$intercept.adalasso
  
  sigma.e2<-as.numeric(crossprod(Y.part.1-X.r%*%b.r-X.intc%*%intc,Y.part.1-X.r%*%b.r-X.intc%*%intc)/length(id.ind))
  
  r2.g<-as.numeric(sum(0.5*opt.im.all.0$effect.a^2+0.25*opt.im.all.0$effect.d^2))
  ge.all<-max(var(Y.part.1),(r2.g+sigma.e2))
  r2.p<-r2.g/ge.all
  r2.new<-(0.5*opt.im.all.0$effect.a^2+0.25*opt.im.all.0$effect.d^2)/ge.all
  
  Genei<-0.5*opt.im.all.0$effect.a^2+0.25*opt.im.all.0$effect.d^2
  vare<- c(round(sigma.e2,4),matrix("",length(Genei)-1,1))
  varp<- c(round(var(Y.part.1),4),matrix("",length(Genei)-1,1))
  
  #######################
  xdata.opt.adalasso<-data.frame(opt.im.all.0,r2=r2.new*100,PVE=r2.p*100,Genei,vare,varp)
  colnames(xdata.opt.adalasso)[c(13,14)]<-c("r^2(%)","PVE(%)")
  aa.adalasso<-data.frame(nametrait=rep(Y.name,times=dim(opt.im.all.0)[1]),xdata.opt.adalasso,stringsAsFactors=F)
  bb.adalasso<-aa.adalasso[,-c(2,5,6,9,13,15)]

  result<-bb.adalasso[,c(1,2,3,6,7,8,4,5,10,9,11,12)]
  rownames(result)<-NULL
  colnames(result)<-c("Trait","Chr","Position(cM)","Effect.a","Effect.d","LOD","Left_marker","right_marker","Var_Genet","r2(%)","Var_Error","Var_Phen(total)")
 
  galaxyy1<-as.matrix(result[,c(2,3,6)])
  # unlink(dir()[dir()=="temp"],recursive = T)
  }
  output<-list(result=result,mxmp=mxmp,galaxyy1=galaxyy1,res1a=res1a,res1d=res1d,chr_name=chr_name)
  return(output)
}









    
    