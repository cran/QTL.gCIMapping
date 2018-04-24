WangF<-function(pheRaw=NULL,genRaw=NULL,mapRaw1=NULL,yygg1=NULL,
             cov_en=NULL,Population=NULL,WalkSpeed=NULL,CriLOD=NULL,dir=NULL){
  
  cl<-WalkSpeed;sLOD<-CriLOD;yygg<-NULL
  mx=NULL;phe=NULL;chr_name=NULL;gen=NULL;v.map=NULL
  
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
    
    if(Population=="BC1"){
      Genotypes=c("A","H")
      Crosstype="bc"
    }else if(Population=="BC2"){
      Genotypes=c("H","B")
      Crosstype="bc"
    }else if(Population=="DH"){
      Genotypes=c("A","B")
      Crosstype="dh"
    }else if(Population=="RIL"){
      Genotypes=c("A","B")
      Crosstype="riself"
    }
    data.qtl<-read.cross("csvr",dir,"listeria_rotY.csv",genotypes=Genotypes,na.strings = "-",crosstype=Crosstype)
    bin<-cl  #user's options 1cM
    data.qtl.2<-calc.genoprob(data.qtl,step=bin,error.prob = 0.001)
    rownames(g1)<-NULL
    map.raw<-as.data.frame(g1[,1:3],stringsAsFactors = F)
    map.raw[,2:3]<-sapply(map.raw[,2:3],as.numeric)
    names(map.raw)<-c("marker","chr","pos")
    gen.all<-NULL
    v.map<-NULL
    nchr<-max(as.numeric(mx[,1]))
    for(ii in 1:nchr){
      map.gen<-data.qtl.2$geno[[ii]]$prob
      insert.gen<-map.gen[,,1]-map.gen[,,2]
      gen.all<-cbind(gen.all,insert.gen)
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
    gen<-cbind(v.map[,2],v.map[,3],t(gen.all))
    
  }
  output<-list(yygg=yygg,mx=mx,phe=phe,chrRaw_name=chrRaw_name,chr_name=chr_name,gen=gen,v.map=v.map)
  return(output)
} 

#######################################################

WangS<-function(flag=NULL,CriLOD=NULL,NUM=NULL,pheRaw=NULL,chrRaw_name=NULL,
                yygg=NULL,mx=NULL,phe=NULL,chr_name=NULL,gen=NULL,v.map=NULL){
  
  sLOD<-CriLOD;result=NULL;mxmp=NULL;galaxyy1=NULL;res11=NULL;
  
  phe<-as.matrix(phe[,NUM])
  deletRow<-which(is.na(phe)==TRUE)
  gentwo<-mx[,1:2]
  t_gen<-mx[,3:ncol(mx)]
  
  if(length(deletRow)>0){
    phe<-as.matrix(phe[-deletRow,])
    t_gen1<-t_gen[,-deletRow]
    if(is.null(yygg)==FALSE){
      yygg<-yygg[-deletRow,]
    }
    phe<-phe
    mx<-cbind(gentwo,t_gen1)
    gen<-gen[,-deletRow]
  }else{
    mx<-mx
    phe<-phe
    if(is.null(yygg)==FALSE){
      yygg<-yygg
    }
    gen<-gen
  }
  
  mx<-as.matrix(mx)
  
  mxmp<-mx[,1:2]
  rownames(mxmp)<-NULL
  mxmp<-as.data.frame(mxmp,stringsAsFactors = F)
  mxmp[,1:2]<-sapply(mxmp[,1:2],as.numeric)
  
  
  map<-mx[,1:2]
  geno<-t(mx[,3:(ncol(mx))])
  n_sam<-nrow(geno)
  
  if(is.null(yygg)==FALSE){
    fx<-cbind(matrix(1,n_sam,1),yygg)
  }
  if(is.null(yygg)==TRUE){
    fx<-matrix(1,n_sam,1)
  }
  
  ori<-NULL
  for (j in 1:(nrow(map)))
  {
    ta<-as.matrix(which(gen[,1]==map[j,1]))
    tb<-gen[ta,2]
    cori<-ta[as.matrix(which(tb==map[j,2])),1]
    ori<-rbind(ori,cori)
  }
  corie<-matrix(0,nrow(gen),1)
  corie[ori,1]<-matrix(1,nrow(ori),1)
  gen<-cbind(gen[,1:2],corie,gen[,3:(ncol(gen))])
  wg<-gen
  iw<-as.matrix(which(gen[,3]==1))
  gk<-t(gen[iw,4:(ncol(gen))])
  m<-ncol(gk)
  n<-nrow(gk)
  kk<-matrix(0,n,n)
  for(k in 1:m){
    z<-as.matrix(gk[,k])
    kk<-kk+z%*%t(z)
  }
  cc<-mean(diag(kk))
  kk<-kk/cc
  gen<-cbind(gen[,1:2],gen[,4:(ncol(gen))])
  
  
  if (flag==1)
  {code<-random(fx=fx,gen=gen,phe=phe,kk=kk)
  tempcode<-code
  }
  if (flag==0)
  {code<-fix(x=fx,gen=gen,y=phe,kk=kk)
  tempcode<-code
  tempcode[,2:3] <- gen[,1:2]
  }
  res1<-tempcode
  
  res11<-res1[,c(2,3,10)]
  colnames(res11)<-NULL

  x0<-t(gen[,3:ncol(gen)])
  y<-phe
  bb<-code
  bb<-as.matrix(bb)
  aa<-numeric()
  for (i in 1:nrow(as.matrix(unique(gen[,1])))){
    
    mc<-which(gen[,1]==i)
    mc<-as.matrix(mc)
    for (j in 1:(nrow(mc)-2))
    {if (bb[mc[j+1],10]<bb[mc[j],10] & bb[mc[j+1],10]<bb[mc[j+2],10]){aa<-rbind(aa,mc[j+1])}
    }
    if (bb[mc[1],10]<bb[mc[2],10]){aa<-rbind(aa,mc[1])}
    if (bb[mc[nrow(mc)],10]<bb[mc[nrow(mc)-1],10]){aa<-rbind(aa,mc[nrow(mc)])}
  }
  
  xx<-x0[,aa]
  par<-ebayes_EM(fx,xx,y)
  selectpos<-which(par[,1]<=0.01)
  
  if(length(selectpos)==0){
    warning("No QTL were detected!")
  }else{
    cc<-which(par[,1]<=0.01)
    name<-as.matrix(aa[cc,1])
    xxx<-as.matrix(x0[,name])
    y<-as.matrix(y)
    lod<-likelihood(fx,xxx,y)
    dd<-as.matrix(which(lod[,1]>=sLOD))
    
    if(length(dd)==0){
      warning("No QTL were detected!")
    }else{ 
      
      if(length(dd)==1){
        na<-matrix(name[dd],1,)
        xxxm<-matrix(xxx[,dd],,1)
        wow<-cbind(fx,xxxm)
        bbbb<-solve(t(wow)%*%wow)%*%t(wow)%*%y
        ef<-bbbb[(ncol(fx)+1):nrow(bbbb),1]
        genna1<-matrix(gen[na,1],,1)
        genna2<-matrix(gen[na,2],,1)
        genna3<-as.matrix(ef)
        genna4<-as.matrix(lod[dd,])
        galaxy<-cbind(genna1,genna2,genna3,genna4)
      }else{
        na<-as.matrix(name[dd])
        wow<-cbind(fx,xxx[,dd])
        bbbb<-solve(t(wow)%*%wow)%*%t(wow)%*%y
        ef<-bbbb[(ncol(fx)+1):nrow(bbbb),1]
        galaxy<-cbind(gen[na,1],gen[na,2],ef,lod[dd,])
      }
      galaxy<-galaxy
      c1<-galaxy
      xx1<-numeric()
      for (i in 1:nrow(c1)){
        ng1<-as.matrix(which(gen[,1]==c1[i,1]))
        ng2<-gen[ng1,]
        ng3<-as.matrix(which(ng2[,2]==c1[i,2]))
        xx1<-rbind(xx1,ng2[ng3,])
      }
      xx1<-as.matrix(xx1)
      x1<-matrix(xx1[,3:(ncol(xx1))],,(ncol(xx1)-2))
      x1<-t(x1)
      x1<-as.matrix(x1)
      wow<-cbind(fx,x1)
      bbbb<-solve(t(wow)%*%wow)%*%t(wow)%*%y
      ef<-as.matrix(bbbb[(ncol(fx)+1):nrow(bbbb),1])
      y<-y-x1%*%ef
      
      if (flag==1)
      {code<-random(fx=fx,gen=gen,phe=y,kk=kk)
      }
      if (flag==0)
      {code<-fix(x=fx,gen=gen,y=y,kk=kk)
      }
      
      x0<-t(gen[,3:(ncol(gen))])
      bb<-code
      bb<-as.matrix(bb)
      aa<-numeric()
      for (i in 1:nrow(as.matrix(unique(gen[,1])))){
        
        mc<-which(gen[,1]==i)
        mc<-as.matrix(mc)
        for (j in 1:(nrow(mc)-2))
        {if (bb[mc[j+1],10]<bb[mc[j],10] & bb[mc[j+1],10]<bb[mc[j+2],10]){aa<-rbind(aa,mc[j+1])}
        }
        if (bb[mc[1],10]<bb[mc[2],10]){aa<-rbind(aa,mc[1])}
        if (bb[mc[nrow(mc)],10]<bb[mc[nrow(mc)-1],10]){aa<-rbind(aa,mc[nrow(mc)])}
      }
      
      mi<-code[aa,2:3]
      style<-numeric()
      for (i in 1:nrow(mi))
      {
        for (j in 1:nrow(xx1))
        {if (mi[i,1]==xx1[j,1] & mi[i,2]==xx1[j,2])
        {style<-rbind(style,aa[i])
        }
        }
      }
      aa<-as.matrix(setdiff(aa,style))
      xx<-x0[,aa]
      par<-ebayes_EM(fx,xx,y)
      cc<-as.matrix(which(par[,1]<=0.01))
      selectpos1<-which(par[,1]<=0.01)
      
      if (nrow(cc)>0)
      {
        name<-as.matrix(aa[cc,1])
        xxx<-as.matrix(x0[,name])
        y<-as.matrix(y)
        lod<-likelihood(fx,xxx,y)
        dd<-as.matrix(which(lod[,1]>=sLOD))
        if (nrow(dd)>0)
        {
          na<-as.matrix(name[dd])
          wow<-cbind(fx,xxx[,dd])
          bbbb<-solve(t(wow)%*%wow)%*%t(wow)%*%y
          ef<-bbbb[(ncol(fx)+1):nrow(bbbb),1]
          galaxy2<-cbind(gen[na,1],gen[na,2],ef,lod[dd,])
          galaxyy<-rbind(galaxy,galaxy2)
          woww<-wow
        }else if (nrow(dd)==0)
        {galaxyy<-galaxy
        woww<-wow
        }
      }
      if (nrow(cc)==0)
      {galaxyy<-galaxy
      woww<-wow
      }
      pp<-as.matrix(phe[,1])
      va<-galaxyy[,3]*galaxyy[,3]
      ve<-(1/(n_sam-1))*t(pp-woww%*%bbbb)%*%(pp-woww%*%bbbb)
      vp<-(1/(n_sam-1))*t(pp-mean(pp))%*%(pp-mean(pp))
      vy<-(sum(va)+ve)
      
      if (vy>=vp){
        heredity<-va/as.vector(vy)
        pv<-vy}
      if (vy<vp){
        heredity<-va/as.vector(vp)
        pv<-vp}
      
      va<-matrix(va,,1)
      va<-round(va,4)
      heredity<-100*heredity
      heredity<-matrix(heredity,,1)
      heredity<-round(heredity,4)
      galaxyy[which(abs(galaxyy)>1e-4)]<-round(galaxyy[which(abs(galaxyy)>1e-4)],4)
      galaxyy[which(abs(galaxyy)<1e-4)]<-as.numeric(sprintf("%.4e",galaxyy[which(abs(galaxyy)<1e-4)]))
      
      markerh<-numeric()
      for (i in 1:nrow(galaxyy))
      {
        for (j in 1:nrow(v.map))
        {
          if (galaxyy[i,1]==v.map[j,2] & galaxyy[i,2]==v.map[j,3])
          {
            markerh<-rbind(markerh,v.map[j,6:7])
          }
        }
      }
      
      if((is.null(chrRaw_name)==FALSE)&&(is.null(chr_name)==FALSE)){
        chr_name<-chr_name
        chrRaw_name<-chrRaw_name
        galaxyysec<-galaxyy[,1]
        galaxyylast<-matrix(galaxyy[,2:ncol(galaxyy)],,(ncol(galaxyy)-1))
        chrName<-numeric()
        for( i in 1:length(galaxyysec)){
          chrLoc<-which(chr_name[]==galaxyysec[i])
          chrName0<-matrix(chrRaw_name[chrLoc],,1)
          chrName<-rbind(chrName,chrName0)
        }
        
        galaxyy_A<-cbind(chrName,galaxyylast)
      }else{
        galaxyy_A<-galaxyy
      }
      
      galaxyy_A<-as.matrix(galaxyy_A)
      vee<-matrix("",nrow(galaxyy_A),1)
      vee[1,1]<-round(ve,4)
      vee<-matrix(vee,,1)
      vpp<-matrix("",nrow(galaxyy_A),1)
      vpp[1,1]<-round(pv,4)
      vpp<-matrix(vpp,,1)
      traitid<-matrix(pheRaw[1,NUM+1],nrow(galaxyy_A),1)
      galaxyy<-galaxyy
      galaxyy1<-galaxyy[,c(1,2,4)]
      result<-cbind(traitid,galaxyy_A,markerh,va,heredity,vee,vpp)
      colnames(result)<-c("Trait","Chr","Position (cM)","Additive Effect","LOD","Left_Marker","Right_Marker","Var_Genet_(i)","r2 (%)","Var_Error",
                          "Var_Phen (total)")
      
    }
  } 
  output<-list(result=result,mxmp=mxmp,galaxyy1=galaxyy1,res11=res11,chr_name=chr_name)
  return(output)
}




