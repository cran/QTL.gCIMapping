#' Process raw data
#'
#' @param fileFormat Format of dataset.
#' @param Population Population type.
#' @param method Method "GCIM" or method "GCIM-QEI"
#' @param Model Random or fixed model.
#' @param readraw Raw data.
#' @param MultiEnv Whether to perform multi-environment analysis
#'
#' @return  a list
#'
#' @examples
#' data(F2data)
#' readraw<-Readdata(file=F2data,fileFormat="GCIM",
#' method="GCIM-QEI",filecov=NULL,
#' MCIMmap=NULL,MultiEnv=TRUE)
#' doda<-Dodata(fileFormat="GCIM",Population="F2",
#' method="GCIM-QEI",Model="Random",
#' readraw,MultiEnv=TRUE)
Dodata<-function(fileFormat=NULL,Population=NULL,method=NULL,Model=NULL,readraw=NULL,MultiEnv=FALSE){
  if(method=="GCIM-QEI"){
    pheRaw=NULL;genRaw=NULL;mapRaw1=NULL
    cov_en=NULL;yygg1=NULL
    geoo<-readraw$geoo;pho<-readraw$pho;poss<-readraw$poss;parm<-readraw$parm
    MCIM_data<-readraw$MCIM_data;map_data<-readraw$map_data
    y_jun3<-readraw$y_jun3
    genRaw1<-readraw$genRaw1;pheRaw<-readraw$pheRaw;mapRaw11<-readraw$mapRaw11
    cov_en<-readraw$cov_en
    if(MultiEnv==TRUE){
      if(fileFormat=="ICIM"){
        ################## fileFormat=="ICIM" ################
        if(parm[3,1]==1){# 1 for intervals
          pos.be<-numeric()
          ChrNum<-unique(poss[,2])
          for(i in 1:length(ChrNum)){
            poss1<-poss[which(poss[,2]==i),]
            positi<-cumsum(poss1[,3,drop=F])
            poss2<-cbind(poss1[,1:2,drop=F],positi)
            pos.be<-rbind(pos.be,poss2)
          }
        }
        if(parm[3,1]==2){# 2 for positions
          pos.be<-poss
        }

        if(parm[4,1]==2){# 2 for Morgan
          pos<-cbind(pos.be[,1:2],matrix(100*pos.be[,3],,1))
        }
        if(parm[4,1]==1){# 1 for centiMorgan
          pos<-pos.be
        }# dim(pos)

        # mapRaw1<-as.matrix(pos)
        mapRaw1<-pos
        colnames(mapRaw1)<-c("marker","chr","pos")# dim(mapRaw1)
        mapRaw1<-as.data.frame(mapRaw1)

        genRaw<-geoo[,-1] # dim(genRaw)  dim(geoo)
        if(Population=="BC1"){
          genRaw<-gsub("-1","99",genRaw)
          genRaw<-gsub("1","-1",genRaw)
          genRaw<-gsub("2","1",genRaw)
        }else if(Population=="BC2"){
          genRaw<-gsub("-1","99",genRaw)
          #gen_2<-gsub("1","1",gen_1)
          genRaw<-gsub("0","-1",genRaw)
        }else if(Population=="DH"){
          genRaw<-gsub("-1","99",genRaw)
          genRaw<-gsub("0","-1",genRaw)
          genRaw<-gsub("2","1",genRaw)
        }else if(Population=="RIL"){
          genRaw<-gsub("-1","99",genRaw)
          genRaw<-gsub("0","-1",genRaw)
          genRaw<-gsub("2","1",genRaw)
        }else if(Population=="F2"){

          genRaw<-gsub("2","A",genRaw)
          genRaw<-gsub("1","H",genRaw)
          genRaw<-gsub("0","B",genRaw)
          genRaw<-gsub("-1","-",genRaw)

          genRaw<-gsub("12","D",genRaw)
          genRaw<-gsub("10","C",genRaw)

          genRaw[genRaw=="XX"]<-"-"
          genRaw[genRaw=="X"]<-"-"
          genRaw[genRaw=="**"]<-"-"
          genRaw[genRaw=="*"]<-"-"# Missing values of marker type are coded as -1, X, XX, *, or **
        }
        genRaw<-as.matrix(genRaw)

        pheRaw<-as.matrix(t(pho)[-1,])
        OrigName<-pho[,1]
        pheRaw<-gsub(-100,NA,pheRaw)# dim(pheRaw)
        pheRaw[pheRaw=="*"]<-"-"
        pheRaw[pheRaw=="."]<-"-"
        EnvNum<-as.numeric(parm[7,1])
        TraitNum<-1
        # traitName<-NULL
        if(dim(pheRaw)[2]!=EnvNum){
          warning("please check the phenotype in file!")
        }

        pheRaw<-as.matrix(pheRaw)
        rownames(pheRaw)<-NULL;colnames(pheRaw)<-OrigName

        if(is.null(cov_en)==FALSE){
          covnum<-t(cov_en[-1,2:ncol(cov_en)])
          yygg1<-numeric()
          for(i in 1:nrow(covnum)){
            otrait_ind<-unique(covnum[i,])
            cov_col<-length(otrait_ind)-1
            col_each<-numeric()
            for(j in 1:length(covnum[i,])){
              if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                cov_0<-matrix(-1,1,cov_col)
              }else{
                cov_0<-matrix(0,1,cov_col)
                covnum_loc<-which(otrait_ind[]==covnum[i,j])
                cov_0[1,covnum_loc]<-1
              }
              col_each<-rbind(col_each,cov_0)

            }
            yygg1<-cbind(yygg1,col_each)# dim(yygg1)
          }
          cov_en<-covnum
        }else{
          yygg1<-NULL
        }



      }else if(fileFormat=="MCIM"){
        ################## fileFormat=="MCIM" ################
        start_dex<-grep("*MapBegin*",map_data,fixed = TRUE)
        stop_dex<-grep("*MapEnd*",map_data,fixed = TRUE)
        chr_dex<-grep("_Chromosomes",map_data,fixed = TRUE)
        chr_num<-as.numeric(map_data[chr_dex+1])
        MarkerNum_dex<-grep("_MarkerNumbers",map_data,fixed = TRUE)
        Marker_num<-as.numeric(map_data[(MarkerNum_dex+1):(MarkerNum_dex+chr_num)])
        maxMarker_num<-max(Marker_num)
        posMatrix<-matrix(NA,nrow=maxMarker_num,ncol=chr_num)
        chrStart_dex<-grep(paste("Ch",chr_num,sep=""),map_data,fixed = TRUE)
        if(length(chrStart_dex)<=0){
          chrStart_dex<-grep(paste("ch",chr_num,sep=""),map_data,fixed = TRUE)
        }
        ChrPos_data<-as.numeric(map_data[(chrStart_dex+1):(stop_dex-1)])
        last_dex<-0
        for(i in 1:maxMarker_num){
          ni<-which(Marker_num>=i)
          posMatrix[i,ni]<-ChrPos_data[last_dex+seq(1,length(ni),1)+1]
          if(i!=ChrPos_data[last_dex+1]){
            warning("please check the map!")
          }
          last_dex<-last_dex+length(ni)+1
        }
        sorted_map0<-cbind(sort(rep(seq(1,chr_num,1),maxMarker_num)),matrix(posMatrix,ncol=1,byrow=F))
        sorted_map<-as.matrix(na.omit(sorted_map0))[,1:2]
        colnames(sorted_map)<-NULL;rownames(sorted_map)<-NULL

        pos.be<-numeric()
        ChrNum<-unique(sorted_map[,1])
        for(i in 1:length(ChrNum)){
          poss1<-sorted_map[which(sorted_map[,1]==i),]
          positi<-cumsum(poss1[,2,drop=F])
          poss2<-cbind(poss1[,1,drop=F],positi)
          pos.be<-rbind(pos.be,poss2)
        }

        indivNum<-as.numeric(MCIM_data[grep("_Genotypes",MCIM_data,fixed = TRUE)+1])
        GenStart_dex<-grep("*MarkerBegin*",MCIM_data,fixed = TRUE)
        GenStop_dex<-grep("*MarkerEnd*",MCIM_data,fixed = TRUE)
        genRaw<-MCIM_data[(GenStart_dex+1):(GenStop_dex-1)]
        semic_dex<-grep(";",genRaw,fixed = TRUE)
        StrNum<-semic_dex[2]-semic_dex[1]
        MarNum<-as.numeric(MCIM_data[grep("_TotalMarker",MCIM_data,fixed = TRUE)+1])
        if((StrNum-2)==MarNum){
          genRaw<-matrix(genRaw,nrow=indivNum+1,ncol=dim(pos.be)[1]+2,byrow=T)
          MarkerName<-genRaw[1,2:(dim(pos.be)[1]+1)]
          mapRaw1<-data.frame(marker=MarkerName,chr=pos.be[,1],pos=pos.be[,2])# dim(mapRaw1)
          genRaw<-genRaw[2:(indivNum+1),2:(dim(mapRaw1)[1]+1)]# dim(genRaw)
          colnames(genRaw)<-NULL;rownames(genRaw)<-NULL
        }else{
          genRaw<-matrix(genRaw,nrow=dim(pos.be)[1]+1,ncol=indivNum+2,byrow=T)
          MarkerName<-genRaw[2:(dim(pos.be)[1]+1),1]
          mapRaw1<-data.frame(marker=MarkerName,chr=pos.be[,1],pos=pos.be[,2])# dim(mapRaw1)
          genRaw<-genRaw[2:(dim(mapRaw1)[1]+1),2:(indivNum+1)]# dim(genRaw)
          colnames(genRaw)<-NULL;rownames(genRaw)<-NULL
        }
        genRaw[genRaw=="."]<-"-"# Missing values of marker type are coded as .
        genRaw<-t(genRaw)

        PheStart_dex<-grep("*TraitBegin*",MCIM_data,fixed = TRUE)
        PheStop_dex<-grep("*TraitEnd*",MCIM_data,fixed = TRUE)
        PheData<-MCIM_data[(PheStart_dex+1):(PheStop_dex-1)]
        semicolon_dex<-grep(";",PheData,fixed = TRUE)
        StringNum<-semicolon_dex[2]-semicolon_dex[1]
        TraitNum<-as.numeric(MCIM_data[grep("_TraitNumber",MCIM_data,fixed = TRUE)+1])
        PheData<-matrix(PheData,ncol=StringNum,byrow=T)
        pheRaw<-PheData[,(StringNum-TraitNum):(StringNum-1),drop=F]
        traitName<-pheRaw[1,]
        pheRaw<-pheRaw[-1,,drop=F]
        # pheRaw<-as.data.frame(pheRaw)
        # pheRaw<-apply(pheRaw,2,as.numeric)
        EnvRaw<-PheData[-1,1,drop=F]
        EnvName<-unique(EnvRaw)
        EnvNum<-length(EnvName)
        bb<-rep(EnvNum,dim(pheRaw)[2])

        aa<-NULL
        for(i in 1:length(traitName)){
          aa<-c(aa,rep(traitName[i],EnvNum))
        }
        traitName<-aa

        if(StringNum>(3+TraitNum)){
          if(MCIM_data[grep("_Replications",MCIM_data,fixed = TRUE)+1]=="yes"){
            newPhe<-NULL
            newGeo<-NULL
            newEN<-NULL
            RepRaw<-PheData[-1,2,drop=F]
            for(i in 1:EnvNum){
              now_id<-which(EnvRaw==EnvName[i])
              now_Rep<-RepRaw[now_id]
              ns<-unique(now_Rep)
              now_Geo<-as.numeric(PheData[-1,,drop=F][now_id,(StringNum-TraitNum-1)])
              now_Phe<-pheRaw[now_id,,drop=F]
              if(length(ns)>1){
                nr<-sort(unique(now_Geo))
                newEN<-c(newEN,rep(i,length(nr)))
                for(j in 1:length(nr)){
                  now_id2<-which(now_Geo==nr[j])
                  meanphe<-now_Phe[now_id2,,drop=F]
                  if(dim(meanphe)[1]>1){
                    meanphe<-apply(meanphe,2,as.numeric)
                    meanphe<-apply(meanphe,2,function(x)mean(x,na.rm=T))
                    meanphe<-matrix(meanphe,nrow=1)
                    newPhe<-rbind(newPhe, meanphe)
                    newGeo<-rbind(newGeo,nr[j])
                  }else{
                    # meanphe<-apply(meanphe,2,as.numeric)
                    # meanphe<-apply(meanphe,2,function(x)mean(x,na.rm=T))
                    # meanphe<-matrix(meanphe,nrow=1)
                    newPhe<-rbind(newPhe, meanphe)
                    newGeo<-rbind(newGeo,nr[j])
                  }

                }
              }else{
                newEN<-c(newEN,rep(i,length(now_Geo)))
                newPhe<-rbind(newPhe,now_Phe)
                newGeo<-rbind(newGeo,matrix(now_Geo,ncol=1))
              }
            }
          }else{
            warning("please check the phenotype or the parameter (_Replications) in file!")
          }

        }else{
          newPhe<-pheRaw
          newGeo<-as.numeric(PheData[-1,2])
          newEN<-NULL
          for(i in 1:EnvNum){

            now_id<-which(EnvRaw==EnvName[i])
            newEN<-c(newEN,rep(i,length(now_id)))
          }
        }
        newID<-(newEN-1)*indivNum+newGeo

        pheRaw_All<-matrix(NA,nrow=indivNum*EnvNum,ncol=dim(pheRaw)[2])
        # pheRaw_All[newID,,drop=F]<-newPhe
        pheRaw_All[newID,]<-newPhe
        # pheRaw<-as.matrix(pheRaw_All)
        pheRaw<-matrix(pheRaw_All,nrow=indivNum,ncol=dim(pheRaw)[2]*EnvNum,byrow=FALSE)
        pheRaw[which(is.na(pheRaw))]<-"-"
        pheRaw[pheRaw=="."]<-"-"

        colnames(pheRaw)<-traitName;rownames(pheRaw)<-NULL# dim(pheRaw)

        EnvNum<-bb

        if(is.null(cov_en)==FALSE){
          covnum<-t(cov_en[-1,2:ncol(cov_en)])
          yygg1<-numeric()
          for(i in 1:nrow(covnum)){
            otrait_ind<-unique(covnum[i,])
            cov_col<-length(otrait_ind)-1
            col_each<-numeric()
            for(j in 1:length(covnum[i,])){
              if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                cov_0<-matrix(-1,1,cov_col)
              }else{
                cov_0<-matrix(0,1,cov_col)
                covnum_loc<-which(otrait_ind[]==covnum[i,j])
                cov_0[1,covnum_loc]<-1
              }
              col_each<-rbind(col_each,cov_0)

            }
            yygg1<-cbind(yygg1,col_each)# dim(yygg1)
          }
          cov_en<-covnum
        }else{
          yygg1<-NULL
        }

      }else if(fileFormat=="Cart"){
        ################## fileFormat=="Cart" ################
        start_dex<-grep("-start",y_jun3,fixed = TRUE)
        stop_dex<-grep("-stop",y_jun3,fixed = TRUE)
        chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)

        chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome")
        chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)

        chrdata_dexlen<-as.numeric(y_jun3[grep("-chromosomes",y_jun3,fixed = TRUE)+1])

        chr_num<-numeric()
        chr_numfirst<-numeric()
        markername0<-numeric()
        chr_pos<-numeric()

        for(i in 1:chrdata_dexlen){
          chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
          marker_name<-numeric()
          marker_pos<-numeric()

          for(j in 0:(trunc(length(chr_each)/2)-1) ){
            marker_name<-cbind(marker_name,chr_each[2*j+1])
            marker_pos<-cbind(marker_pos,suppressWarnings(as.numeric(chr_each[2*(j+1)])))

            if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){

              marker_posm<-100*((-0.5)*log(1-2*marker_pos))
              markerlen<-length(marker_posm)
              marker_pos1<-marker_posm[1:(markerlen-1)]
              marker_pos2<-c(0,marker_pos1)
              marker_possum<-cumsum(marker_pos2)
            }
            if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
              marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))
              markerlen<-length(marker_posm)
              marker_pos1<-marker_posm[1:(markerlen-1)]
              marker_pos2<-c(0,marker_pos1)
              marker_possum<-cumsum(marker_pos2)
            }
            if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
              # y_jun3[9]=="M"  Morgan. The distance over which, on average, one crossover occurs per meiosis.
              marker_posm<-100*marker_pos
              markerlen<-length(marker_posm)
              marker_pos1<-marker_posm[1:(markerlen-1)]
              marker_pos2<-c(0,marker_pos1)
              marker_possum<-cumsum(marker_pos2)
            }
            if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
              marker_posm<-marker_pos
              markerlen<-length(marker_posm)
              marker_pos1<-marker_posm[1:(markerlen-1)]
              marker_pos2<-c(0,marker_pos1)
              marker_possum<-cumsum(marker_pos2)
            }
            if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
              marker_possum<-100*marker_pos
            }
            if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
              marker_possum<-marker_pos
            }
          }
          markername0<-cbind(markername0,marker_name)

          chr_pos<-rbind(chr_pos,matrix(marker_possum,,1))
          chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
          chr_data<-na.omit(chr_a)
          chr_num<-rbind(chr_num,length(chr_data))

          chr_numfirst<-rbind(chr_numfirst,matrix(rep(i,chr_num[i]),,1))
        }

        chr_leng<-length(chr_pos)

        mapRaw1<-data.frame(marker=matrix(markername0,ncol=1),chr=chr_numfirst,pos=chr_pos)# dim(mapRaw1)

        marker_dex<-grep("markers",y_jun3,fixed = TRUE)
        marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]

        indi_num<-as.numeric(y_jun3[grep("-SampleSize",y_jun3,fixed = TRUE)+1])

        if_indi<-y_jun3[marker_dex[1]-1]
        if(if_indi=="individuals"){
          genRaw<-matrix(marker_snp[-c(1:3,length(marker_snp))],nrow=indi_num,byrow = T)
          genRaw<-genRaw[,-1,drop=F]
          genRaw<-t(genRaw)# dim(genRaw)
        }else{
          genRaw<-matrix(marker_snp[-c(1:2,length(marker_snp))],ncol=indi_num+1,byrow = T)
          genRaw<-genRaw[,-1,drop=F]
        }# dim(genRaw)
        rownames(genRaw)<-NULL;colnames(genRaw)<-NULL# class(genRaw)
        genRaw<-apply(genRaw,2,as.character)# class(genRaw)

        trait_total<-y_jun3[start_dex[3]:stop_dex[3]]
        trait_total[which(trait_total==".")]<-"-"
        trait_total[which(trait_total=="*")]<-"-"

        trait_dex<-grep("traits",trait_total)
        iftrait_indi<-trait_total[trait_dex[1]-1]
        if(iftrait_indi=="individuals"){

          name_dex<-grep("named",trait_total)
          traitdata<-trait_total[(name_dex+1):(length(trait_total)-1)]
          pheRaw<-matrix(traitdata,nrow=indi_num,byrow = T)
          pheRaw<-pheRaw[,-1,drop=F]# dim(pheRaw)
          trait_num<-dim(pheRaw)[2]

          traitName<-trait_total[(trait_dex+2):(name_dex-1)]

          TraitNum<-0
          for(i in 1:length(traitName)){
            traitStr<-strsplit(traitName[i],"E")[[1]]
            # EnvNum<-max(as.numeric(traitStr[2]),EnvNum)
            TraitNum<-max(as.numeric(strsplit(traitStr[1],"t")[[1]][2]),TraitNum)
          }
          EnvNum<-rep(0,TraitNum)
          traitStrAll<-strsplit(traitName,"E")
          HalfName<-paste("t",seq(1,TraitNum),sep="")
          for(i in 1:length(traitName)){
            EnvID<-which(HalfName%in%traitStrAll[[i]][1])
            EnvNum[EnvID]<-EnvNum[EnvID]+1
          }

          newtraitName<-NULL
          for(i in 1:TraitNum){
            for(j in 1:EnvNum[i]){
              newtraitName<-c(newtraitName,paste("t",i,"E",j,sep=""))
            }
          }
          NameID<-match(newtraitName,traitName)
          pheRaw<-pheRaw[,NameID]# dim(pheRaw)
          traitName<-traitName[NameID]
          # rownames(pheRaw)<-NULL;colnames(pheRaw)<-NULL

        }else{
          traitdata<-trait_total[(trait_dex[1]+1):(length(trait_total)-1)]
          pheRaw<-matrix(traitdata,ncol=indi_num+1,byrow = T)
          traitName<-c(pheRaw[,1])
          pheRaw<-t(pheRaw[,-1,drop=F])# dim(pheRaw)
          trait_num<-dim(pheRaw)[2]

          TraitNum<-0
          for(i in 1:length(traitName)){
            traitStr<-strsplit(traitName[i],"E")[[1]]
            # EnvNum<-max(as.numeric(traitStr[2]),EnvNum)
            TraitNum<-max(as.numeric(strsplit(traitStr[1],"t")[[1]][2]),TraitNum)
          }
          EnvNum<-rep(0,TraitNum)
          traitStrAll<-strsplit(traitName,"E")
          HalfName<-paste("t",seq(1,TraitNum),sep="")
          for(i in 1:length(traitName)){
            EnvID<-which(HalfName%in%traitStrAll[[i]][1])
            EnvNum[EnvID]<-EnvNum[EnvID]+1
          }

          newtraitName<-NULL
          for(i in 1:TraitNum){
            for(j in 1:EnvNum[i]){
              newtraitName<-c(newtraitName,paste("t",i,"E",j,sep=""))
            }
          }
          NameID<-match(newtraitName,traitName)
          pheRaw<-pheRaw[,NameID]# dim(pheRaw)
          traitName<-traitName[NameID]
        }
        rownames(pheRaw)<-NULL;colnames(pheRaw)<-traitName

        if(length(start_dex)==3){yygg1<-NULL}
        if(length(start_dex)==4){
          if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
            cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
            cov_dex<-grep("otraits",cov_total)
            cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
            bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
            otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
            if(y_jun3[start_dex[4]+1]=="otraits"){

              covnumonly<-numeric()
              for( i in 0:(otrait_indi-1)){
                cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
                covnumonly<-rbind(covnumonly,cov_each)
              }
              covnum<-covnumonly[,-1]
              yygg1<-numeric()
              for(i in 1:nrow(covnum)){

                otrait_ind<-unique(covnum[i,])
                cov_col<-length(otrait_ind)-1

                col_each<-numeric()
                for(j in 1:length(covnum[i,])){

                  if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                    cov_0<-matrix(-1,1,cov_col)

                  }else{
                    cov_0<-matrix(0,1,cov_col)
                    covnum_loc<-which(otrait_ind[]==covnum[i,j])
                    cov_0[1,covnum_loc]<-1
                  }
                  col_each<-rbind(col_each,cov_0)

                }
                yygg1<-cbind(yygg1,col_each)
              }

            }
            if(y_jun3[start_dex[4]+1]=="individuals"){
              named_dex<-grep("named",cov_only)
              covdata<-cov_only[(named_dex+1):length(cov_only)]
              covnum<-numeric()
              # otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)

              for(n in 0:(indi_num-1)){
                cov_each<-matrix(covdata[(2+(otrait_indi+1)*n):(1+(otrait_indi+1)*n+otrait_indi)],otrait_indi,1)
                covnum<-cbind(covnum,cov_each)
              }
              yygg1<-numeric()
              for(i in 1:nrow(covnum)){

                otrait_ind<-unique(covnum[i,])
                cov_col<-length(otrait_ind)-1

                col_each<-numeric()
                for(j in 1:length(covnum[i,])){

                  if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                    cov_0<-matrix(-1,1,cov_col)

                  }else{
                    cov_0<-matrix(0,1,cov_col)
                    covnum_loc<-which(otrait_ind[]==covnum[i,j])
                    cov_0[1,covnum_loc]<-1
                  }
                  col_each<-rbind(col_each,cov_0)

                }
                yygg1<-cbind(yygg1,col_each)

              }
            }

          }
        }

        if(Population=="BC1"){
          genRaw<-gsub("1","-1",genRaw)
          genRaw<-gsub("2","1",genRaw)
          genRaw<-gsub("*","99",genRaw,fixed = TRUE)# dim(genRaw)
        }else if(Population=="BC2"){
          genRaw<-gsub("0","-1",genRaw)
          #marker_snp2<-gsub("1","1",marker_snp1)
          genRaw<-gsub("*","99",genRaw,fixed = TRUE)
        }else if(Population=="DH"){
          genRaw<-gsub("0","-1",genRaw)
          genRaw<-gsub("2","1",genRaw)
          genRaw<-gsub("*","99",genRaw,fixed = TRUE)
        }else if(Population=="RIL"){
          genRaw<-gsub("1","99",genRaw)
          genRaw<-gsub("0","-1",genRaw)
          genRaw<-gsub("2","1",genRaw)
          genRaw<-gsub("*","99",genRaw,fixed = TRUE)
        }else if(Population=="F2"){
          genRaw<-gsub("12","D",genRaw)
          genRaw<-gsub("10","C",genRaw)
          genRaw<-gsub("0","B",genRaw)
          genRaw<-gsub("2","A",genRaw)
          genRaw<-gsub("-1","-",genRaw)
          genRaw<-gsub("1","H",genRaw,fixed = TRUE)
        }

        if(is.null(yygg1)==FALSE){
          cov_en<-covnum
        }else{
          cov_en<-NULL
        }

      }else if(fileFormat=="GCIM"){
        ################## fileFormat=="GCIM" ################
        mapRaw1<-data.frame(marker=mapRaw11[,1],chr=mapRaw11[,2],pos=mapRaw11[,3])# dim(mapRaw1)

        pheRaw<-as.matrix(pheRaw)
        TraitName<-unique(pheRaw[2,-1])
        TraitNum<-length(TraitName)# Trait_Name<-unique(pheRaw[2,-1])
        EnvNum<-NULL
        for(i in 1:length(TraitName)){
          EnvNum<-c(EnvNum,length(which(pheRaw[2,-1]==TraitName[i])))
        }
        OrigName<-pheRaw[3,-1]

        TraitColname<-NULL
        pheorder<-NULL
        for(i in 1:TraitNum){
          for(j in 1:EnvNum[i]){
            TraitColname<-c(TraitColname,paste("t",i,"E",j,sep=""))
            pheorder<-c(pheorder,which(((pheRaw[2,-1]==paste("trait",i,sep=""))|((pheRaw[2,-1]==paste("Trait",i,sep=""))))
                                       &((pheRaw[1,-1]==paste("Env",j,sep=""))|(pheRaw[1,-1]==paste("env",j,sep=""))))
            )
          }
        }

        pheRaw<-pheRaw[4:dim(pheRaw)[1],-1]# dim(pheRaw)
        pheRaw[which(is.na(pheRaw))]<-"-"
        pheRaw<-pheRaw[,pheorder,drop=F]
        # pheRaw<-as.data.frame(pheRaw)
        pheRaw<-as.matrix(pheRaw)
        colnames(pheRaw)<-OrigName[pheorder];rownames(pheRaw)<-NULL# dim(pheRaw)

        genRaw<-as.matrix(genRaw1)

        # if(Population=="BC1"){
        #   genRaw<-gsub("-","99",genRaw)
        #   genRaw<-gsub("H","-1",genRaw)
        #   genRaw<-gsub("A","1",genRaw)
        # }else if(Population=="BC2"){
        #   genRaw<-gsub("-","99",genRaw)
        #   genRaw<-gsub("B","-1",genRaw)
        #   genRaw<-gsub("H","1",genRaw)
        # }else if(Population=="DH"){
        #   genRaw<-gsub("-","99",genRaw)
        #   genRaw<-gsub("B","-1",genRaw)
        #   genRaw<-gsub("A","1",genRaw)
        # }else if(Population=="RIL"){
        #   genRaw<-gsub("-","99",genRaw)
        #   genRaw<-gsub("B","-1",genRaw)
        #   genRaw<-gsub("A","1",genRaw)
        # }else if(Population=="F2"){
        #
        # }

        if(is.null(cov_en)==FALSE){
          covnum<-t(cov_en)
          yygg1<-numeric()
          for(i in 1:nrow(covnum)){
            otrait_ind<-unique(covnum[i,])
            cov_col<-length(otrait_ind)-1
            col_each<-numeric()
            for(j in 1:length(covnum[i,])){
              if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                cov_0<-matrix(-1,1,cov_col)
              }else{
                cov_0<-matrix(0,1,cov_col)
                covnum_loc<-which(otrait_ind[]==covnum[i,j])
                cov_0[1,covnum_loc]<-1
              }
              col_each<-rbind(col_each,cov_0)

            }
            yygg1<-cbind(yygg1,col_each)
          }
        }else{
          yygg1<-NULL
        }
      }

      result<-list(pheRaw=pheRaw,genRaw=genRaw,mapRaw1=mapRaw1,yygg1=yygg1,EnvNum=EnvNum)# cov_en=cov_en,

    }else{
      if(fileFormat=="ICIM"){
        ################## fileFormat=="ICIM" ################
        if(parm[3,1]==1){# 1 for intervals
          pos.be<-numeric()
          ChrNum<-unique(poss[,2])
          for(i in 1:length(ChrNum)){
            poss1<-poss[which(poss[,2]==i),]
            positi<-cumsum(poss1[,3,drop=F])
            poss2<-cbind(poss1[,1:2,drop=F],positi)
            pos.be<-rbind(pos.be,poss2)
          }
        }
        if(parm[3,1]==2){
          pos.be<-poss
        }
        if(parm[4,1]==2){# 2 for Morgan
          pos<-cbind(pos.be[,1:2],matrix(100*pos.be[,3],,1))
        }
        if(parm[4,1]==1){# 1 for centiMorgan
          pos<-pos.be
        }# dim(pos)

        # mapRaw1<-data.frame(marker=pos[,1],chr=pos[,2],pos=pos[,3])# dim(mapRaw1)
        mapRaw1<-pos
        colnames(mapRaw1)<-c("marker","chr","pos")# dim(mapRaw1)
        mapRaw1<-as.data.frame(mapRaw1)

        genRaw<-geoo[,-1] # dim(genRaw)  dim(geoo)
        if(Population=="BC1"){
          genRaw<-gsub("-1","99",genRaw)
          genRaw<-gsub("1","-1",genRaw)
          genRaw<-gsub("2","1",genRaw)
        }else if(Population=="BC2"){
          genRaw<-gsub("-1","99",genRaw)
          #gen_2<-gsub("1","1",gen_1)
          genRaw<-gsub("0","-1",genRaw)
        }else if(Population=="DH"){
          genRaw<-gsub("-1","99",genRaw)
          genRaw<-gsub("0","-1",genRaw)
          genRaw<-gsub("2","1",genRaw)
        }else if(Population=="RIL"){
          genRaw<-gsub("-1","99",genRaw)
          genRaw<-gsub("0","-1",genRaw)
          genRaw<-gsub("2","1",genRaw)
        }else if(Population=="F2"){
          genRaw<-gsub("2","A",genRaw)
          genRaw<-gsub("1","H",genRaw)
          genRaw<-gsub("0","B",genRaw)
          genRaw<-gsub("-1","-",genRaw)

          genRaw<-gsub("12","D",genRaw)
          genRaw<-gsub("10","C",genRaw)

          genRaw[genRaw=="XX"]<-"-"
          genRaw[genRaw=="X"]<-"-"
          genRaw[genRaw=="**"]<-"-"
          genRaw[genRaw=="*"]<-"-"# Missing values of marker type are coded as -1, X, XX, *, or **
        }
        genRaw<-as.matrix(genRaw)

        pheRaw<-as.matrix(t(pho)[-1,])
        OrigName<-pho[,1]
        pheRaw<-gsub(-100,NA,pheRaw)# dim(pheRaw)
        pheRaw[pheRaw=="*"]<-"-"
        pheRaw[pheRaw=="."]<-"-"
        EnvNum<-as.numeric(parm[7,1])
        TraitNum<-1
        # traitName<-NULL
        if(dim(pheRaw)[2]!=EnvNum){
          warning("please check the phenotype in file!")
        }
        # for(i in 1:EnvNum){
        #   traitName<-c(traitName,paste("t1E",i,sep=""))
        #
        # }
        pheRaw<-as.matrix(pheRaw)
        rownames(pheRaw)<-NULL;colnames(pheRaw)<-OrigName

        if(is.null(cov_en)==FALSE){
          covnum<-t(cov_en[-1,2:ncol(cov_en)])
          yygg1<-numeric()
          for(i in 1:nrow(covnum)){
            otrait_ind<-unique(covnum[i,])
            cov_col<-length(otrait_ind)-1
            col_each<-numeric()
            for(j in 1:length(covnum[i,])){
              if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                cov_0<-matrix(-1,1,cov_col)
              }else{
                cov_0<-matrix(0,1,cov_col)
                covnum_loc<-which(otrait_ind[]==covnum[i,j])
                cov_0[1,covnum_loc]<-1
              }
              col_each<-rbind(col_each,cov_0)
            }
            yygg1<-cbind(yygg1,col_each)# dim(yygg1)
          }
          cov_en<-covnum
        }else{
          yygg1<-NULL
        }

      }else if(fileFormat=="MCIM"){
        ################## fileFormat=="MCIM" ################
        start_dex<-grep("*MapBegin*",map_data,fixed = TRUE)
        stop_dex<-grep("*MapEnd*",map_data,fixed = TRUE)
        chr_dex<-grep("_Chromosomes",map_data,fixed = TRUE)
        chr_num<-as.numeric(map_data[chr_dex+1])
        MarkerNum_dex<-grep("_MarkerNumbers",map_data,fixed = TRUE)
        Marker_num<-as.numeric(map_data[(MarkerNum_dex+1):(MarkerNum_dex+chr_num)])
        maxMarker_num<-max(Marker_num)
        posMatrix<-matrix(NA,nrow=maxMarker_num,ncol=chr_num)
        chrStart_dex<-grep(paste("Ch",chr_num,sep=""),map_data,fixed = TRUE)
        if(length(chrStart_dex)<=0){
          chrStart_dex<-grep(paste("ch",chr_num,sep=""),map_data,fixed = TRUE)
        }
        ChrPos_data<-as.numeric(map_data[(chrStart_dex+1):(stop_dex-1)])
        last_dex<-0
        for(i in 1:maxMarker_num){
          ni<-which(Marker_num>=i)
          posMatrix[i,ni]<-ChrPos_data[last_dex+seq(1,length(ni),1)+1]
          if(i!=ChrPos_data[last_dex+1]){
            warning("please check the map!")
          }
          last_dex<-last_dex+length(ni)+1
        }
        sorted_map0<-cbind(sort(rep(seq(1,chr_num,1),maxMarker_num)),matrix(posMatrix,ncol=1,byrow=F))
        sorted_map<-as.matrix(na.omit(sorted_map0))[,1:2]
        colnames(sorted_map)<-NULL;rownames(sorted_map)<-NULL

        # NaNum<-which(is.na(sorted_map0[,2]))
        # if(length(NaNum)>0){
        #   sorted_map<-sorted_map0[-NaNum,]
        # }else{
        #   sorted_map<-sorted_map0
        # }
        pos.be<-numeric()
        ChrNum<-unique(sorted_map[,1])
        for(i in 1:length(ChrNum)){
          poss1<-sorted_map[which(sorted_map[,1]==i),]
          positi<-cumsum(poss1[,2,drop=F])#
          poss2<-cbind(poss1[,1,drop=F],positi)
          pos.be<-rbind(pos.be,poss2)
        }
        # mapRaw1<-data.frame(marker=paste("Bin",seq(1,dim(pos.be)[1]),sep=""),chr=pos.be[,1],pos=pos.be[,2])# dim(mapRaw1)

        indivNum<-as.numeric(MCIM_data[grep("_Genotypes",MCIM_data,fixed = TRUE)+1])
        GenStart_dex<-grep("*MarkerBegin*",MCIM_data,fixed = TRUE)
        GenStop_dex<-grep("*MarkerEnd*",MCIM_data,fixed = TRUE)
        genRaw<-MCIM_data[(GenStart_dex+1):(GenStop_dex-1)]
        semic_dex<-grep(";",genRaw,fixed = TRUE)
        StrNum<-semic_dex[2]-semic_dex[1]
        MarNum<-as.numeric(MCIM_data[grep("_TotalMarker",MCIM_data,fixed = TRUE)+1])
        if((StrNum-2)==MarNum){
          genRaw<-matrix(genRaw,nrow=indivNum+1,ncol=dim(pos.be)[1]+2,byrow=T)
          MarkerName<-genRaw[1,2:(dim(pos.be)[1]+1)]
          mapRaw1<-data.frame(marker=MarkerName,chr=pos.be[,1],pos=pos.be[,2])# dim(mapRaw1)
          genRaw<-genRaw[2:(indivNum+1),2:(dim(mapRaw1)[1]+1)]# dim(genRaw)
          colnames(genRaw)<-NULL;rownames(genRaw)<-NULL
        }else{
          genRaw<-matrix(genRaw,nrow=dim(pos.be)[1]+1,ncol=indivNum+2,byrow=T)
          MarkerName<-genRaw[2:(dim(pos.be)[1]+1),1]
          mapRaw1<-data.frame(marker=MarkerName,chr=pos.be[,1],pos=pos.be[,2])# dim(mapRaw1)
          genRaw<-genRaw[2:(dim(mapRaw1)[1]+1),2:(indivNum+1)]# dim(genRaw)
          colnames(genRaw)<-NULL;rownames(genRaw)<-NULL
        }
        genRaw[genRaw=="."]<-"-"# Missing values of marker type are coded as .
        # if(Population=="F2"){flagRIL<-0}
        genRaw<-t(genRaw)

        PheStart_dex<-grep("*TraitBegin*",MCIM_data,fixed = TRUE)
        PheStop_dex<-grep("*TraitEnd*",MCIM_data,fixed = TRUE)
        PheData<-MCIM_data[(PheStart_dex+1):(PheStop_dex-1)]
        semicolon_dex<-grep(";",PheData,fixed = TRUE)
        StringNum<-semicolon_dex[2]-semicolon_dex[1]
        TraitNum<-as.numeric(MCIM_data[grep("_TraitNumber",MCIM_data,fixed = TRUE)+1])
        PheData<-matrix(PheData,ncol=StringNum,byrow=T)
        pheRaw<-PheData[,(StringNum-TraitNum):(StringNum-1),drop=F]
        traitName<-pheRaw[1,]
        pheRaw<-pheRaw[-1,,drop=F]
        # pheRaw<-as.data.frame(pheRaw)
        # pheRaw<-apply(pheRaw,2,as.numeric)
        EnvRaw<-PheData[-1,1,drop=F]
        EnvName<-unique(EnvRaw)
        EnvNum<-length(EnvName)

        pheRaw_All<-matrix(NA,nrow=indivNum*EnvNum,ncol=dim(pheRaw)[2])
        if(StringNum>(3+TraitNum)){
          if(MCIM_data[grep("_Replications",MCIM_data,fixed = TRUE)+1]=="yes"){
            newPhe<-NULL
            newGeo<-NULL
            newEN<-NULL
            RepRaw<-PheData[-1,2,drop=F]
            for(i in 1:EnvNum){
              now_id<-which(EnvRaw==EnvName[i])
              now_Rep<-RepRaw[now_id]
              ns<-unique(now_Rep)
              now_Geo<-as.numeric(PheData[-1,,drop=F][now_id,(StringNum-TraitNum-1)])
              now_Phe<-pheRaw[now_id,,drop=F]
              if(length(ns)>1){
                nr<-sort(unique(now_Geo))
                newEN<-c(newEN,rep(i,length(nr)))
                for(j in 1:length(nr)){
                  now_id2<-which(now_Geo==nr[j])
                  meanphe<-now_Phe[now_id2,,drop=F]
                  if(dim(meanphe)[1]>1){
                    meanphe<-apply(meanphe,2,as.numeric)
                    meanphe<-apply(meanphe,2,function(x)mean(x,na.rm=T))
                    meanphe<-matrix(meanphe,nrow=1)
                    newPhe<-rbind(newPhe, meanphe)
                    newGeo<-rbind(newGeo,nr[j])
                  }else{
                    # meanphe<-apply(meanphe,2,as.numeric)
                    # meanphe<-apply(meanphe,2,function(x)mean(x,na.rm=T))
                    # meanphe<-matrix(meanphe,nrow=1)
                    newPhe<-rbind(newPhe, meanphe)
                    newGeo<-rbind(newGeo,nr[j])
                  }

                }
              }else{
                newEN<-c(newEN,rep(i,length(now_Geo)))
                newPhe<-rbind(newPhe,now_Phe)
                newGeo<-rbind(newGeo,matrix(now_Geo,ncol=1))
              }
            }
          }else{
            warning("please check the phenotype or the parameter (_Replications) in file!")
          }

        }else{
          newPhe<-pheRaw
          newGeo<-as.numeric(PheData[-1,2])
          newEN<-NULL
          for(i in 1:EnvNum){

            now_id<-which(EnvRaw==EnvName[i])
            newEN<-c(newEN,rep(i,length(now_id)))
          }
        }
        newID<-(newEN-1)*indivNum+newGeo
        # pheRaw_All[newID,,drop=F]<-newPhe
        pheRaw_All[newID,]<-newPhe

        if(EnvNum>1){
          pheRaw<-matrix(pheRaw_All,nrow=indivNum,byrow=FALSE)
        }
        pheRaw[which(is.na(pheRaw))]<-"-"
        pheRaw[pheRaw=="."]<-"-"


        aa<-NULL
        for(i in 1:length(traitName)){
          aa<-c(aa,rep(traitName[i],EnvNum))
        }
        traitName<-aa

        colnames(pheRaw)<-traitName;rownames(pheRaw)<-NULL# dim(pheRaw)


        if(is.null(cov_en)==FALSE){
          covnum<-t(cov_en[-1,2:ncol(cov_en)])
          yygg1<-numeric()
          for(i in 1:nrow(covnum)){
            otrait_ind<-unique(covnum[i,])
            cov_col<-length(otrait_ind)-1
            col_each<-numeric()
            for(j in 1:length(covnum[i,])){
              if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                cov_0<-matrix(-1,1,cov_col)
              }else{
                cov_0<-matrix(0,1,cov_col)
                covnum_loc<-which(otrait_ind[]==covnum[i,j])
                cov_0[1,covnum_loc]<-1
              }
              col_each<-rbind(col_each,cov_0)

            }
            yygg1<-cbind(yygg1,col_each)# dim(yygg1)
          }
          cov_en<-covnum
        }else{
          yygg1<-NULL
        }


      }else if(fileFormat=="Cart"){
        ################## fileFormat=="Cart" ################
        start_dex<-grep("-start",y_jun3,fixed = TRUE)
        stop_dex<-grep("-stop",y_jun3,fixed = TRUE)
        chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)

        chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome")
        chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)

        chrdata_dexlen<-as.numeric(y_jun3[grep("-chromosomes",y_jun3,fixed = TRUE)+1])

        chr_num<-numeric()
        chr_numfirst<-numeric()
        markername0<-numeric()
        chr_pos<-numeric()

        for(i in 1:chrdata_dexlen){
          chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
          marker_name<-numeric()
          marker_pos<-numeric()

          for(j in 0:(trunc(length(chr_each)/2)-1) ){# trunc()
            marker_name<-cbind(marker_name,chr_each[2*j+1])
            marker_pos<-cbind(marker_pos,suppressWarnings(as.numeric(chr_each[2*(j+1)])))

            if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){
              marker_posm<-100*((-0.5)*log(1-2*marker_pos))
              markerlen<-length(marker_posm)
              marker_pos1<-marker_posm[1:(markerlen-1)]
              marker_pos2<-c(0,marker_pos1)
              marker_possum<-cumsum(marker_pos2)
            }
            if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
              marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))
              markerlen<-length(marker_posm)
              marker_pos1<-marker_posm[1:(markerlen-1)]
              marker_pos2<-c(0,marker_pos1)
              marker_possum<-cumsum(marker_pos2)
            }
            if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
              # y_jun3[9]=="M"  Morgan. The distance over which, on average, one crossover occurs per meiosis.
              marker_posm<-100*marker_pos
              markerlen<-length(marker_posm)
              marker_pos1<-marker_posm[1:(markerlen-1)]
              marker_pos2<-c(0,marker_pos1)
              marker_possum<-cumsum(marker_pos2)
            }
            if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
              marker_posm<-marker_pos
              markerlen<-length(marker_posm)
              marker_pos1<-marker_posm[1:(markerlen-1)]
              marker_pos2<-c(0,marker_pos1)
              marker_possum<-cumsum(marker_pos2)
            }
            if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
              marker_possum<-100*marker_pos
            }
            if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
              marker_possum<-marker_pos
            }
          }
          markername0<-cbind(markername0,marker_name)# length(markername0) class(markername0)
          chr_pos<-rbind(chr_pos,matrix(marker_possum,,1))
          chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
          chr_data<-na.omit(chr_a)
          chr_num<-rbind(chr_num,length(chr_data))

          chr_numfirst<-rbind(chr_numfirst,matrix(rep(i,chr_num[i]),,1))
        }

        chr_leng<-length(chr_pos)
        mapRaw1<-data.frame(marker=matrix(markername0,,1),chr=chr_numfirst,pos=chr_pos)# dim(mapRaw1)

        marker_dex<-grep("markers",y_jun3,fixed = TRUE)
        marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]

        indi_num<-as.numeric(y_jun3[grep("-SampleSize",y_jun3,fixed = TRUE)+1])

        if_indi<-y_jun3[marker_dex[1]-1]
        if(if_indi=="individuals"){
          genRaw<-matrix(marker_snp[-c(1:3,length(marker_snp))],nrow=indi_num,byrow = T)
          genRaw<-genRaw[,-1,drop=F]
          genRaw<-t(genRaw)# dim(genRaw)

        }else{
          genRaw<-matrix(marker_snp[-c(1:2,length(marker_snp))],ncol=indi_num+1,byrow = T)
          genRaw<-genRaw[,-1,drop=F]# dim(genRaw)
        }# dim(genRaw)
        rownames(genRaw)<-NULL;colnames(genRaw)<-NULL# class(genRaw)
        genRaw<-apply(genRaw,2,as.character)# class(genRaw)

        trait_total<-y_jun3[start_dex[3]:stop_dex[3]]#
        trait_total[which(trait_total==".")]<-"-"
        trait_total[which(trait_total=="*")]<-"-"


        trait_dex<-grep("traits",trait_total)
        iftrait_indi<-trait_total[trait_dex[1]-1]
        if(iftrait_indi=="individuals"){

          name_dex<-grep("named",trait_total)
          traitdata<-trait_total[(name_dex+1):(length(trait_total)-1)]
          pheRaw<-matrix(traitdata,nrow=indi_num,byrow = T)
          pheRaw<-pheRaw[,-1,drop=F]# dim(pheRaw)
          trait_num<-dim(pheRaw)[2]

          traitName<-trait_total[(trait_dex+2):(name_dex-1)]
          EnvNum<-0;TraitNum<-0
          for(i in length(traitName)){
            traitStr<-strsplit(traitName[i],"E")[[1]]
            EnvNum<-max(as.numeric(traitStr[2]),EnvNum)
            TraitNum<-max(as.numeric(strsplit(traitStr[1],"t")[[1]][2]),TraitNum)
          }
          newtraitName<-NULL
          for(i in 1:TraitNum){
            for(j in 1:EnvNum){
              newtraitName<-c(newtraitName,paste("t",i,"E",j,sep=""))
            }
          }
          NameID<-match(newtraitName,traitName)
          pheRaw<-pheRaw[,NameID]# dim(pheRaw)
          traitName<-traitName[NameID]
          # rownames(pheRaw)<-NULL;colnames(pheRaw)<-NULL

        }else{
          traitdata<-trait_total[(trait_dex[1]+1):(length(trait_total)-1)]
          pheRaw<-matrix(traitdata,ncol=indi_num+1,byrow = T)
          traitName<-c(pheRaw[,1])
          pheRaw<-t(pheRaw[,-1,drop=F])
          trait_num<-dim(pheRaw)[2]

          EnvNum<-0;TraitNum<-0
          for(i in length(traitName)){
            traitStr<-strsplit(traitName[i],"E")[[1]]
            EnvNum<-max(as.numeric(traitStr[2]),EnvNum)
            TraitNum<-max(as.numeric(strsplit(traitStr[1],"t")[[1]][2]),TraitNum)
          }
          newtraitName<-NULL
          for(i in 1:TraitNum){
            for(j in 1:EnvNum){
              newtraitName<-c(newtraitName,paste("t",i,"E",j,sep=""))
            }
          }
          NameID<-match(newtraitName,traitName)
          pheRaw<-pheRaw[,NameID]# dim(pheRaw)
          traitName<-traitName[NameID]
        }
        rownames(pheRaw)<-NULL;colnames(pheRaw)<-traitName

        if(length(start_dex)==3){yygg1<-NULL}
        if(length(start_dex)==4){
          if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
            cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
            cov_dex<-grep("otraits",cov_total)
            cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
            bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
            otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
            if(y_jun3[start_dex[4]+1]=="otraits"){

              covnumonly<-numeric()
              for( i in 0:(otrait_indi-1)){
                cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
                covnumonly<-rbind(covnumonly,cov_each)
              }
              covnum<-covnumonly[,-1]
              yygg1<-numeric()
              for(i in 1:nrow(covnum)){

                otrait_ind<-unique(covnum[i,])
                cov_col<-length(otrait_ind)-1

                col_each<-numeric()
                for(j in 1:length(covnum[i,])){

                  if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                    cov_0<-matrix(-1,1,cov_col)

                  }else{
                    cov_0<-matrix(0,1,cov_col)
                    covnum_loc<-which(otrait_ind[]==covnum[i,j])
                    cov_0[1,covnum_loc]<-1
                  }
                  col_each<-rbind(col_each,cov_0)

                }
                yygg1<-cbind(yygg1,col_each)
              }
            }
            if(y_jun3[start_dex[4]+1]=="individuals"){
              named_dex<-grep("named",cov_only)
              covdata<-cov_only[(named_dex+1):length(cov_only)]
              covnum<-numeric()
              # otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)

              for(n in 0:(indi_num-1)){
                cov_each<-matrix(covdata[(2+(otrait_indi+1)*n):(1+(otrait_indi+1)*n+otrait_indi)],otrait_indi,1)
                covnum<-cbind(covnum,cov_each)
              }
              yygg1<-numeric()
              for(i in 1:nrow(covnum)){

                otrait_ind<-unique(covnum[i,])
                cov_col<-length(otrait_ind)-1

                col_each<-numeric()
                for(j in 1:length(covnum[i,])){

                  if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                    cov_0<-matrix(-1,1,cov_col)

                  }else{
                    cov_0<-matrix(0,1,cov_col)
                    covnum_loc<-which(otrait_ind[]==covnum[i,j])
                    cov_0[1,covnum_loc]<-1
                  }
                  col_each<-rbind(col_each,cov_0)

                }
                yygg1<-cbind(yygg1,col_each)

              }
            }

          }
        }

        if(Population=="BC1"){
          genRaw<-gsub("1","-1",genRaw)
          genRaw<-gsub("2","1",genRaw)
          genRaw<-gsub("*","99",genRaw,fixed = TRUE)# dim(genRaw)
          # flagRIL<-0
        }else if(Population=="BC2"){
          genRaw<-gsub("0","-1",genRaw)
          #marker_snp2<-gsub("1","1",marker_snp1)
          genRaw<-gsub("*","99",genRaw,fixed = TRUE)
          # flagRIL<-0
        }else if(Population=="DH"){
          genRaw<-gsub("0","-1",genRaw)
          genRaw<-gsub("2","1",genRaw)
          genRaw<-gsub("*","99",genRaw,fixed = TRUE)
          # flagRIL<-0
        }else if(Population=="RIL"){
          genRaw<-gsub("1","99",genRaw)
          genRaw<-gsub("0","-1",genRaw)
          genRaw<-gsub("2","1",genRaw)
          genRaw<-gsub("*","99",genRaw,fixed = TRUE)
          # flagRIL<-1
        }else if(Population=="F2"){
          genRaw<-gsub("12","D",genRaw)
          genRaw<-gsub("10","C",genRaw)
          genRaw<-gsub("0","B",genRaw)
          genRaw<-gsub("2","A",genRaw)
          genRaw<-gsub("-1","-",genRaw)
          genRaw<-gsub("1","H",genRaw,fixed = TRUE)
          # flagRIL<-0
        }

        # seq_indi3<-c("covariate",seq_indi)
        # seq_indi3<-matrix(seq_indi3,,1)
        if(is.null(yygg1)==FALSE){
          cov_en<-covnum
        }else{
          cov_en<-NULL
        }

      }else if(fileFormat=="GCIM"){
        ################## fileFormat=="GCIM" ################
        mapRaw1<-data.frame(marker=mapRaw11[,1],chr=mapRaw11[,2],pos=mapRaw11[,3])# dim(mapRaw1)
        pheRaw<-as.matrix(pheRaw)
        TraitName<-pheRaw[1,-1]
        TraitNum<-length(TraitName)# Trait_Name<-unique(pheRaw[2,-1])
        OrigName<-pheRaw[2,-1]
        pheRaw<-pheRaw[3:dim(pheRaw)[1],-1]# dim(pheRaw)
        pheRaw[which(is.na(pheRaw))]<-"-"
        pheRaw<-as.matrix(pheRaw)
        colnames(pheRaw)<-OrigName;rownames(pheRaw)<-NULL# dim(pheRaw)
        genRaw<-as.matrix(genRaw1)

        if(is.null(cov_en)==FALSE){
          covnum<-t(cov_en)
          yygg1<-numeric()
          for(i in 1:nrow(covnum)){
            otrait_ind<-unique(covnum[i,])
            cov_col<-length(otrait_ind)-1
            col_each<-numeric()
            for(j in 1:length(covnum[i,])){
              if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                cov_0<-matrix(-1,1,cov_col)
              }else{
                cov_0<-matrix(0,1,cov_col)
                covnum_loc<-which(otrait_ind[]==covnum[i,j])
                cov_0[1,covnum_loc]<-1
              }
              col_each<-rbind(col_each,cov_0)

            }
            yygg1<-cbind(yygg1,col_each)
          }
        }else{
          yygg1<-NULL
        }
      }
      result<-list(pheRaw=pheRaw,genRaw=genRaw,mapRaw1=mapRaw1,yygg1=yygg1)# cov_en=cov_en,
      # dim(pheRaw)  dim(genRaw)  dim(mapRaw1)  flagRIL  dim(yygg1)  dim(cov_en)
    }
  }else if(method=="GCIM"){
    pheRaw=NULL;genRaw=NULL;mapRaw1=NULL;flag=NULL;yygg1=NULL;cov_en=NULL
    geoo<-readraw$geoo;pho<-readraw$pho;poss<-readraw$poss;parm<-readraw$parm;y_jun3<-readraw$y_jun3
    genRaw1<-readraw$genRaw1;pheRaw<-readraw$pheRaw;mapRaw11<-readraw$mapRaw11;cov_en<-readraw$cov_en

    if(fileFormat=="ICIM"){
      if(parm[4,1]==1){
        pos.be<-numeric()
        for(i in 1:10){
          pos1<-poss[which(poss[,2]==i),]
          poss1<-pos1
          positi<-as.matrix(cumsum(poss1[,3]))
          chrr<-as.matrix(poss1[,1:2])
          poss2<-cbind(chrr,positi)
          pos.be<-rbind(pos.be,poss2)
        }
      }
      if(parm[4,1]==2){
        pos.be<-poss
      }
      if(parm[5,1]==2){
        posthree<-matrix(100*pos.be[,3],,1)
        postwo<-pos.be[,1:2]
        pos<-cbind(postwo,posthree)
      }
      if(parm[5,1]==1){
        pos<-as.matrix(pos.be)
      }
      pos<-as.matrix(pos)
      geo<-geoo

      if(Model=="Random"){
        flag<-1
      }else if(Model=="Fixed"){
        flag<-0
      }
      gen_0<-geo[,-1]
      gen_0<-as.data.frame(gen_0,stringsAsFactors = F)
      gen_0<-sapply(gen_0,as.numeric)
      gen_0<-as.matrix(gen_0)
      gen_0<-matrix(as.character(gen_0),nrow(gen_0),ncol(gen_0))

      if(Population=="BC1"){
        gen_1<-gsub("-1","99",gen_0)
        gen_2<-gsub("1","-1",gen_1)
        gen_11<-gsub("2","1",gen_2)
        flagRIL<-0
      }else if(Population=="BC2"){
        gen_1<-gsub("-1","99",gen_0)
        #gen_2<-gsub("1","1",gen_1)
        gen_11<-gsub("0","-1",gen_1)
        flagRIL<-0
      }else if(Population=="DH"){
        gen_1<-gsub("-1","99",gen_0)
        gen_2<-gsub("0","-1",gen_1)
        gen_11<-gsub("2","1",gen_2)
        flagRIL<-0
      }else if(Population=="RIL"){
        gen_1<-gsub("-1","99",gen_0)
        gen_2<-gsub("0","-1",gen_1)
        gen_11<-gsub("2","1",gen_2)
        flagRIL<-1
      }else if(Population=="F2"){
        gen_1<-gsub("12","D",gen_0)
        gen_2<-gsub("10","C",gen_1)
        gen_3<-gsub("0","B",gen_2)
        gen_4<-gsub("2","A",gen_3)
        gen_5<-gsub("-1","-",gen_4)
        gen_11<-gsub("1","H",gen_5)
        flagRIL<-0
      }

      phett<-t(pho)
      phe_m<-as.matrix(phett[-1,])
      phe_00<-gsub(-100,NA,phe_m)

      seq_indiv<-seq(1,nrow(phe_00))
      seq_indiv1<-c("genotype",seq_indiv)
      seq_indiv1<-matrix(seq_indiv1,1,)
      geo1<-cbind(geo[,1],gen_11)
      genRaw<-rbind(seq_indiv1,geo1)

      seq_indiv2<-c("phenotype",seq_indiv)
      seq_indiv2<-matrix(seq_indiv2,,1)
      phename<-matrix(phett[1,],1,)
      phe<-rbind(phename,phe_00)
      pheRaw<-cbind(seq_indiv2,phe)

      colname_mapRaw1<-c("marker","chr","pos")
      colname_mapRaw1<-matrix(colname_mapRaw1,1,)
      mapRaw1<-rbind(colname_mapRaw1,pos)

      if(is.null(cov_en)==FALSE){
        cov_en1<-cov_en[-1,2:ncol(cov_en)]
        covnum<-t(cov_en1)
        yygg1<-numeric()
        for(i in 1:nrow(covnum)){
          otrait_ind<-unique(covnum[i,])
          cov_col<-length(otrait_ind)-1
          col_each<-numeric()
          for(j in 1:length(covnum[i,])){
            if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
              cov_0<-matrix(-1,1,cov_col)
            }else{
              cov_0<-matrix(0,1,cov_col)
              covnum_loc<-which(otrait_ind[]==covnum[i,j])
              cov_0[1,covnum_loc]<-1
            }
            col_each<-rbind(col_each,cov_0)

          }
          yygg1<-cbind(yygg1,col_each)
        }
      }else{
        yygg1<-NULL
      }

    }else if(fileFormat=="Cart"){
      if(Model=="Random"){
        flag<-1
      }else if(Model=="Fixed"){
        flag<-0
      }
      start_dex<-grep("-start",y_jun3,fixed = TRUE)
      stop_dex<-grep("-stop",y_jun3,fixed = TRUE)

      chr_dex<-grep("-Chromosome",y_jun3,fixed = TRUE)

      chrdata<-c(y_jun3[chr_dex[1]:(stop_dex[1]-1)],"-Chromosome",y_jun3[stop_dex[1]])
      chrdata_dex<-grep("-Chromosome",chrdata,fixed = TRUE)
      chrdata_dexlen<-length(chrdata_dex)

      chr_num<-numeric()
      chrname_num<-numeric()
      chr_numfirst<-numeric()
      markername0<-numeric()
      chr_pos<-numeric()
      chrRaw_name<-numeric()
      chr_Rawnumfirst<-numeric()
      for(i in 1:(chrdata_dexlen-1)){
        chr_each<-chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]
        marker_name<-numeric()
        marker_pos<-numeric()

        for(j in 0:(trunc(length(chr_each)/2)-1) ){
          marker_name0<-chr_each[2*j+1]
          marker_name<-cbind(marker_name,marker_name0)
          marker_pos0<-suppressWarnings(as.numeric(chr_each[2*(j+1)]))
          marker_pos<-cbind(marker_pos,marker_pos0)
          if((y_jun3[5]=="interval")&&(y_jun3[7]=="1")&&(y_jun3[9]=="r")){
            marker_posm<-100*((-0.5)*log(1-2*marker_pos))

            markerlen<-length(marker_posm)
            marker_pos1<-marker_posm[1:(markerlen-1)]
            marker_pos2<-c(0,marker_pos1)
            marker_possum<-cumsum(marker_pos2)
          }
          if((y_jun3[5]=="interval")&&(y_jun3[7]=="2")&&(y_jun3[9]=="r")){
            marker_posm<-100*(0.25*log((1+2*marker_pos)/(1-2*marker_pos)))

            markerlen<-length(marker_posm)
            marker_pos1<-marker_posm[1:(markerlen-1)]
            marker_pos2<-c(0,marker_pos1)
            marker_possum<-cumsum(marker_pos2)

          }
          if((y_jun3[5]=="interval")&&(y_jun3[9]=="M")){
            marker_posm<-100*marker_pos

            markerlen<-length(marker_posm)
            marker_pos1<-marker_posm[1:(markerlen-1)]
            marker_pos2<-c(0,marker_pos1)
            marker_possum<-cumsum(marker_pos2)
          }
          if((y_jun3[5]=="interval")&&(y_jun3[9]=="cM")){
            marker_posm<-marker_pos

            markerlen<-length(marker_posm)
            marker_pos1<-marker_posm[1:(markerlen-1)]
            marker_pos2<-c(0,marker_pos1)
            marker_possum<-cumsum(marker_pos2)
          }
          if((y_jun3[5]=="position")&&(y_jun3[9]=="M")){
            marker_possum<-100*marker_pos
          }
          if((y_jun3[5]=="position")&&(y_jun3[9]=="cM")){
            marker_possum<-marker_pos
          }
        }
        markername0<-cbind(markername0,marker_name)
        markername<-matrix(markername0,,1)
        marker_possum0<-matrix(marker_possum,,1)
        chr_pos<-rbind(chr_pos,marker_possum0)
        chr_a<-suppressWarnings(as.numeric(chrdata[(chrdata_dex[i]+2):(chrdata_dex[i+1]-1)]))
        chr_data<-na.omit(chr_a)
        chr_datalen<-length(chr_data)
        chr_num<-rbind(chr_num,chr_datalen)

        chrRawname<-chrdata[(chrdata_dex[i]+1)]
        chrname<-str_extract_all(chrRawname,"[0-9]+")
        chrRawname<-matrix(chrRawname,,1)
        chrRaw_name<-rbind(chrRaw_name,chrRawname)

        chrname_num<-rbind(chrname_num,chrname)
        chr_numxx<-rep(as.numeric(chrname_num[i]),chr_num[i])
        chr_numfirst<-rbind(chr_numfirst,matrix(chr_numxx,,1))

        chr_Rawnumxx<-rep(chrRaw_name[i],chr_num[i])
        chr_Rawnumfirst<-rbind(chr_Rawnumfirst,matrix(chr_Rawnumxx,,1))
      }

      chr_leng<-length(chr_pos)
      chr_numtwo<-cbind(chr_numfirst,chr_pos)
      marker_dex<-grep("markers",y_jun3,fixed = TRUE)
      marker_snp<-y_jun3[start_dex[2]:stop_dex[2]]

      marker_snpnum<-marker_snp
      snpa<-suppressWarnings(as.numeric(marker_snpnum))


      snpdata<-na.omit(snpa)
      indi_num<-length(snpdata)/chr_leng
      snp_data<-numeric()

      if_indi<-y_jun3[marker_dex[1]-1]
      if(if_indi=="individuals"){

        for(i in 0:(indi_num-1)){
          snp_eve<-matrix(snpdata[(chr_leng*i+1):(chr_leng*(i+1))],1,)
          snp_data<-rbind(snp_data,snp_eve)
        }
        snp_data1<-t(snp_data)
      }else{
        for(i in 0:(chr_leng-1)){
          snp_eve<-matrix(snpdata[(indi_num*i+1):(indi_num*(i+1))],1,)
          snp_data<-rbind(snp_data,snp_eve)
        }
        snp_data1<-snp_data
      }

      trait_total<-y_jun3[start_dex[3]:stop_dex[3]]

      for(i in 1:length(trait_total)){
        if(trait_total[i]=="."){
          trait_total[i]<-"0"
        }
      }
      trait_dex<-grep("traits",trait_total)
      traita<-suppressWarnings(as.numeric(trait_total))
      traitdata<-na.omit(traita)
      trait_num<-length(traitdata)/indi_num
      trait_data<-numeric()

      iftrait_indi<-trait_total[trait_dex[1]-1]
      if(iftrait_indi=="individuals"){

        for(i in 0:(indi_num-1)){
          trait_bbb<-traitdata[(trait_num*i+1):(trait_num*(i+1))]
          for(j in 1:length( trait_bbb)){
            if(trait_bbb[j]==0){
              trait_bbb[j]<-NA
            }
          }
          trait_eve<-matrix(trait_bbb,1,)
          trait_data<-rbind(trait_data,trait_eve)
        }

      }else{

        for(i in 0:(trait_num-1)){
          trait_aaa<-traitdata[(indi_num*i+1):(indi_num*(i+1))]
          for(j in 1:length(trait_aaa)){
            if(trait_aaa[j]==0){
              trait_aaa[j]<-NA
            }
          }
          trait_eve<-matrix(trait_aaa,,1)
          trait_data<-cbind(trait_data,trait_eve)
        }

      }
      if(length(start_dex)==3){yygg1<-NULL}
      if(length(start_dex)==4){
        if(y_jun3[start_dex[4]+1]=="otraits"||y_jun3[start_dex[4]+1]=="individuals"){
          cov_total<-y_jun3[start_dex[4]:stop_dex[4]]
          cov_dex<-grep("otraits",cov_total)
          cov_only<-y_jun3[(start_dex[4]+2):(stop_dex[4]-1)]
          bycross_dex<-grep("#bycross",y_jun3,fixed = TRUE)
          otrait_indi<-as.numeric(y_jun3[bycross_dex+8])
          if(y_jun3[start_dex[4]+1]=="otraits"){

            covnumonly<-numeric()
            for( i in 0:(otrait_indi-1)){
              cov_each<-matrix(cov_only[(i*(indi_num+1)+1):((indi_num+1)*(i+1))],1,)
              covnumonly<-rbind(covnumonly,cov_each)
            }
            covnum<-covnumonly[,-1]
            yygg1<-numeric()
            for(i in 1:nrow(covnum)){

              otrait_ind<-unique(covnum[i,])
              cov_col<-length(otrait_ind)-1

              col_each<-numeric()
              for(j in 1:length(covnum[i,])){

                if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                  cov_0<-matrix(-1,1,cov_col)

                }else{
                  cov_0<-matrix(0,1,cov_col)
                  covnum_loc<-which(otrait_ind[]==covnum[i,j])
                  cov_0[1,covnum_loc]<-1
                }
                col_each<-rbind(col_each,cov_0)

              }
              yygg1<-cbind(yygg1,col_each)

            }

          }
          if(y_jun3[start_dex[4]+1]=="individuals"){
            covdata<-cov_only[(2+otrait_indi):length(cov_only)]
            covnum<-numeric()
            otrait_name<-matrix(cov_only[2:2+(otrait_indi-1)],otrait_indi,1)
            for(m in 1:otrait_indi){

              coveach<-numeric()
              for(n in 0:(indi_num-1)){
                cov_each<-matrix(covdata[m+otrait_indi*n],1,1)

                coveach<-cbind(coveach,cov_each)
              }
              covnum<-rbind(covnum,coveach)
            }
            covnumonly<-cbind(otrait_name,covnum)
            yygg1<-numeric()
            for(i in 1:nrow(covnum)){

              otrait_ind<-unique(covnum[i,])
              cov_col<-length(otrait_ind)-1

              col_each<-numeric()
              for(j in 1:length(covnum[i,])){

                if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
                  cov_0<-matrix(-1,1,cov_col)

                }else{
                  cov_0<-matrix(0,1,cov_col)
                  covnum_loc<-which(otrait_ind[]==covnum[i,j])
                  cov_0[1,covnum_loc]<-1
                }
                col_each<-rbind(col_each,cov_0)

              }
              yygg1<-cbind(yygg1,col_each)

            }
          }

        }
      }

      seq_indi<-seq(1,nrow(trait_data))
      seq_indi1<-c("genotype",seq_indi)
      seq_indi1<-matrix(seq_indi1,1,)
      snp1<-cbind(markername,snp_data1)
      seq_indi2<-c("phenotype",seq_indi)
      seq_indi2<-matrix(seq_indi2,,1)
      num_trait<-ncol(trait_data)
      seq_trait<-seq(1,num_trait)
      seq_trait<-matrix(seq_trait,1,)
      trait_data00<-rbind(seq_trait,trait_data)

      colnames_mapname<-c("marker","chr","pos")
      colnames_mapname<-matrix(colnames_mapname,1,)
      mapRaw1<-cbind(markername,chr_Rawnumfirst,chr_pos)

      mapRaw1<-rbind(colnames_mapname,mapRaw1)
      genRawqq<-rbind(seq_indi1,snp1)
      genRawq<-genRawqq[-1,]

      if(Population=="BC1"){
        marker_snp1<-gsub("1","-1",genRawq)
        marker_snp2<-gsub("2","1",marker_snp1)
        genRawh<-gsub("*","99",marker_snp2,fixed = TRUE)
        genRaw<-rbind(genRawqq[1,],genRawh)
        flagRIL<-0
      }else if(Population=="BC2"){
        marker_snp1<-gsub("0","-1",genRawq)
        #marker_snp2<-gsub("1","1",marker_snp1)
        genRawh<-gsub("*","99",marker_snp1,fixed = TRUE)
        genRaw<-rbind(genRawqq[1,],genRawh)
        flagRIL<-0
      }else if(Population=="DH"){
        marker_snp1<-gsub("0","-1",genRawq)
        marker_snp2<-gsub("2","1",marker_snp1)
        genRawh<-gsub("*","99",marker_snp2,fixed = TRUE)
        genRaw<-rbind(genRawqq[1,],genRawh)
        flagRIL<-0
      }else if(Population=="RIL"){
        marker_snp10<-gsub("1","99",genRawq)
        marker_snp1<-gsub("0","-1",marker_snp10)
        marker_snp2<-gsub("2","1",marker_snp1)
        genRawh<-gsub("*","99",marker_snp2,fixed = TRUE)
        genRaw<-rbind(genRawqq[1,],genRawh)
        flagRIL<-1
      }else if(Population=="F2"){
        marker_snp1<-gsub("12","D",genRawq)
        marker_snp2<-gsub("10","C",marker_snp1)
        marker_snp3<-gsub("0","B",marker_snp2)
        marker_snp4<-gsub("2","A",marker_snp3)
        marker_snp5<-gsub("-1","-",marker_snp4)
        genRawh<-gsub("1","H",marker_snp5,fixed = TRUE)
        genRaw<-rbind(genRawqq[1,],genRawh)
        flagRIL<-0
      }
      pheRaw<-cbind(seq_indi2,trait_data00)

      seq_indi3<-c("covariate",seq_indi)
      seq_indi3<-matrix(seq_indi3,,1)
      if(is.null(yygg1)==FALSE){
        cov_en<-cbind(seq_indi3,t(covnumonly))
      }else{
        cov_en<-NULL
      }
    }else if(fileFormat=="GCIM"){
      genRaw1qq<-as.matrix(genRaw1)
      genRaw1<-genRaw1qq[-1,-1]
      pheRaw<-as.matrix(pheRaw)
      mapRaw1<-as.matrix(mapRaw11)

      if(Population=="BC1"){
        genRaw_<-gsub("-","99",genRaw1)
        genRaw_Aa<-gsub("H","-1",genRaw_)
        genRaw_AA<-gsub("A","1",genRaw_Aa)

        genRaw<-rbind(genRaw1qq[1,],cbind(genRaw1qq[-1,1],genRaw_AA))
        flagRIL<-0
      }else if(Population=="BC2"){
        genRaw_<-gsub("-","99",genRaw1)
        genRaw_aa<-gsub("B","-1",genRaw_)
        genRaw_Aa<-gsub("H","1",genRaw_aa)
        genRaw<-rbind(genRaw1qq[1,],cbind(genRaw1qq[-1,1],genRaw_Aa))
        flagRIL<-0
      }else if(Population=="DH"){
        genRaw_<-gsub("-","99",genRaw1)
        genRaw_aa<-gsub("B","-1",genRaw_)
        genRaw_AA<-gsub("A","1",genRaw_aa)

        genRaw<-rbind(genRaw1qq[1,],cbind(genRaw1qq[-1,1],genRaw_AA))
        flagRIL<-0
      }else if(Population=="RIL"){
        genRaw_<-gsub("-","99",genRaw1)
        genRaw_aa<-gsub("B","-1",genRaw_)
        genRaw_AA<-gsub("A","1",genRaw_aa)
        genRaw<-rbind(genRaw1qq[1,],cbind(genRaw1qq[-1,1],genRaw_AA))
        flagRIL<-1
      }else if(Population=="F2"){
        genRaw<-genRaw1qq
        flagRIL<-0
      }

      if(Model=="Random"){
        flag<-1
      }else if(Model=="Fixed"){
        flag<-0
      }


      if(is.null(cov_en)==FALSE){
        cov_en1<-cov_en[-1,2:ncol(cov_en)]
        covnum<-t(cov_en1)
        yygg1<-numeric()
        for(i in 1:nrow(covnum)){
          otrait_ind<-unique(covnum[i,])
          cov_col<-length(otrait_ind)-1
          col_each<-numeric()
          for(j in 1:length(covnum[i,])){
            if(covnum[i,j]==otrait_ind[length(otrait_ind)]){
              cov_0<-matrix(-1,1,cov_col)
            }else{
              cov_0<-matrix(0,1,cov_col)
              covnum_loc<-which(otrait_ind[]==covnum[i,j])
              cov_0[1,covnum_loc]<-1
            }
            col_each<-rbind(col_each,cov_0)

          }
          yygg1<-cbind(yygg1,col_each)
        }
      }else{
        yygg1<-NULL
      }
    }
    result<-list(pheRaw=pheRaw,genRaw=genRaw,mapRaw1=mapRaw1,flag=flag,flagRIL=flagRIL,yygg1=yygg1,cov_en=cov_en)
  }else{
    warning("Please specify the method GCIM or GCIM-QEI!")

  }



  return(result)
}
