#' Read raw data
#'
#' @param file Dataset input
#' @param fileFormat Format of dataset.
#' @param method Method "GCIM" or method "GCIM-QEI"
#' @param filecov Covariate file of QTLIciMapping or QTLNetwork.
#' @param MCIMmap Map file of QTLNetwork.
#' @param MultiEnv Whether to perform multi-environment analysis
#'
#' @return a list
#' @export
#'
#' @examples
#' data(F2data)
#' Readdata(file=F2data,fileFormat="GCIM",
#' method="GCIM-QEI",filecov=NULL,
#' MCIMmap=NULL,MultiEnv=TRUE)
Readdata<-function(file=NULL,fileFormat=NULL,method=NULL,filecov=NULL,MCIMmap=NULL,MultiEnv=FALSE){
  if(method=="GCIM-QEI"){
    geoo=NULL;pho=NULL;poss=NULL;parm=NULL;
    MCIM_data=NULL;map_data=NULL;
    y_jun3=NULL;
    genRaw1=NULL;pheRaw=NULL;mapRaw11=NULL;
    cov_en=NULL

    if(MultiEnv==TRUE){
      if(fileFormat=="ICIM"){
        ################## fileFormat=="ICIM" ##################
        suppressMessages(geoo<-as.matrix(read_excel(file,sheet = "Genotype",col_names = FALSE)))
        suppressMessages(pho<-as.matrix(read_excel(file,sheet = "Phenotype",col_names = FALSE)))
        suppressMessages(poss<-read_excel(file,sheet = "LinkageMap",col_names = FALSE))
        suppressMessages(parm<-read_excel(file,sheet = "GeneralInfo",col_names = FALSE))

        if(is.null(filecov)==FALSE){
          cov_en<-fread(filecov,header = FALSE,stringsAsFactors=T)
          cov_en<-as.matrix(cov_en)
        }
      }else if(fileFormat=="MCIM"){
        ################## fileFormat=="MCIM" ##################
        MCIM_data<-scan(file,what="",sep="")
        # read.table(file,sep=",",header=T,stringsAsFactors=F)
        if(length(MCIMmap)>0){
          map_data<- scan(MCIMmap,what = "",sep = "")
        }else{
          warning("Please enter MCIMmap parameter!")
        }
        # read.table(MCIMmap,sep=",",header=T,stringsAsFactors=F)
        if(is.null(filecov)==FALSE){
          cov_en<-fread(filecov,header = FALSE,stringsAsFactors=T)
          cov_en<-as.matrix(cov_en)
        }
      }else if(fileFormat=="Cart"){
        ################## fileFormat=="Cart" ##################
        y_jun3<-scan(file,what = "",sep = "")
      }else if(fileFormat=="GCIM"){
        ################## fileFormat=="GCIM" ##################
        if(is.character(file)==TRUE){
          genRaw<-fread(file,header = T,stringsAsFactors=T)
          titlenameGen<-colnames(genRaw)[1:3]
        }else{
          genRaw<-file
          titlenameGen<-genRaw[1,1:3]
        }

        hapName<-c("marker","chr","pos")
        if(all(titlenameGen==hapName)==FALSE){
          warning("please check the Linkage map's name in file!")
        }

        if(is.character(file)==FALSE){
          genRaw_colname<-genRaw[1,]
          genRaw<-genRaw[-1,]
          colnames(genRaw)<-as.character(genRaw_colname)
        }

        traitloc<-which((genRaw[,1]=="Env1")|(genRaw[,1]=="env1"))
        if(length(traitloc)==0){
          warning("please check the phenotype in file!")
        }

        envirloc<-which((genRaw[,2]=="Covar1")|(genRaw[,2]=="covar1"))#
        indv_name<-colnames(genRaw)
        indv_name[1:3]<-c("Env:","Trait:","TraitName:")
        if(length(envirloc)!=0){
          pheRaw<-cbind(indv_name,t(genRaw[traitloc[1]:(envirloc-1),]))

          cov_en<-t(genRaw[envirloc:nrow(genRaw),-c(1:3)])
          colnames(cov_en)<-NULL;rownames(cov_en)<-NULL
        }else{
          pheRaw<-cbind(indv_name,t(genRaw[traitloc[1]:nrow(genRaw),]))
        }

        pheRaw<-as.matrix(pheRaw)
        colnames(pheRaw)<-NULL;rownames(pheRaw)<-NULL# dim(pheRaw)
        genRaw1<-as.matrix(genRaw[c(1:(traitloc[1]-1)),-c(1:3)])
        colnames(genRaw1)<-NULL;rownames(genRaw1)<-NULL# dim(genRaw1)
        mapRaw11<-as.matrix(genRaw[1:(traitloc[1]-1),1:3])
        colnames(mapRaw11)<-NULL;rownames(mapRaw11)<-NULL# dim(mapRaw11)
      }
      result<-list(geoo=geoo,pho=pho,poss=poss,parm=parm,
                   MCIM_data=MCIM_data,map_data=map_data,
                   y_jun3=y_jun3,
                   genRaw1=genRaw1,pheRaw=pheRaw,mapRaw11=mapRaw11,cov_en=cov_en)
    }else{
      if(fileFormat=="ICIM"){
        ################## fileFormat=="ICIM" ##################
        suppressMessages(geoo<-as.matrix(read_excel(file,sheet = "Genotype",col_names = FALSE)))
        suppressMessages(pho<-as.matrix(read_excel(file,sheet = "Phenotype",col_names = FALSE)))
        suppressMessages(poss<-read_excel(file,sheet = "LinkageMap",col_names = FALSE))
        suppressMessages(parm<-read_excel(file,sheet = "GeneralInfo",col_names = FALSE))
        if(is.null(filecov)==FALSE){
          cov_en<-fread(filecov,header = FALSE,stringsAsFactors=T)
          cov_en<-as.matrix(cov_en)
        }
      }else if(fileFormat=="MCIM"){
        ################## fileFormat=="MCIM" ##################
        MCIM_data<-scan(file,what = "",sep = "")
        if(length(MCIMmap)>0){
          map_data<- scan(MCIMmap,what = "",sep = "")
        }else{
          warning("Please enter MCIMmap parameter!")
        }
        if(is.null(filecov)==FALSE){
          cov_en<-fread(filecov,header = FALSE,stringsAsFactors=T)
          cov_en<-as.matrix(cov_en)
        }
      }else if(fileFormat=="Cart"){
        ################## fileFormat=="Cart" ##################
        y_jun3<-scan(file,what = "",sep = "")
      }else if(fileFormat=="GCIM"){
        ################## fileFormat=="GCIM" ##################
        if(is.character(file)==TRUE){
          genRaw<-fread(file,header = T,stringsAsFactors=T)
          titlenameGen<-colnames(genRaw)[1:3]
        }else{
          genRaw<-file
          titlenameGen<-genRaw[1,1:3]
        }

        hapName<-c("marker","chr","pos")
        if(all(titlenameGen==hapName)==FALSE){
          warning("please check the Linkage map's name in file!")
        }

        if(is.character(file)==FALSE){
          genRaw_colname<-genRaw[1,]
          genRaw<-genRaw[-1,]
          colnames(genRaw)<-as.character(genRaw_colname)
        }

        traitloc<-which((genRaw[,2]=="trait1")|(genRaw[,2]=="Trait1"))
        if(length(traitloc)==0){
          warning("please check the phenotype in file!")
        }

        envirloc<-which((genRaw[,2]=="Covar1")|(genRaw[,2]=="covar1"))#
        indv_name<-c("Trait:","TraitName:",colnames(genRaw)[-c(1:3)])
        if(length(envirloc)!=0){
          pheRaw<-cbind(indv_name,t(genRaw[traitloc[1]:(envirloc-1),-1,drop=F]))

          cov_en<-t(genRaw[envirloc:nrow(genRaw),-c(1:3)])
          colnames(cov_en)<-NULL;rownames(cov_en)<-NULL
        }else{
          pheRaw<-cbind(indv_name,t(genRaw[traitloc[1]:nrow(genRaw),-1,drop=F]))
        }

        pheRaw<-as.matrix(pheRaw)
        colnames(pheRaw)<-NULL;rownames(pheRaw)<-NULL# dim(pheRaw)
        genRaw1<-as.matrix(genRaw[c(1:(traitloc[1]-1)),-c(1:3)])
        colnames(genRaw1)<-NULL;rownames(genRaw1)<-NULL# dim(genRaw1)
        mapRaw11<-as.matrix(genRaw[1:(traitloc[1]-1),1:3])
        colnames(mapRaw11)<-NULL;rownames(mapRaw11)<-NULL# dim(mapRaw11)

      }

      result<-list(geoo=geoo,pho=pho,poss=poss,parm=parm,
                   MCIM_data=MCIM_data,map_data=map_data,
                   y_jun3=y_jun3,
                   genRaw1=genRaw1,pheRaw=pheRaw,mapRaw11=mapRaw11,cov_en=cov_en)#
    }

  }else if(method=="GCIM"){
    geoo=NULL;pho=NULL;poss=NULL;parm=NULL;y_jun3=NULL;genRaw1=NULL;pheRaw=NULL;mapRaw11=NULL;cov_en=NULL
    if(fileFormat=="ICIM"){
      geoo<-as.matrix(read.xlsx(file,sheet = "Genotype",colNames = FALSE))
      pho<-as.matrix(read.xlsx(file,sheet = "Phenotype",colNames = FALSE))
      poss<-read.xlsx(file,sheet = "LinkageMap",colNames = FALSE)
      parm<-read.xlsx(file,sheet = "GeneralInfo",colNames = FALSE)
      if(is.null(filecov)==FALSE){
        cov_en<-fread(filecov,header = FALSE,stringsAsFactors=T)
        cov_en<-as.matrix(cov_en)
      }
    }else if(fileFormat=="Cart"){
      y_jun3<-scan(file,what = "",sep = "")

    }else if(fileFormat=="GCIM"){
      if(is.character(file)==TRUE){
        genRaw<-fread(file,header = FALSE,stringsAsFactors=T)
      }else{
        genRaw<-file
      }

      titlenameGen<-genRaw[1,1:3]
      hapName<-c("marker","chr","pos")

      if(all(titlenameGen==hapName)==FALSE){
        warning("please check the Linkage map's name in file!")
      }
      traitloc<-which(genRaw[,2]=="trait1")[1]

      if(length(traitloc)==0){
        warning("please check the phenotype in file!")
      }
      envirloc<-which(genRaw[,2]=="Covar1")
      if(length(envirloc)!=0){
        pheRaw<-t(rbind(genRaw[1,][,-(1:2)],genRaw[traitloc:(envirloc-1),][,-(1:2)]))
        cov_en<-cbind(t(genRaw[1,][,-(1:2)]),t(genRaw[envirloc:nrow(genRaw),][,-(1:2)]))
        colnames(cov_en)<-NULL;rownames(cov_en)<-NULL
      }else{
        pheRaw<-t(rbind(genRaw[1,][,-(1:2)],genRaw[traitloc:nrow(genRaw),][,-(1:2)]))
      }
      colnames(pheRaw)<-NULL;rownames(pheRaw)<-NULL
      genRaw1<-as.matrix(genRaw[1:(traitloc-1),-c(2,3)])
      mapRaw11<-as.matrix(genRaw[1:(traitloc-1),1:3])
    }
    result<-list(geoo=geoo,pho=pho,poss=poss,parm=parm,y_jun3=y_jun3,genRaw1=genRaw1,pheRaw=pheRaw,mapRaw11=mapRaw11,cov_en=cov_en)
  }else{
    warning("Please specify the method GCIM or GCIM-QEI!")
  }

  return(result)
}
