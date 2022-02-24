#' QTL Genome-Wide Composite Interval Mapping
#'
#' @param file File path and name in your computer.
#' @param fileFormat Format for input file: GCIM, ICIM, Cart, or MCIM.
#' @param filecov Covariate file of QTLIciMapping or QTLNetwork.
#' @param MCIMmap Map file of QTLNetwork.
#' @param Population Population type: F2, BC1, BC2, DH, RIL.
#' @param method Method "GCIM" or method "GCIM-QEI".
#' @param MultiEnv Whether to perform multi-environment analysis.
#' @param Model Random or fixed model.
#' @param WalkSpeed Walk speed for Genome-wide Scanning.
#' @param CriLOD Critical LOD scores for significant QTL.
#' @param CriDis The distance of optimization.
#' @param Likelihood This parameter is only for F2 population,including REML (restricted maximum likelihood ) and ML(maximum likelihood).
#' @param SetSeed In which the cross validation experiment is needed. Generally speaking, the random seed in the cross-validation experiment was set as 11001. If some known genes are not identified by the seed, users may try to use some new random seeds. At this case, one better result may be obtained.
#' @param flagrqtl This parameter is only for F2 population,   flagrqtl="FALSE" in the first running. If the other software detects only one QTL in a neighborhood but this software finds two linked QTLs (one with additive effect and another with dominant effect) in the region, let flagrqtl="TRUE"
#' @param DrawPlot This parameter is for all the populations, including FALSE and TRUE, DrawPlot=FALSE indicates no figure output, DrawPlot=TRUE indicates the output of the figure against genome position.
#' @param PlotFormat This parameter is for all the figure files, including *.jpeg, *.png, *.tiff and *.pdf.
#' @param Resolution This parameter is for all the figure files, including Low and High.
#' @param Trait Trait=1:3 indicates the analysis from the first trait to the third trait.
#' @param dir This parameter is for the save path.
#' @param CLO Number of CPUs.
#'
#' @export
#'
#' @examples
#' data(F2data)
#' QTL.gCIMapping(file=F2data,Population="F2",
#' MultiEnv=TRUE,Model="Random",CriLOD=3,
#' Trait=1,dir=tempdir(),CLO=2)
QTL.gCIMapping<-function(file=NULL,fileFormat="GCIM",filecov=NULL,MCIMmap=NULL,Population=NULL,
                         method="GCIM-QEI",MultiEnv=FALSE,Model="Random",
                         WalkSpeed=1,CriLOD=3,CriDis=5,Likelihood="REML",SetSeed=11001,flagrqtl=FALSE,
                         DrawPlot=TRUE,PlotFormat="jpeg",Resolution="Low",Trait=NULL,dir=NULL,CLO=NULL){
  # library(data.table)
  # library(qtl)
  # library(foreach)
  # library(doParallel)
  # library(lars)
  # library(MASS)
  # library(Rcpp)
  # library(readxl)
  # library(openxlsx)
  GCIM_QEI_plotF2<-function(PlotFormat,Resolution,OUTresult,moreMAP,P_Q,P_QE,CriLOD,NUM,MultiEnv){
    if(MultiEnv==TRUE){
      ###########################
      moreCHR<-as.numeric(moreMAP[,1])
      morePOS<-as.numeric(moreMAP[,2])

      original_Chr=unique(mapRaw1$chr)
      max_pos<-NULL
      chr_ID<-NULL;s<-0
      for(i in 1:max(moreCHR)){
        s<-s+length(which(moreCHR==i))
        chr_ID<-c(chr_ID,s)
        max_pos<-c(max_pos,max(morePOS[which(moreCHR==i)]))
      };rm(s)
      halfmax_pos<-max_pos*0.5
      u<-sum(max_pos)
      max_pos<-c(0,max_pos[-max(moreCHR)])
      sum_pos<-NULL;big_sumPOS<-NULL;s<-0
      for(i in 1:max(moreCHR)){
        s<-s+max_pos[i]
        sum_pos<-c(sum_pos,s)
        big_sumPOS<-c(big_sumPOS,morePOS[which(moreCHR==i)]+sum_pos[i])
      }
      halfsum_pos<-sum_pos+halfmax_pos

      ID_Q<-which(OUTresult$LOD_QTL>=CriLOD)
      ID_QE<-which(OUTresult$LOD_QEI>=CriLOD)
      ResultFinal_Q<-OUTresult[ID_Q,]
      ResultFinal_QE<-OUTresult[ID_QE,]

      Result_P_Q<-ResultFinal_Q$LOD_QTL
      Result_P_QE<-ResultFinal_QE$LOD_QEI
      if(dim(ResultFinal_Q)[1]>0){
        t<-NULL
        for(i in 1:dim(ResultFinal_Q)[1]){t<-c(t,sum_pos[ResultFinal_Q$Chr[i]])}
        ResultFinal_Q$pos_all<-ResultFinal_Q[,3]+t
      }
      if(dim(ResultFinal_QE)[1]>0){
        t<-NULL
        for(i in 1:dim(ResultFinal_QE)[1]){t<-c(t,sum_pos[ResultFinal_QE$Chr[i]])}
        ResultFinal_QE$pos_all<-ResultFinal_QE[,3]+t
      }
      #########################################
      margin_space<-1.5
      axis_space<-1
      unit_value<-"mm"
      pointsizeqqvalue<-30
      if(Resolution=="High"){
        widqqvalue<-1600
        heightqqvalue<-400
        resppi<-300
      }else{
        widqqvalue<-1600
        heightqqvalue<-400
        resppi<-100
      }
      #########################################
      try({
        if(PlotFormat=="png"){
          png(paste(dir,"/","Trait",NUM,"_resF2.png",sep=""), width=widqqvalue, height=heightqqvalue, units= unit_value, pointsize = pointsizeqqvalue,res=resppi)
        }else if(PlotFormat=="tiff"){
          tiff(paste(dir,"/","Trait",NUM,"_resF2.tiff",sep=""), width=widqqvalue, height=heightqqvalue, units= unit_value, pointsize = pointsizeqqvalue,res=resppi,compression = c("zip"))
        }else if(PlotFormat=="jpeg"){
          jpeg(paste(dir,"/","Trait",NUM,"_resF2.jpeg",sep=""), width=widqqvalue, height=heightqqvalue, units= unit_value, pointsize = pointsizeqqvalue,res=resppi)
        }else if(PlotFormat=="pdf"){
          pdf(paste(dir,"/","Trait",NUM,"_resF2.pdf",sep=""), width=16,height=6)
        }

        par(mar=c(2*margin_space+2,2*margin_space+2,0.5*margin_space,2*margin_space+2)+margin_space+2,mgp=c(3*axis_space,axis_space,2))
        max_ayis1<-ceiling(max(P_Q,P_QE)*2)

        id<-which(c(1,2,3,4,6,8,9,10)%in%max_ayis1)
        if(length(id)>0){
          Div<-c(0.5,1,1,2,2,2,3,2)[id]
        }else{
          idx<-which(c(max_ayis1%%3,max_ayis1%%4,max_ayis1%%5)==0)
          if(length(idx)>0){
            section<-c(3,4,5)[idx][1]
            Div<-max_ayis1/section
          }else{
            max_ayis1<-max_ayis1+1
            idx<-which(c(max_ayis1%%3,max_ayis1%%4,max_ayis1%%5)==0)
            if(length(idx)>0){
              section<-c(3,4,5)[which(c(max_ayis1%%3,max_ayis1%%4,max_ayis1%%5)==0)][1]
              Div<-max_ayis1/section
            }else{
              max_ayis1<-max_ayis1+1
              section<-c(3,4,5)[which(c(max_ayis1%%3,max_ayis1%%4,max_ayis1%%5)==0)][1]
              Div<-max_ayis1/section
            }
          }
        }
        plot(big_sumPOS,P_Q,type="l",col="gray60",xlim = c(0,max(big_sumPOS)),ylim = c(0,max_ayis1),
             xlab="",ylab ="",xaxt="n",yaxt="n",
             yaxs="i",xaxs="i",
             axes = F,bty="l",lwd=3,cex.axis=0.8,cex.main=0.8,cex.lab=0.8,font.lab=1.3)#,col.axis="black"
        par(new=T)
        plot(big_sumPOS,P_QE,type="l",col="lightblue",xlim = c(0,max(big_sumPOS)),ylim = c(0,max_ayis1),
             xlab="",ylab ="",xaxt="n",yaxt="n",
             yaxs="i",xaxs="i",
             axes = F,bty="l",lwd=3,cex.axis=0.8,cex.main=0.8,cex.lab=0.8,font.lab=1.3)#,col.axis="black"
        lines(c(sum_pos[-1],u),rep(max_ayis1,times=length(chr_ID)),type="h",col="gray50",lty=2,lwd=3)
        axis(4,at=seq(0,max_ayis1,Div),
             cex.axis=1.8,col.axis="black",las=0,col.ticks = "black",padj=2.7,tck=-0.05,lwd=3)
        mtext(expression(-'log'[10]*'(P-value)'),side=4,line=2.5,cex.lab=1,las=0,cex=2,padj=2.5)# cex: Coordinate label size,at=2,adj=0.5

        if( (length(Result_P_Q)>0) | (length(Result_P_QE)>0)){

          max_ayis2<-ceiling(max(Result_P_Q,Result_P_QE))

          id<-which(c(1,2,3,4,6,8,9,10)%in%max_ayis2)
          if(length(id)>0){
            Div<-c(0.5,1,1,2,2,2,3,2)[id]
          }else{
            idx<-which(c(max_ayis2%%3,max_ayis2%%4,max_ayis2%%5)==0)
            if(length(idx)>0){
              section<-c(3,4,5)[idx][1]
              Div<-max_ayis2/section
            }else{
              max_ayis2<-max_ayis2+1
              idx<-which(c(max_ayis2%%3,max_ayis2%%4,max_ayis2%%5)==0)
              if(length(idx)>0){
                section<-c(3,4,5)[which(c(max_ayis2%%3,max_ayis2%%4,max_ayis2%%5)==0)][1]
                Div<-max_ayis2/section
              }else{
                max_ayis2<-max_ayis2+1
                section<-c(3,4,5)[which(c(max_ayis2%%3,max_ayis2%%4,max_ayis2%%5)==0)][1]
                Div<-max_ayis2/section
              }
            }
          }

          if(length(Result_P_Q)>0){
            par(new=T)
            plot(ResultFinal_Q$pos_all,Result_P_Q,type="h",col="red",xlim = c(0,max(big_sumPOS)),ylim = c(0,max_ayis2),
                 xlab="",ylab ="",xaxt="n",yaxt="n",
                 yaxs="i",xaxs="i",
                 axes = F,bty="l",lwd=6,cex.axis=0.8,cex.main=1,cex.lab=1.3,font.lab=1)
          }
          if(length(Result_P_QE)>0){
            par(new=T)
            plot(ResultFinal_QE$pos_all,Result_P_QE,type="h",col="magenta",xlim = c(0,max(big_sumPOS)),ylim = c(0,max_ayis2),
                 xlab="",ylab ="",xaxt="n",yaxt="n",
                 yaxs="i",xaxs="i",
                 axes = F,bty="l",lwd=6,cex.axis=0.8,cex.main=1,cex.lab=1.3,font.lab=1)
          }

          axis(2,at=seq(0,max_ayis2,Div),
               cex.axis=1.8,col.axis="black",las=0,col.ticks = "black",padj=-2.5,tck=-0.05,lwd=3)
        }

        abline(h=CriLOD,col="black",lwd=3,lty=2)
        mtext("LOD score",side=2,cex.lab=1,las=0,cex=2,line=6.4,padj=0.5)#adj=0.5,
        mtext("Chromosomes",side=1,cex.lab=1,las=0,cex=2,line=4.5,padj=0.5)#adj=0.5,
        axis(1,at=halfsum_pos,labels=original_Chr,
             cex.axis=2,col.axis="black",col.ticks="black",# col="gray60",
             las=0,padj=0.7,tck=-0.05,lwd=3.5,lwd.ticks=3.5,line=0.75)

        dev.off()
      })
    }else{
      ###########################
      moreCHR<-as.numeric(moreMAP[,1])
      morePOS<-as.numeric(moreMAP[,2])

      original_Chr=unique(mapRaw1$chr)
      max_pos<-NULL
      chr_ID<-NULL;s<-0
      for(i in 1:max(moreCHR)){
        s<-s+length(which(moreCHR==i))
        chr_ID<-c(chr_ID,s)
        max_pos<-c(max_pos,max(morePOS[which(moreCHR==i)]))
      };rm(s)
      halfmax_pos<-max_pos*0.5
      u<-sum(max_pos)
      max_pos<-c(0,max_pos[-max(moreCHR)])
      sum_pos<-NULL;big_sumPOS<-NULL;s<-0
      for(i in 1:max(moreCHR)){
        s<-s+max_pos[i]
        sum_pos<-c(sum_pos,s)
        big_sumPOS<-c(big_sumPOS,morePOS[which(moreCHR==i)]+sum_pos[i])
      }
      halfsum_pos<-sum_pos+halfmax_pos

      ID_Q<-which(OUTresult$LOD>=CriLOD)
      ResultFinal_Q<-OUTresult[ID_Q,]
      Result_P_Q<-ResultFinal_Q$LOD

      if(dim(ResultFinal_Q)[1]>0){
        t<-NULL
        for(i in 1:dim(ResultFinal_Q)[1]){t<-c(t,sum_pos[ResultFinal_Q$Chr[i]])}
        ResultFinal_Q$pos_all<-ResultFinal_Q[,3]+t
      }
      #########################################
      margin_space<-1.5
      axis_space<-1
      unit_value<-"mm"
      pointsizeqqvalue<-30
      if(Resolution=="High"){
        widqqvalue<-1600
        heightqqvalue<-400
        resppi<-300
      }else{
        widqqvalue<-1600
        heightqqvalue<-400
        resppi<-100
      }
      #########################################
      try({
        if(PlotFormat=="png"){
          png(paste(dir,"/","Trait",NUM,"_resF2.png",sep=""), width=widqqvalue, height=heightqqvalue, units= unit_value, pointsize = pointsizeqqvalue,res=resppi)
        }else if(PlotFormat=="tiff"){
          tiff(paste(dir,"/","Trait",NUM,"_resF2.tiff",sep=""), width=widqqvalue, height=heightqqvalue, units= unit_value, pointsize = pointsizeqqvalue,res=resppi,compression = c("zip"))
        }else if(PlotFormat=="jpeg"){
          jpeg(paste(dir,"/","Trait",NUM,"_resF2.jpeg",sep=""), width=widqqvalue, height=heightqqvalue, units= unit_value, pointsize = pointsizeqqvalue,res=resppi)
        }else if(PlotFormat=="pdf"){
          pdf(paste(dir,"/","Trait",NUM,"_resF2.pdf",sep=""), width=16,height=6)
        }

        par(mar=c(2*margin_space+2,2*margin_space+2,0.5*margin_space,2*margin_space+2)+margin_space+2,mgp=c(3*axis_space,axis_space,2))
        max_ayis1<-ceiling(max(P_Q)*2)

        id<-which(c(1,2,3,4,6,8,9,10)%in%max_ayis1)
        if(length(id)>0){
          Div<-c(0.5,1,1,2,2,2,3,2)[id]
        }else{
          idx<-which(c(max_ayis1%%3,max_ayis1%%4,max_ayis1%%5)==0)
          if(length(idx)>0){
            section<-c(3,4,5)[idx][1]
            Div<-max_ayis1/section
          }else{
            max_ayis1<-max_ayis1+1
            idx<-which(c(max_ayis1%%3,max_ayis1%%4,max_ayis1%%5)==0)
            if(length(idx)>0){
              section<-c(3,4,5)[which(c(max_ayis1%%3,max_ayis1%%4,max_ayis1%%5)==0)][1]
              Div<-max_ayis1/section
            }else{
              max_ayis1<-max_ayis1+1
              section<-c(3,4,5)[which(c(max_ayis1%%3,max_ayis1%%4,max_ayis1%%5)==0)][1]
              Div<-max_ayis1/section
            }

          }
        }
        plot(big_sumPOS,P_Q,type="l",col="gray60",xlim = c(0,max(big_sumPOS)),ylim = c(0,max_ayis1),
             xlab="",ylab ="",xaxt="n",yaxt="n",
             yaxs="i",xaxs="i",
             axes = F,bty="l",lwd=3,cex.axis=0.8,cex.main=0.8,cex.lab=0.8,font.lab=1.3)#,col.axis="black"
        lines(c(sum_pos[-1],u),rep(max_ayis1,times=length(chr_ID)),type="h",col="gray50",lty=2,lwd=3)
        axis(4,at=seq(0,max_ayis1,Div),
             cex.axis=1.8,col.axis="black",las=0,col.ticks = "black",padj=2.7,tck=-0.05,lwd=3)
        mtext(expression(-'log'[10]*'(P-value)'),side=4,line=2.5,cex.lab=1,las=0,cex=2,padj=2.5)

        if(length(Result_P_Q)>0){

          max_ayis2<-ceiling(max(Result_P_Q))
          id<-which(c(1,2,3,4,6,8,9,10)%in%max_ayis2)
          if(length(id)>0){
            Div<-c(0.5,1,1,2,2,2,3,2)[id]
          }else{
            idx<-which(c(max_ayis2%%3,max_ayis2%%4,max_ayis2%%5)==0)
            if(length(idx)>0){
              section<-c(3,4,5)[idx][1]
              Div<-max_ayis2/section
            }else{
              max_ayis2<-max_ayis2+1
              idx<-which(c(max_ayis2%%3,max_ayis2%%4,max_ayis2%%5)==0)
              if(length(idx)>0){
                section<-c(3,4,5)[which(c(max_ayis2%%3,max_ayis2%%4,max_ayis2%%5)==0)][1]
                Div<-max_ayis2/section
              }else{
                max_ayis2<-max_ayis2+1
                section<-c(3,4,5)[which(c(max_ayis2%%3,max_ayis2%%4,max_ayis2%%5)==0)][1]
                Div<-max_ayis2/section
              }
            }
          }
          par(new=T)
          plot(ResultFinal_Q$pos_all,Result_P_Q,type="h",col="red",xlim = c(0,max(big_sumPOS)),ylim = c(0,max_ayis2),
               xlab="",ylab ="",xaxt="n",yaxt="n",
               yaxs="i",xaxs="i",
               axes = F,bty="l",lwd=6,cex.axis=0.8,cex.main=1,cex.lab=1.3,font.lab=1)
          axis(2,at=seq(0,max_ayis2,Div),
               cex.axis=1.8,col.axis="black",las=0,col.ticks = "black",padj=-2.5,tck=-0.05,lwd=3)
        }
        abline(h=CriLOD,col="black",lwd=3,lty=2)
        mtext("LOD score",side=2,cex.lab=1,las=0,cex=2,line=6.4,padj=0.5)#adj=0.5,
        mtext("Chromosomes",side=1,cex.lab=1,las=0,cex=2,line=4.5,padj=0.5)#adj=0.5,
        axis(1,at=halfsum_pos,labels=original_Chr,
             cex.axis=2,col.axis="black",col.ticks="black",# col="gray60",
             las=0,padj=0.7,tck=-0.05,lwd=3.5,lwd.ticks=3.5,line=0.75)

        dev.off()
      })
    }

  }
  if(method=="GCIM-QEI"){
    readraw<-Readdata(file,fileFormat,method,filecov,MCIMmap,MultiEnv)
    DoResult<-Dodata(fileFormat,Population,method,Model,readraw,MultiEnv)
    print("Running in progress, please be patient...")
    if(MultiEnv==TRUE){
      pheRaw<-DoResult$pheRaw;genRaw<-DoResult$genRaw;mapRaw1<-DoResult$mapRaw1
      yygg<-DoResult$yygg1;# cov_en<-DoResult$cov_en
      EnvNum<-DoResult$EnvNum

      if(Population=="F2"){

        ZhouMatrices<-ZhouF(pheRaw,genRaw,mapRaw1,WalkSpeed,CriLOD,dir)
        small_map<-ZhouMatrices$mapRaw
        big_map<-ZhouMatrices$genoname

        if(length(Trait)<1){
          Trait<-seq(1,length(EnvNum))
        }

        for(NUM in Trait){
          TRY1<-try({
            OutputZhou<-ZhouMethod(Model=Model,pheRaw=pheRaw,genRaw=genRaw,mapRaw=small_map,
                                   CriLOD=CriLOD,NUM=NUM,EnvNum=EnvNum,yygg=yygg,genoname=big_map,
                                   Ax0=ZhouMatrices$Ax0,Hx0=ZhouMatrices$Hx0,Bx0=ZhouMatrices$Bx0,Ax=ZhouMatrices$Ax,Hx=ZhouMatrices$Hx,Bx=ZhouMatrices$Bx,
                                   dir=dir,CriDis=CriDis,CLO=CLO)

          },silent=FALSE)
          if ('try-error' %in% class(TRY1)|| !('try-error' %in% class(TRY1))){
            TRY2<-try({
              # setwd(dir)
              en<-EnvNum[NUM]
              sum_en<-sum(EnvNum[0:(NUM-1)])
              TraitName<-data.frame(Trait=rep(paste("Trait",NUM,sep=""),dim(OutputZhou$result)[1]))
              OUTresult<-cbind(TraitName,OutputZhou$result)
              write.table(OUTresult,paste(dir,"/","Trait",NUM,"_GCIM-QEI result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
              print(paste("Trait ",NUM," has been analyzed!",sep = ""))
              if(DrawPlot==TRUE){
                GCIM_QEI_plotF2(PlotFormat=PlotFormat,Resolution=Resolution,OUTresult=OUTresult,moreMAP=big_map[,-1,drop=F],P_Q=OutputZhou$p_Q,P_QE=OutputZhou$p_QE,
                                CriLOD=CriLOD,NUM=NUM,MultiEnv)
              }
            },silent=FALSE)
          }
        }

      }else{

        warning("Please enter F2 population parameters!")

      }
    }else{
      pheRaw<-DoResult$pheRaw;genRaw<-DoResult$genRaw;mapRaw1<-DoResult$mapRaw1
      yygg<-DoResult$yygg1;#cov_en<-DoResult$cov_en

      if(Population=="F2"){

        ZhouMatrices<-ZhouF(pheRaw,genRaw,mapRaw1,WalkSpeed,CriLOD,dir)
        small_map<-ZhouMatrices$mapRaw
        big_map<-ZhouMatrices$genoname

        if(length(Trait)<1){
          Trait<-seq(1,dim(pheRaw)[2])
        }

        for(NUM in Trait){
          TRY1<-try({
            OutputZhou<-ZhouMethod_single_env(Model=Model,pheRaw=pheRaw,genRaw=genRaw,mapRaw=small_map,
                                              CriLOD=CriLOD,NUM=NUM,yygg=yygg,genoname=big_map,
                                              Ax0=ZhouMatrices$Ax0,Hx0=ZhouMatrices$Hx0,Bx0=ZhouMatrices$Bx0,Ax=ZhouMatrices$Ax,Hx=ZhouMatrices$Hx,Bx=ZhouMatrices$Bx,
                                              dir=dir,CriDis=CriDis,CLO=CLO)

          },silent=FALSE)

          if ('try-error' %in% class(TRY1)|| !('try-error' %in% class(TRY1))){
            TRY2<-try({
              # setwd(dir)
              TraitName<-data.frame(Trait=matrix(colnames(pheRaw)[NUM],nrow=dim(OutputZhou$result)[1],ncol=1))
              OUTresult<-cbind(TraitName,OutputZhou$result)
              write.table(OUTresult,paste(dir,"/","Trait",NUM,"_GCIM-QEI result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)
              print(paste("Trait ",NUM," has been analyzed!",sep = ""))
              if(DrawPlot==TRUE){
                GCIM_QEI_plotF2(PlotFormat=PlotFormat,Resolution=Resolution,OUTresult=OUTresult,moreMAP=big_map[,-1,drop=F],P_Q=OutputZhou$p_Q,P_QE=NULL,
                                CriLOD=CriLOD,NUM=NUM,MultiEnv)
              }
            },silent=FALSE)
          }
        }

      }else{

        warning("Please enter F2 population parameters!")

      }
    }

  }else if(method=="GCIM"){
    Resolution="High"
    WEN1re<-NULL; W1re<-NULL;readraw<-NULL;DoResult<-NULL;CLO<-NULL;

    readraw<-Readdata(file,fileFormat,method,filecov,MCIMmap,MultiEnv)
    DoResult<-Dodata(fileFormat,Population,method,Model,readraw,MultiEnv)

    print("Running in progress, please be patient...")

    if(is.character(file)==FALSE){CLO<-1}

    pheRaw<-DoResult$pheRaw;genRaw<-DoResult$genRaw;mapRaw1<-DoResult$mapRaw1
    flag<-DoResult$flag;flagRIL<-DoResult$flagRIL;yygg1<-DoResult$yygg1;cov_en<-DoResult$cov_en

    if(Resolution=="Low"){
      widqqvalue<-1500
      heightqqvalue<-600
      pointsizeqqvalue<-12
      resppi<-72
      ####################
      if(Population=="F2"){
        legend_size=1.2
        mainline_size=2.5
        backline_size=1.5
        margin_space=1.5
        axis_space=1.0
        logPCoff=1.5
        lodthred=2.5
      }else{
        legend_size=1.2
        mainline_size=2.5
        backline_size=1.5
        margin_space=1.5
        axis_space=1.0
        logPCoff=1.5
        lodthred=2.5
      }
    }else if(Resolution=="High"){
      widqqvalue<-10000
      heightqqvalue<-4000
      pointsizeqqvalue<-30
      resppi<-300
      ###########################
      if(Population=="F2"){
        legend_size=0.8
        mainline_size=2.5
        backline_size=0.8
        margin_space=1.5
        axis_space=1.0
        logPCoff=1.5
        lodthred=2.5
      }else{
        legend_size=0.8
        mainline_size=2.5
        backline_size=2.5
        margin_space=1.5
        axis_space=1.0
        logPCoff=1.5
        lodthred=2.5
      }
    }
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
      ################# "pos" last ID, "pos" last value

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
        ###################################change 20200914 two y axis
        par(mar=c(2*margin_space,2*margin_space,0.5*margin_space,2*margin_space)+margin_space,mgp=c(3*axis_space,axis_space,0))
        plot(newposadd,negloP,type="l",col=color2,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,
             xlim=c(0,max(newposadd)),ylim=c(0,logPCoff*max(negloP)))
        axis(side=4,cex.axis=legend_size)
        abline(v=pos_acc,lty=2,col="gray")
        par(new=TRUE)
        plot(pospic,lodpic,type="h",col=color1,yaxt="n",xlab="Genome position (cM)",ylab="",
             cex.axis=legend_size,cex.lab=legend_size+0.4,lwd=mainline_size,xlim=c(0,max(newposadd)),
             ylim=c(0,max(lodpic)))
        axis(side=2,cex.axis=legend_size)
        mtext("LOD score",side=2,line=3*axis_space,cex=legend_size+0.4,col=color1)
        abline(h=lodthred,lty=5,col=color2)
        mtext(expression(-'log'[10]*'(P-value)'),side=4,line=3*axis_space,cex=legend_size+0.4,col = color2)

      }else{
        par(mar=c(2*margin_space,2*margin_space,0.5*margin_space,2*margin_space)+margin_space,mgp=c(3*axis_space,axis_space,0))
        plot(newposadd,negloP,type="l",col=color1,xlab="Genome position (cM)",ylab=expression('Expected -log'[10]*'(P)'),cex.axis=legend_size+0.4,cex.lab=legend_size,lwd=mainline_size,xlim=c(0,max(newposadd)),ylim=c(0,logPCoff*max(negloP)))
        abline(h=lodthred,lty=5,col=color1)
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
      # firFild <- res1d[,1:2]
      # newposaddd <- as.matrix(firFild[,2])
      # for(i in 1:chr_num)
      # {
      #   temp1d <- numeric()
      #   temp1d <- which(firFild[,1]==i)
      #   if(i>1)
      #   {
      #     newposaddd[temp1d] <- newposaddd[temp1d]+pos_acc[i-1]
      #   }
      # }
      newposaddd<-newposadda
      #############newposaddd==newposadda
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
        #####################################change 20200914 two y axis
        par(mar=c(2*margin_space,2*margin_space,0.5*margin_space,2*margin_space)+margin_space,mgp=c(3*axis_space,axis_space,0))
        plot(newposadda,negloPa,type="l",col=color3,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,
             xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPa,negloPd)))
        par(new=TRUE)
        plot(newposaddd,negloPd,type="l",col=color2,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,
             xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPa,negloPd)))
        axis(side=4,at=seq(0,logPCoff*max(negloPa,negloPd),ceiling(logPCoff*max(negloPa,negloPd)/5)),cex.axis=legend_size)
        abline(v=pos_acc,lty=2,col="gray")
        par(new=TRUE)
        plot(pospic,lodpic,type="h",col=color1,yaxt="n",xlab="Genome position (cM)",ylab="",
             cex.axis=legend_size,cex.lab=legend_size+0.4,lwd=mainline_size,xlim=c(0,max(newposadda,newposaddd)),
             ylim=c(0,max(lodpic)))
        axis(side=2,cex.axis=legend_size)
        mtext("LOD score",side=2,line=3*axis_space,cex=legend_size+0.4,col=color1)
        abline(h=lodthred,lty=5,col=color2)
        mtext(expression(-'log'[10]*'(P-value)'),side=4,line=3*axis_space,cex=legend_size+0.4,col = color2)
      }else{
        ##########################change 20200914
        par(mar=c(2*margin_space,2*margin_space,0.5*margin_space,2*margin_space)+margin_space,mgp=c(3*axis_space,axis_space,0))
        plot(newposadda,negloPa,type="l",col=color3,xaxt="n",yaxt="n",xlab="",ylab="",lwd=backline_size,
             xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPa,negloPd)))
        par(new=TRUE)
        plot(newposaddd,negloPd,type="l",col=color2,yaxt="n",xlab="Genome position (cM)",ylab="",lwd=backline_size,
             xlim=c(0,max(newposadda,newposaddd)),ylim=c(0,logPCoff*max(negloPa,negloPd)))
        axis(side=4,at=seq(0,logPCoff*max(negloPa,negloPd),ceiling(logPCoff*max(negloPa,negloPd)/5)),cex.axis=legend_size)
        abline(v=pos_acc,lty=2,col="gray")
        mtext(expression(-'log'[10]*'(P-value)'),side=4,line=3*axis_space,cex=legend_size+0.4,col = color2)

      }
    }

    if(Population=="F2"){
      WEN1re<-WenF(pheRaw,genRaw,mapRaw1,yygg1,cov_en,WalkSpeed,CriLOD,dir)
      for(NUM in Trait){
        rewen<-NULL;mxmp=NULL;galaxyy1<-NULL;res1a=NULL;res1d=NULL;chr_name=NULL
        TRY1<-try({
          outWEN<-WenS(flag,CriLOD,NUM,pheRaw,Likelihood,SetSeed,flagrqtl,WEN1re$yygg,WEN1re$mx,WEN1re$phe,WEN1re$chr_name,
                       WEN1re$v.map,WEN1re$gen.raw,WEN1re$a.gen.orig,WEN1re$d.gen.orig,WEN1re$n,WEN1re$names.insert2,WEN1re$X.ad.tran.data,WEN1re$X.ad.t4,dir)
          rewen<-outWEN$result
          mxmp<-outWEN$mxmp;galaxyy1<-outWEN$galaxyy1;res1a<-outWEN$res1a;res1d<-outWEN$res1d;chr_name<-outWEN$chr_name

        },silent=FALSE)

        if ('try-error' %in% class(TRY1)|| !('try-error' %in% class(TRY1))){
          TRY2<-try({

            write.table(rewen,paste(dir,"/",NUM,"_GCIM result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)

            if(DrawPlot==TRUE){
              if(PlotFormat=="png")
              {
                png(paste(dir,"/",NUM,"_resF2.png",sep=""), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
              }else if(PlotFormat=="tiff"){
                tiff(paste(dir,"/",NUM,"_resF2.tiff",sep=""), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
              }else if(PlotFormat=="jpeg"){
                jpeg(paste(dir,"/",NUM,"_resF2.jpeg",sep=""), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
              }else if(PlotFormat=="pdf"){
                pdf(paste(dir,"/",NUM,"_resF2.pdf",sep=""), width=16)
              }
              gcimFuncF2(mxmp,galaxyy1,res1a,res1d,chr_name,legend_size,mainline_size,backline_size,margin_space,axis_space,logPCoff,"red","gray50","green",lodthred)

              dev.off()
            }
          },silent=FALSE)
        }
      }

    }else{

      W1re<-WangF(pheRaw,genRaw,mapRaw1,yygg1,flagRIL,cov_en,Population,WalkSpeed,CriLOD)

      for(NUM in Trait){
        rew<-NULL;mxmp=NULL;galaxyy1<-NULL;res11=NULL;chr_name=NULL
        TRY1<-try({
          outW<-WangS(flag,CriLOD,NUM,pheRaw,W1re$chrRaw_name,W1re$yygg,W1re$mx,W1re$phe,W1re$chr_name,W1re$gen,W1re$mapname,CLO)
          rew<-outW$result
          mxmp<-outW$mxmp;galaxyy1<-outW$galaxyy1;res11<-outW$res11;chr_name<-outW$chr_name
        },silent=FALSE)


        if ('try-error' %in% class(TRY1)|| !('try-error' %in% class(TRY1))){
          TRY2<-try({
            write.table(rew,paste(dir,"/",NUM,"_GCIM result.csv",sep=""),sep=",",row.names=FALSE,col.names = T)

            if(DrawPlot==TRUE){
              if(PlotFormat=="png")
              {
                png(paste(dir,"/",NUM,"_res.png",sep=""), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
              }else if(PlotFormat=="tiff"){
                tiff(paste(dir,"/",NUM,"_res.tiff",sep=""), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
              }else if(PlotFormat=="jpeg"){
                jpeg(paste(dir,"/",NUM,"_res.jpeg",sep=""), width=widqqvalue, height=heightqqvalue, units= "px", pointsize = pointsizeqqvalue,res=resppi)
              }else if(PlotFormat=="pdf"){
                pdf(paste(dir,"/",NUM,"_res.pdf",sep=""), width=16)
              }

              gcimFunc(mxmp,galaxyy1,res11,chr_name,legend_size,mainline_size,backline_size,margin_space,axis_space,logPCoff,"red","gray50",lodthred)

              dev.off()
            }
          },silent=FALSE)
        }
      }
    }
  }else{
    warning("Please specify the method GCIM or GCIM-QEI!")
  }

}
