\name{WenS}
\alias{WenS}
\title{The second step of Wen method}
\description{
An efficient multi-locus mixed model framework for the detection of small and linked QTLs in F2}
\usage{
WenS(flag,CriLOD,NUM,pheRaw,Likelihood,flagrqtl,yygg,mx,phe,chr_name,v.map,
     gen.raw,a.gen.orig,d.gen.orig,n,names.insert2,X.ad.tran.data,X.ad.t4,dir)
}
\arguments{
      \item{flag}{random or fix model.}
	\item{CriLOD}{LOD score.}
	\item{NUM}{the number of trait.}
	\item{pheRaw}{raw phenotype matrix .}
	\item{Likelihood}{likelihood function.}
	\item{flagrqtl}{do CIM or not.}
      \item{yygg}{covariate matrix.}
      \item{mx}{raw genotype matrix.}
      \item{phe}{phenotype matrix.}
	\item{chr_name}{chromosome name.}
	\item{v.map}{linkage map matrix.}
	\item{gen.raw}{raw genotype matrix.}
	\item{a.gen.orig}{additive genotype matrix.}
	\item{d.gen.orig}{dominant genotype matrix.}
      \item{n}{number of individual.}
      \item{names.insert2}{linkage map after insert.}
      \item{X.ad.tran.data}{genotype matrix after insert.}
	\item{X.ad.t4}{genotype matrix.}
      \item{dir}{file storage path.}
       }
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, and Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\examples{
F2=data(F2data)
readraw<-Readdata(file=F2data,fileFormat="GCIM",fileICIMcov=NULL)
DoResult<-Dodata(fileFormat="GCIM",Population="F2",Model="Random",readraw)
WEN1re<-WenF(DoResult$pheRaw,DoResult$genRaw,DoResult$mapRaw1,
DoResult$yygg1,DoResult$cov_en,WalkSpeed=1,CriLOD=2.5,
dir=tempdir())
WenS(DoResult$flag,CriLOD=2.5,NUM=1,DoResult$pheRaw,Likelihood="REML",
flagrqtl=FALSE,WEN1re$yygg,WEN1re$mx,WEN1re$phe,WEN1re$chr_name,
WEN1re$v.map,WEN1re$gen.raw,WEN1re$a.gen.orig,WEN1re$d.gen.orig,
WEN1re$n,WEN1re$names.insert2,WEN1re$X.ad.tran.data,WEN1re$X.ad.t4,
dir=tempdir())
}




