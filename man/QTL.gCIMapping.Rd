\name{QTL.gCIMapping}
\alias{QTL.gCIMapping}
\alias{QTL.gCIMapping}
\docType{package}
\title{QTL Genome-Wide Composite Interval Mapping}
\description{
Conduct multiple quantitative trait loci (QTL) mapping under the framework of random-QTL-effect mixed linear model. First, each position on the genome is detected in order to construct a negative logarithm P-value curve against genome position. Then, all the peaks on each effect (additive or dominant) curve are viewed as potential QTL, all the effects of the potential QTL are included in a multi-QTL model, their effects are estimated by empirical Bayes in doubled haploid or by adaptive lasso in F2, and true QTL are identified by likelihood radio test.}
\details{
\tabular{ll}{
Package: \tab QTL.gCIMapping\cr
Type: \tab Package\cr
Version: \tab 3.1\cr
Date: \tab 2018-10-16\cr
Depends: \tab MASS,dplyr,parcor,qtl,doParallel\cr
Imports: \tab methods,openxlsx,stringr,Rcpp\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
}
\usage{
 QTL.gCIMapping(file,fileFormat,fileICIMcov,Population,Model,WalkSpeed,CriLOD,
 Likelihood,flagrqtl,DrawPlot,PlotFormat,Resolution,Trait,dir) 
 }
\arguments{
	\item{file}{File path and name in your computer.}
	\item{fileFormat}{Format for input file (GCIM, ICIM, Cart).}
	\item{fileICIMcov}{File path and name in your computer.}
	\item{Population}{BC1, BC2, DH, RIL, F2.}
	\item{Model}{Random (random model) or Fixed (fixed model) for QTL effects.}
	\item{WalkSpeed}{Walk speed for Genome-wide Scanning.(WalkSpeed=1)}
       \item{CriLOD}{Critical LOD scores for significant QTL (CriLOD=2.5).}
       \item{Likelihood}{This parameter is only for F2 population, including restricted maximum likelihood (REML) and maximum likelihood (ML).}
       \item{flagrqtl}{This parameter is only for F2 population,   flagrqtl="FALSE" in the first running. If the other software detects only one QTL in a neighborhood but this software finds two linked QTLs (one with additive effect and another with dominant effect) in the region, let flagrqtl="TRUE"}
      \item{DrawPlot}{This parameter is for all the populations, including FALSE and TRUE, DrawPlot=FALSE indicates no figure output, DrawPlot=TRUE indicates the output of the figure against genome position.}
	\item{PlotFormat}{This parameter is for all the figure files, including *.jpeg, *.png, *.tiff and *.pdf.}
	\item{Resolution}{This parameter is for all the figure files, including Low and High.}
       \item{Trait}{Trait=1:3 indicates the analysis from the first trait to the third trait.}
       \item{dir}{This parameter is for the save path.}

}
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, and Zhang Yuan-Ming\cr
Maintainer: Yuanming Zhang<soyzhang@mail.hzau.edu.cn>
}
\references{
Mapping small-effect and linked quantitative trait loci for complex traits in backcross or DH populations via a multi-locus GWAS methodology.
Wang Shi-Bo,Wen Yang-Jun,Ren Wen-Long,Ni Yuan-Li,Zhang Jin,Feng Jian-Ying,Zhang Yuan-Ming*
}
\keyword{
 R
QTL.gCIMapping
package
} 
\examples{
G=data(f2data)
QTL.gCIMapping(file=f2data,fileFormat="GCIM",fileICIMcov=NULL,Population="F2",
Model="Random",WalkSpeed=1,CriLOD=2.5,Likelihood="REML",flagrqtl="FALSE",
DrawPlot="FALSE",PlotFormat="png",Resolution="Low",Trait=1,dir=tempdir()) 
}