\name{QTL.gCIMapping-package}
\alias{QTL.gCIMapping-package}
\alias{QTL.gCIMapping}
\docType{package}
\title{QTL Genome-Wide Composite Interval Mapping}
\description{
Conduct multi-Quantitative trait locus (QTL) mapping under the framework of random-QTL-effect mixed linear model. First, each position on the genome is detected in order to construct a negative logarithm P-value curve against genome position. Then, all the peaks on each effect (additive or dominant) curve are viewed as potential QTL, all the effects of the potential QTL are included in a multi-QTL model, their effects are estimated by empirical Bayes in doubled haploid or by adaptive lasso in F2, and true QTL are identified by likelihood radio test.
}
\details{
\tabular{ll}{
Package: \tab QTL.gCIMapping\cr
Type: \tab Package\cr
Version: \tab 2.0\cr
Date: \tab 2018-4-26\cr
Depends: \tab shiny,MASS,dplyr,parcor,qtl,doParallel\cr
Imports: \tab methods,openxlsx,stringr\cr
License: \tab GPL version 2 or newer\cr
LazyLoad: \tab yes\cr
}
}
\author{
Zhang Ya-Wen, Wen Yang-Jun, Wang Shi-Bo, Zhang Yuan-Ming\cr
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
\dontrun{QTL.gCIMapping()}
}