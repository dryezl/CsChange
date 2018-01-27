\name{CsChange}
\alias{CsChange}
\title{Testing for Change in C-Statistic}
\usage{
CsChange(fit1,fit2,data,nb,signif)
}
\description{
Calculate the confidence interval and p value for change in C-statistic. The adjusted C-statistic is calculated by using formula as "Somers' Dxy rank correlation"/2+0.5. The confidence interval was calculated by using the bootstrap method. The p value was calculated by using the Z testing method. Please refer to the article of Peter Ganz et al. (2016) <doi:10.1001/jama.2016.5951>.
}
\arguments{
  \item{fit1}{an object from 'cph' or 'coxph' model}
  \item{fit2}{another object from 'cph' or 'coxph' model}
  \item{data}{a data frame used in the fit1 or fit2}
  \item{nb}{the number of bootstrap replicate, with a default of 100}
  \item{signif}{the significant level of confidence interval, with a default of 0.05 and a two-sided test}
}
\value{
  \item{change}{in list 1, the change of C-statistic from fit1 to fit2}
  \item{low,up}{in list 1, the confidence interval of change}
  \item{p}{in list 1, the p value of testing for change}
  \item{c}{in list 2, the C-statistic of fit1 and fit2}
  \item{low,up}{in list 2, the confidence interval of change}
}
\author{
Zhicheng Du<dgdzc@hotmail.com>, Yuantao Hao<haoyt@mail.sysu.edu.cn>
}
\references{
Peter Ganz, Bettina Heidecker, Kristian Hveem, Christian Jonasson, Shintaro Kato, Mark R. Segal, David G. Sterling, Stephen A. Williams. Development and Validation of a Protein-Based Risk Score for Cardiovascular Outcomes Among Patients With Stable Coronary Heart Disease. JAMA. 2016; 315(23):2532-2541. doi:10.1001/jama.2016.5951
}
\note{
Please feel free to contact us, if you have any advice and find any bug!

Update description:

Version 0.1.1: Fix the error of "variables not found" for the 'coxph' model.

Version 0.1.2: Fix the error of data with missing value.

Version 0.1.3: a) Add the c statistic of fit1 and fit2 to the output. b) Fix the error of model with multiple variables. c) Fix the error of independent variable with functions, i.e., rcc(), ploy().

more functions will be included in 'CsChange' package!
}
\examples{
require("rms")
set.seed(123)
n=50
age=50+12*rnorm(n)
sex=factor(sample(c('Male','Female'), n,rep=TRUE, prob=c(.6, .4)))
cens=15*runif(n)
h=.02*exp(.04*(age-50)+.8*(sex=='Female'))
dt=-log(runif(n))/h
e=ifelse(dt <= cens,1,0)
dt=pmin(dt, cens)
units(dt)="Year"
data=data.frame(dt,e,age,sex)
dd=datadist(age, sex)
options(datadist='dd')

#for 'cph' model
fit1=cph(Surv(dt,e)~age,data)
fit2=cph(Surv(dt,e)~age+sex,data)
CsChange(fit1,fit2,data,nb=20)

#for 'coxph' model
fit1=coxph(Surv(dt,e)~age,data)
fit2=coxph(Surv(dt,e)~age+sex,data)
CsChange(fit1,fit2,data,nb=20)

#for 'coxph' model
fit1=coxph(Surv(dt,e)~age+sex,data)
fit2=coxph(Surv(dt,e)~rcs(age)+sex,data)
CsChange(fit1,fit2,data,nb=20)#ignore the warnings messages
}
