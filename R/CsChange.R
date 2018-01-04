CsChange=function(fit1,fit2,data,nb=100,signif=0.05){

  if.stop=c(as.character(fit1$call)[1],as.character(fit2$call)[1]) %in% c("coxph","cph")
  if(any(!if.stop)){
    stop(paste("Only 'coxph' model in 'survival' package and 'cph' model in 'rms'",
              "package are supported currently!"))
  }

  surv=sapply(strsplit(as.character(fit1$call)[2],split="~"),"[",1)
  surv=gsub("Surv","",surv)
  surv=gsub("\\(","",surv)
  surv=gsub("\\)","",surv)
  surv=gsub(" ","",surv)
  surv.t=sapply(strsplit(surv,split=","),"[",1)
  surv.s=sapply(strsplit(surv,split=","),"[",2)

  x.fit1=sapply(strsplit(as.character(fit1$call)[2],split="~"),"[",2)
  x.fit1=gsub(" ","",x.fit1)
  x.fit2=sapply(strsplit(as.character(fit2$call)[2],split="~"),"[",2)
  x.fit2=gsub(" ","",x.fit2)
  x.fit2=strsplit(x.fit2,"\\+")[[1]]
  x.fit=c(x.fit1,x.fit2)
  x.fit=x.fit[!duplicated(x.fit)]

  all=c(surv.t,surv.s,x.fit)

  data=data[,all]
  nrows1=nrow(data)
  data=na.omit(data)
  nrows2=nrow(data)

  if(nrows1!=nrows2){message("Note: some cases with missing value were removed.")}

  if(as.character(fit1$call)[1]=="coxph" & as.character(fit2$call)[1]=="coxph"){
    change.boot=function(data,indices){
      data.boot=data[indices,]
      fit1.boot=coxph(formula(as.character(fit1$call)[2]),data=data.boot)
      fit2.boot=coxph(formula(as.character(fit2$call)[2]),data=data.boot)
      w1=rcorrcens(formula(paste("Surv(data.boot$",surv.t,",","data.boot$",surv.s,") ~",
                                 "predict(fit1.boot)")))
      w2=rcorrcens(formula(paste("Surv(data.boot$",surv.t,",","data.boot$",surv.s,") ~",
                                 "predict(fit2.boot)")))
      c1=w1[3]/2+0.5
      c2=w2[3]/2+0.5
      as.numeric(c2-c1)
    }
  }

  if(as.character(fit1$call)[1]=="cph" & as.character(fit2$call)[1]=="cph"){
    change.boot=function(data,indices){
      data.boot=data[indices,]
      fit1.boot=cph(formula(as.character(fit1$call)[2]),data=data.boot)
      fit2.boot=cph(formula(as.character(fit2$call)[2]),data=data.boot)
      c1=fit1.boot$stats["Dxy"]/2+0.5
      c2=fit2.boot$stats["Dxy"]/2+0.5
      as.numeric(c2-c1)
    }
  }

  rstb <- boot(data, change.boot, R=nb)
  low=rstb$t0-qnorm(1-signif*0.5)*sd(rstb$t)
  up=rstb$t0+qnorm(1-signif*0.5)*sd(rstb$t)
  rst=data.frame(change=rstb$t0,low=low,up=up)
  z=rst$change/((rst$up-rst$low)/(2*qnorm(1-signif*0.5)))
  p=2*pnorm(abs(z),lower.tail=F)
  rst$p=p

  return(rst)
}


