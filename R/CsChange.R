CsChange=function(fit1,fit2,form1=NULL,form2=NULL,data,nb=200,signif=0.05,seed=123){

  set.seed(seed)

  if.stop=c(as.character(fit1$call)[1],as.character(fit2$call)[1]) %in% c("glm","lrm","coxph","cph")
  if(any(!if.stop)){
    stop(cat(paste("The models supported currently are:",
                   "'coxph' model in 'survival' package",
                   "'cph' model in 'rms' package",
                   "'lrm' model in 'rms' package",
                   "'glm' model in 'stats' package",sep="\n")))
  }

  #############################
  ## for coxph and cph model ##
  #############################

  if(as.character(fit1$call)[1]=="coxph" | as.character(fit1$call)[1]=="cph"){

    if(is.null(fit1$formula)){fit1.formula=fit1$sformula}
    if(is.null(fit1$sformula)){fit1.formula=fit1$formula}
    if(is.null(fit2$formula)){fit2.formula=fit2$sformula}
    if(is.null(fit2$sformula)){fit2.formula=fit2$formula}

    surv=as.character(fit1.formula)[2]
    surv=gsub("Surv\\(","",surv)
    surv=gsub("\\)","",surv)
    surv=gsub(" ","",surv)
    surv.t=sapply(strsplit(surv,split=","),"[",1)
    surv.s=sapply(strsplit(surv,split=","),"[",2)
    if(length(grep("\\=",surv.t))==1){surv.t=sapply(strsplit(surv.t,split="="),"[",2)}
    if(length(grep("\\=",surv.s))==1){surv.s=sapply(strsplit(surv.s,split="="),"[",2)}
    surv.t1=surv.t
    surv.s1=surv.s

    surv=as.character(fit2.formula)[2]
    surv=gsub("Surv\\(","",surv)
    surv=gsub("\\)","",surv)
    surv=gsub(" ","",surv)
    surv.t=sapply(strsplit(surv,split=","),"[",1)
    surv.s=sapply(strsplit(surv,split=","),"[",2)
    if(length(grep("\\=",surv.t))==1){surv.t=sapply(strsplit(surv.t,split="="),"[",2)}
    if(length(grep("\\=",surv.s))==1){surv.s=sapply(strsplit(surv.s,split="="),"[",2)}
    surv.t2=surv.t
    surv.s2=surv.s

    surv.t=c(surv.t1,surv.t2)
    surv.s=c(surv.s1,surv.s2)
    surv.t=surv.t[!duplicated(surv.t)]
    surv.s=surv.s[!duplicated(surv.s)]

    x.fit1=as.character(fit1.formula)[3]
    x.fit1=gsub(" |\\n","",x.fit1)
    x.fit1=strsplit(x.fit1,"\\+")[[1]]
    if(length(grep("\\(",x.fit1))>=1){
      x.fit1[grep("\\(",x.fit1)]=sapply(strsplit(x.fit1[grep("\\(",x.fit1)],"\\("),"[",2)
      x.fit1[grep("\\)",x.fit1)]=gsub("\\)","",x.fit1[grep("\\)",x.fit1)])
    }
    x.fit2=as.character(fit2.formula)[3]
    x.fit2=gsub(" |\\n","",x.fit2)
    x.fit2=strsplit(x.fit2,"\\+")[[1]]
    if(length(grep("\\(",x.fit2))>=1){
      x.fit2[grep("\\(",x.fit2)]=sapply(strsplit(x.fit2[grep("\\(",x.fit2)],"\\("),"[",2)
      x.fit2[grep("\\)",x.fit2)]=gsub("\\)","",x.fit2[grep("\\)",x.fit2)])
    }
    x.fit=c(x.fit1,x.fit2)
    x.fit=x.fit[!duplicated(x.fit)]

    all=c(surv.t,surv.s,x.fit)

    data=data[,all]
    nrows1=nrow(data)
    data=na.omit(data)
    nrows2=nrow(data)

    if(nrows1!=nrows2){message("Note: some cases with missing value were removed.")}

    fit1.c=coxph(formula(fit1.formula),data=data)
    fit2.c=coxph(formula(fit2.formula),data=data)

    w1=rcorrcens(formula(paste("Surv(data$",surv.t1,",","data$",surv.s1,") ~",
                               "predict(fit1.c)")))
    w2=rcorrcens(formula(paste("Surv(data$",surv.t2,",","data$",surv.s2,") ~",
                               "predict(fit2.c)")))
    c1=w1[3]/2+0.5
    se1=w1[4]/2;low1=c1-1.96*se1;up1=c1+1.96*se1
    c2=w2[3]/2+0.5
    se2=w2[4]/2;low2=c2-1.96*se2;up2=c2+1.96*se2
    c12=data.frame(c=c(c1,c2),low=c(low1,low2),up=c(up1,up2))
    z=(c12$c-0.5)/((c12$up-c12$low)/(2*qnorm(1-signif*0.5)))
    p=2*pnorm(abs(z),lower.tail=F)
    c12$p=p
    row.names(c12)=c("fit1","fit2")

    if(as.character(fit1$call)[1]=="coxph" & as.character(fit2$call)[1]=="coxph"){
      change.boot=function(data,indices){
        data.boot=data[indices,]
        fit1.boot=coxph(formula(fit1.formula),data=data.boot)
        fit2.boot=coxph(formula(fit2.formula),data=data.boot)
        w1=rcorrcens(formula(paste("Surv(data.boot$",surv.t1,",","data.boot$",surv.s1,") ~",
                                   "predict(fit1.boot)")))
        w2=rcorrcens(formula(paste("Surv(data.boot$",surv.t2,",","data.boot$",surv.s2,") ~",
                                   "predict(fit2.boot)")))
        c1=w1[3]/2+0.5
        c2=w2[3]/2+0.5
        as.numeric(c2-c1)
      }
    }

    if(as.character(fit1$call)[1]=="cph" & as.character(fit2$call)[1]=="cph"){
      change.boot=function(data,indices){
        data.boot=data[indices,]
        fit1.boot=cph(formula(fit1.formula),data=data.boot)
        fit2.boot=cph(formula(fit2.formula),data=data.boot)
        c1=fit1.boot$stats["Dxy"]/2+0.5
        c2=fit2.boot$stats["Dxy"]/2+0.5
        as.numeric(c2-c1)
      }
    }

    rstb=boot(data, change.boot, R=nb)
    ci=boot.ci(rstb,type="norm")
    low=ci$normal[2]
    up=ci$normal[3]
    rst=data.frame(change=rstb$t0,low=low,up=up)
    z=rst$change/((rst$up-rst$low)/(2*qnorm(1-signif*0.5)))
    p=2*pnorm(abs(z),lower.tail=F)
    rst$p=p
    row.names(rst)="fit2-fit1"

    return(list(rst,c12))

  }#coxph or cph

  #######################
  ## for  or lrm model ##
  #######################

  if(as.character(fit1$call)[1]=="lrm"){

    if(is.null(form1)){
      form1=as.character(fit1$call)[2]
      form2=as.character(fit2$call)[2]
    }

    if(!is.null(form1)){
      y1=as.character(form1)[2]
      y2=as.character(form2)[2]
      y=c(y1,y2)
      y=y[!duplicated(y)]

      x.fit1=as.character(form1)[3]
      x.fit1=strsplit(x.fit1,"\\+")[[1]]
      x.fit1=gsub(" |\\n","",x.fit1)
      if(length(grep("\\(",x.fit1))>=1){
        x.fit1[grep("\\(",x.fit1)]=sapply(strsplit(x.fit1[grep("\\(",x.fit1)],"\\("),"[",2)
        x.fit1[grep("\\)",x.fit1)]=gsub("\\)","",x.fit1[grep("\\)",x.fit1)])
      }

      x.fit2=as.character(form2)[3]
      x.fit2=strsplit(x.fit2,"\\+")[[1]]
      x.fit2=gsub(" |\\n","",x.fit2)
      if(length(grep("\\(",x.fit2))>=1){
        x.fit2[grep("\\(",x.fit2)]=sapply(strsplit(x.fit2[grep("\\(",x.fit2)],"\\("),"[",2)
        x.fit2[grep("\\)",x.fit2)]=gsub("\\)","",x.fit2[grep("\\)",x.fit2)])
      }

      x.fit=c(x.fit1,x.fit2)
      x.fit=x.fit[!duplicated(x.fit)]
      all=c(y,x.fit)
    }

    data=data[,all]
    nrows1=nrow(data)
    data=na.omit(data)
    nrows2=nrow(data)

    if(nrows1!=nrows2){message("Note: some cases with missing value were removed.")}

    fit1.c=lrm(form1,data=data)
    fit2.c=lrm(form2,data=data)

    w1=rcorrcens(formula(paste("data$",y1,"~","predict(fit1.c)",sep="")))
    w2=rcorrcens(formula(paste("data$",y2,"~","predict(fit2.c)",sep="")))
    c1=w1[3]/2+0.5
    se1=w1[4]/2;low1=c1-1.96*se1;up1=c1+1.96*se1
    c2=w2[3]/2+0.5
    se2=w2[4]/2;low2=c2-1.96*se2;up2=c2+1.96*se2
    c12=data.frame(c=c(c1,c2),low=c(low1,low2),up=c(up1,up2))
    z=(c12$c-0.5)/((c12$up-c12$low)/(2*qnorm(1-signif*0.5)))
    p=2*pnorm(abs(z),lower.tail=F)
    c12$p=p
    row.names(c12)=c("fit1","fit2")

    change.boot=function(data,indices){
      data.boot=data[indices,]
      fit1.boot=lrm(form1,data=data.boot)
      fit2.boot=lrm(form2,data=data.boot)
      w1=rcorrcens(formula(paste("data.boot$",y1,"~","predict(fit1.boot)")))
      w2=rcorrcens(formula(paste("data.boot$",y2,"~","predict(fit2.boot)")))
      c1=w1[3]/2+0.5
      c2=w2[3]/2+0.5
      as.numeric(c2-c1)
    }

    tryb=try(boot(data, change.boot, R=nb),TRUE)
    if(inherits(tryb,"try-error")){
      message("There is problem in the 'boot' function! Try the glm model again!")
      rst=NULL
    }else{
      rstb=boot(data, change.boot, R=nb)
      if(rstb$t0==0){message("These two models are equal!");rst=NULL}
      if(rstb$t0!=0){
        ci=boot.ci(rstb,type="norm")
        low=ci$normal[2]
        up=ci$normal[3]
        rst=data.frame(change=rstb$t0,low=low,up=up)
        z=rst$change/((rst$up-rst$low)/(2*qnorm(1-signif*0.5)))
        p=2*pnorm(abs(z),lower.tail=F)
        rst$p=p
        row.names(rst)="fit2-fit1"
      }
    }

    return(list(rst,c12))

  }#lrm

  #############################
  ## for glm(binomial) model ##
  #############################

  if(as.character(fit1$call)[1]=="glm" & as.character(fit1$family)[1]=="binomial"){

    if(is.null(form1)){
      form1=as.character(fit1$call)[2]
      form2=as.character(fit2$call)[2]
    }

    if(!is.null(form1)){
      y1=as.character(form1)[2]
      y2=as.character(form2)[2]
      y=c(y1,y2)
      y=y[!duplicated(y)]

      x.fit1=as.character(form1)[3]
      x.fit1=strsplit(x.fit1,"\\+")[[1]]
      x.fit1=gsub(" ","",x.fit1)
      if(length(grep("\\(",x.fit1))>=1){
        x.fit1[grep("\\(",x.fit1)]=sapply(strsplit(x.fit1[grep("\\(",x.fit1)],"\\("),"[",2)
        x.fit1[grep("\\)",x.fit1)]=gsub("\\)","",x.fit1[grep("\\)",x.fit1)])
      }

      x.fit2=as.character(form2)[3]
      x.fit2=strsplit(x.fit2,"\\+")[[1]]
      x.fit2=gsub(" ","",x.fit2)
      if(length(grep("\\(",x.fit2))>=1){
        x.fit2[grep("\\(",x.fit2)]=sapply(strsplit(x.fit2[grep("\\(",x.fit2)],"\\("),"[",2)
        x.fit2[grep("\\)",x.fit2)]=gsub("\\)","",x.fit2[grep("\\)",x.fit2)])
      }

      x.fit=c(x.fit1,x.fit2)
      x.fit=x.fit[!duplicated(x.fit)]
      all=c(y,x.fit)
    }

    data=data[,all]
    nrows1=nrow(data)
    data=na.omit(data)
    nrows2=nrow(data)

    if(nrows1!=nrows2){message("Note: some cases with missing value were removed.")}

    fit1.c=glm(form1,data=data,family = "binomial")
    fit2.c=glm(form2,data=data,family = "binomial")

    w1=rcorrcens(formula(paste("data$",y1,"~","predict(fit1.c)",sep="")))
    w2=rcorrcens(formula(paste("data$",y2,"~","predict(fit2.c)",sep="")))
    c1=w1[3]/2+0.5
    se1=w1[4]/2;low1=c1-1.96*se1;up1=c1+1.96*se1
    c2=w2[3]/2+0.5
    se2=w2[4]/2;low2=c2-1.96*se2;up2=c2+1.96*se2
    c12=data.frame(c=c(c1,c2),low=c(low1,low2),up=c(up1,up2))
    z=(c12$c-0.5)/((c12$up-c12$low)/(2*qnorm(1-signif*0.5)))
    p=2*pnorm(abs(z),lower.tail=F)
    c12$p=p
    row.names(c12)=c("fit1","fit2")

    change.boot=function(data,indices){
      data.boot=data[indices,]
      fit1.boot=glm(form1,data=data.boot,family = "binomial")
      fit2.boot=glm(form2,data=data.boot,family = "binomial")
      w1=rcorrcens(formula(paste("data.boot$",y1,"~","predict(fit1.boot)")))
      w2=rcorrcens(formula(paste("data.boot$",y2,"~","predict(fit2.boot)")))
      c1=w1[3]/2+0.5
      c2=w2[3]/2+0.5
      as.numeric(c2-c1)
    }

    tryb=try(boot(data, change.boot, R=nb),TRUE)
    if(inherits(tryb,"try-error")){
      message("There is problem in the 'boot' function!")
      rst=NULL
    }else{
      rstb=boot(data, change.boot, R=nb)
      if(rstb$t0==0){message("These two models are equal!");rst=NULL}
      if(rstb$t0!=0){
        ci=boot.ci(rstb,type="norm")
        low=ci$normal[2]
        up=ci$normal[3]
        rst=data.frame(change=rstb$t0,low=low,up=up)
        z=rst$change/((rst$up-rst$low)/(2*qnorm(1-signif*0.5)))
        p=2*pnorm(abs(z),lower.tail=F)
        rst$p=p
        row.names(rst)="fit2-fit1"
      }
    }

    return(list(rst,c12))

  }#glm

}
