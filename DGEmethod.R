rm(list=setdiff(ls(), "x"))
#---------------------all packages I needed-------------------------
P.clean<-c("dplyr", "Rmisc", "reshape2","stringr")
P.plot <-c("ggplot2", "vioplot", "vioplot","ggpubr","RColorBrewer","ggmap","devEMF")
P.analy<-c("mclust","dplR")
lapply(P.clean, require, character.only = TRUE)
lapply(P.plot,  require, character.only = TRUE)
lapply(P.analy, require, character.only = TRUE)

filedir<-"C:/03_FEcology/02_Disturbances_Global/01_DGE_disturbances/Global_TRW/"

DGE<-function(trw,trw.fn0exten){
  # Line 1: chronology -index [crn.NegExpHUR.cut$HURstd] ------------------------
  rwi.NegExp<-trw
  #!!!select for real data and simulation data
  for (i in 1:length(trw)) {
    ser<-trw[,i]# 
    # ii<-try(detrend.series(ser,make.plot = F,method = "Friedman",
    #                        verbose = T,return.info = T, make.plot = T),silent=T)
    ii<-detrend.series(ser,method = "Friedman",verbose = F,return.info = T)
    rwi.NegExp[,i]<-ii[[1]]
  }
  # rwi.NegExp<-i.detrend(raw)
  #!!!select for real data and simulation data
  crn.NegExpHUR     <- chron(rwi.NegExp, prefix="HUR")
  # crn.NegExpHUR.cut <- crn.NegExpHUR[(dim(crn.NegExpHUR)[1]-99):dim(crn.NegExpHUR)[1],] 
  crn.NegExpHUR.cut <- subset(crn.NegExpHUR, samp.depth>30)
  #1:names=HURstd,samp.depth. #2:1940-2010, 71 numbers,
  # crn.NegExpHUR.cut <- crn.NegExpHUR[(dim(crn.NegExpHUR)[1]-99):dim(crn.NegExpHUR)[1],] #1:names=HURstd,samp.depth. #2:1940-2010, 71 numbers,
  # rawr<-subset(chro, samp.depth>10)
  # start<-as.numeric(row.names(rawr)[1])
  # length<-dim(rawr)[1]
  # end<-start+length-1
  # chrc<-subset(chro,rownames(chro) %in% c(start:end))
  # raww<-subset(raw,rownames(raw)   %in% c(start:end))
  # Line 2: chronology -sample size [crn.NegExpHUR.cut$samp.depth]--------------------------------------  
  # Line 3: density - density [rwi.NegExp]--------------------------------------------------------------
  # Line 4: barplot of GMM - [GMM$"Cm1g3","Cm2g3","Cm3g3","Vr1g3","Vr2g3","Vr3g3",]---------------------
  
  start<-as.numeric(row.names(crn.NegExpHUR.cut)[1])
  end  <-start+dim(crn.NegExpHUR.cut)[1]-1  # end<-start+length-1
  rwi  <-subset(rwi.NegExp,rownames(rwi.NegExp)   %in% c(start:end))
  x<-rwi[12,]
  den<-apply(rwi,1,function(x){
    # density =try(Mclust(data=x[!is.na(x)]))
    # summ=try(summary(density,parameters=TRUE))
    density =try(densityMclust(data=x[!is.na(x)]))#modelNames = "V"
    summ=try(summary(density,parameters=TRUE))
    G.optimal =try(summ$G)
    loglik.op=try(summ$loglik)
    bic.op=try(summ$bic)   # they use bic as the standard to select the number of cluster...
    icl.op=try(summ$icl)
    
    density1 =try(densityMclust(data=x[!is.na(x)],G=1))
    summ1=try(summary(density1,parameters=TRUE))
    loglik_g1=try(summ1$loglik)    # larger is better 
    bic_g1=try(summ1$bic)          # smaller better
    icl_g1=try(summ1$icl)          # smaller better
    Cm1g1=try(summ1$mean[1])
    Vr1g1=try(summ1$variance[1])
    Pr1g1=try(summ1$pro[1])
    
    density2 =try(densityMclust(data=x[!is.na(x)],G=2))
    summ2=try(summary(density2,parameters=TRUE))
    loglik_g2=try(summ2$loglik) # larger is better 
    bic_g2=try(summ2$bic)       # smaller better
    icl_g2=try(summ2$icl)       # smaller better
    Cm1g2=try(summ2$mean[1])
    Cm2g2=try(summ2$mean[2])
    Vr1g2=try(summ2$variance[1])
    Vr2g2=try(summ2$variance[2])
    Pr1g2=try(summ2$pro[1])
    Pr2g2=try(summ2$pro[2])
    
    density3 =try(densityMclust(data=x[!is.na(x)],G=3))
    summ3=try(summary(density3,parameters=TRUE))
    loglik_g3=try(summ3$loglik) # larger is better 
    bic_g3=try(summ3$bic)       # smaller better
    icl_g3=try(summ3$icl)       # smaller better
    Cm1g3=try(summ3$mean[1])
    Cm2g3=try(summ3$mean[2])
    Cm3g3=try(summ3$mean[3])
    Vr1g3=try(summ3$variance[1])
    Vr2g3=try(summ3$variance[2])
    Vr3g3=try(summ3$variance[3])
    Pr1g3=try(summ3$pro[1])
    Pr2g3=try(summ3$pro[2])
    Pr3g3=try(summ3$pro[3])
    
    bic.min=max(bic_g1,bic_g2,bic_g3)
    # bic.max=max(icl_g1,icl_g2,icl_g3) # maximum of negative number and min of positive number
    # bic.max=max(icl_g1,icl_g2,icl_g3) # maximum of negative number and min of positive number
    # if (loglik_g1==loglik.max){ # for using loglik as justic rule, loglik not default
    if (bic_g1==bic.min){
      G.123 =try(summ1$G)
      Cm1=try(summ1$mean[1])
      Cm2=NA
      Cm3=NA
      Vr1=try(summ1$variance[1])
      Vr2=NA
      Vr3=NA
      Pr1=try(summ1$pro[1])
      Pr2=NA
      Pr3=NA
      loglik=try(-1*summ1$loglik) # larger is better 
      bic=try(-1*summ1$bic)       # smaller better
      icl=try(-1*summ1$icl)       # smaller better
    } else if (bic_g2==bic.min){
      #} else if (icl_g2==bic.max){
      G.123 =try(summ2$G)
      Cm1=try(summ2$mean[1])
      Cm2=try(summ2$mean[2])
      Cm3=NA
      Vr1=try(summ2$variance[1])
      Vr2=try(summ2$variance[2])
      Vr3=NA
      Pr1=try(summ2$pro[1])
      Pr2=try(summ2$pro[2])
      Pr3=NA
      loglik=try(-1*summ2$loglik) # larger is better 
      bic=try(-1*summ2$bic)       # smaller better
      icl=try(-1*summ2$icl)       # smaller better
    } else if (bic_g3==bic.min){  
      #} else if (icl_1.2==bic.max){
      G.123 =try(summ3$G)
      Cm1=try(summ3$mean[1])
      Cm2=try(summ3$mean[2])
      Cm3=try(summ3$mean[3])
      Vr1=try(summ3$variance[1])
      Vr2=try(summ3$variance[2])
      Vr3=try(summ3$variance[3])
      Pr1=try(summ3$pro[1])
      Pr2=try(summ3$pro[2])
      Pr3=try(summ3$pro[3])
      loglik=try(-1*summ3$loglik)  # larger is better 
      bic=try(-1*summ3$bic)        # smaller better
      icl=try(-1*summ3$icl)        # smaller better
    } else{    
      G.123 =NA
      Cm1=NA
      Cm2=NA
      Cm3=NA
      Vr1=NA
      Vr2=NA
      Vr3=NA
      Pr1=NA
      Pr2=NA
      Pr3=NA
      loglik=NA   # larger is better 
      bic=NA      # smaller better
      icl=NA      # smaller better
    }
    
    if (!is.na(Pr1)==T && !is.na(Pr2)==F && !is.na(Pr3)==F){
      C.R=NA     # 1 cluster to 1 class
      C.H=Cm1
      C.D=NA
      P.R=NA
      P.H=Pr1
      P.D=NA
      V.R=NA
      V.H=Vr1
      V.D=NA
    }else if (!is.na(Pr1)==T && !is.na(Pr2)==T && !is.na(Pr3)==F && Pr2<0.3){
      G.123 =2  
      C.R=Cm2   # 2 clusters to 2 classes, release class less than 15%.
      C.H=Cm1
      C.D=NA
      P.R=Pr2
      P.H=Pr1
      P.D=NA
      V.R=Vr2
      V.H=Vr1
      V.D=NA
    }else if (!is.na(Pr1)==T && !is.na(Pr2)==T && !is.na(Pr3)==F && Pr2>0.3){
      G.123 =2  
      C.R=NA    # 2 clusters to 2 classes
      C.H=Cm2
      C.D=Cm1
      P.R=NA
      P.H=Pr2
      P.D=Pr1
      V.R=NA
      V.H=Vr2
      V.D=Vr1
    }else if (!is.na(Pr1)==T && !is.na(Pr2)==T && !is.na(Pr3)==T && Pr3+Pr2<0.6 && bic_g1>bic_g2){
      G.123 =1    # 3 clusters to 1 class
      C.R=NA   
      C.H=Cm1g1
      C.D=NA
      P.R=NA
      P.H=Pr1g1
      P.D=NA
      V.R=NA
      V.H=Vr1g1
      V.D=NA 
    }else if (!is.na(Pr1)==T && !is.na(Pr2)==T && !is.na(Pr3)==T && Pr3+Pr2<0.6 && bic_g2>bic_g1 && Pr2g2<0.3){
      G.123 =2
      C.R=Cm2g2   # 3 clusters to 2 classes, with release class
      C.H=Cm1g2
      C.D=NA
      P.R=Pr2g2
      P.H=Pr1g2
      P.D=NA
      V.R=Vr2g2
      V.H=Vr1g2
      V.D=NA 
    }else if (!is.na(Pr1)==T && !is.na(Pr2)==T && !is.na(Pr3)==T && Pr3+Pr2<0.6 && bic_g2>bic_g1 && Pr2g2>0.3){
      G.123 =2
      C.R=NA      # 3 clusters to 2 classes, without release class
      C.H=Cm2g2
      C.D=Cm1g2
      P.R=NA
      P.H=Pr2g2
      P.D=Pr1g2
      V.R=NA
      V.H=Vr2g2
      V.D=Vr1g2 
    }else if (!is.na(Pr1)==T && !is.na(Pr2)==T && !is.na(Pr3)==T && Pr1<0.333 && Pr3<0.3){
      G.123 =3
      C.R=Cm3   # optimal situation...
      C.H=Cm2
      C.D=Cm1
      P.R=Pr3
      P.H=Pr2
      P.D=Pr1
      V.R=Vr3
      V.H=Vr2
      V.D=Vr1
    }else if (!is.na(Pr1)==T && !is.na(Pr2)==T && !is.na(Pr3)==T && Pr1+Pr2<0.6&& bic_g1>bic_g2){
      G.123 =1  # 3 clusters to 1 classes
      C.R=NA
      C.H=Cm1g1
      C.D=NA
      P.R=NA
      P.H=Pr1g1
      P.D=NA
      V.R=NA
      V.H=Vr1g1
      V.D=NA 
    }else if (!is.na(Pr1)==T && !is.na(Pr2)==T && !is.na(Pr3)==T && Pr1+Pr2<0.6 && bic_g2>bic_g1 && Pr2g2<0.3){
      G.123 =2   # 3 clusters to 2 classes, without release class
      C.R=Cm2g2
      C.H=Cm1g2
      C.D=NA
      P.R=Pr2g2
      P.H=Pr1g2
      P.D=NA
      V.R=Vr2g2
      V.H=Vr1g2
      V.D=NA 
    }else if (!is.na(Pr1)==T && !is.na(Pr2)==T && !is.na(Pr3)==T && Pr1+Pr2<0.6&& bic_g2>bic_g1 && Pr2g2>0.3){
      G.123 =2    # 3 clusters to 2 classes, without release class
      C.R=NA
      C.H=Cm2g2
      C.D=Cm1g2
      P.R=NA
      P.H=Pr2g2
      P.D=Pr1g2
      V.R=NA
      V.H=Vr2g2
      V.D=Vr1g2  
    }else{ 
      C.R=NA  
      C.H=NA
      C.D=NA
      P.R=NA
      P.H=NA
      P.D=NA
      V.R=NA
      V.H=NA
      V.D=NA
    }
    
    c(G.optimal,G.123,C.D,C.H,C.R,P.D,P.H,P.R,V.D,V.H,V.R,
      Cm1,Cm2,Cm3,Pr1,Pr2,Pr3,Vr1,Vr2,Vr3,loglik.op,bic.op,icl.op,
      Cm1g1,Pr1g1,Vr1g1,Cm1g2,Cm2g2,Pr1g2,Pr2g2,Vr1g2,Vr2g2,
      Cm1g3,Cm2g3,Cm3g3,Pr1g3,Pr2g3,Pr3g3,Vr1g3,Vr2g3,Vr3g3,
      loglik_g1,loglik_g2,loglik_g3,bic_g1,bic_g2,bic_g3,
      icl_g1,icl_g2,icl_g3)
  })
  densit<-as.data.frame(t(den))
  GMMEM.RHD<-cbind(crn.NegExpHUR.cut,densit)
  names(GMMEM.RHD) <-c("chro","Samdepth","G.optimal","G.123",
                       "C.D","C.H","C.R","P.D","P.H","P.R","V.D","V.H","V.R",
                       "Cm1","Cm2","Cm3","Pr1","Pr2","Pr3","Vr1","Vr2","Vr3","loglik.op","bic.op","icl.op",
                       "Cm1g1","Pr1g1","Vr1g1","Cm1g2","Cm2g2","Pr1g2","Pr2g2","Vr1g2","Vr2g2",
                       "Cm1g3","Cm2g3","Cm3g3","Pr1g3","Pr2g3","Pr3g3","Vr1g3","Vr2g3","Vr3g3",
                       "loglik.1","loglik.2","loglik.3","bic.1","bic.2","bic.3","icl.1","icl.2","icl.3")
  GMMEM.filename <- paste(filedir,trw.fn0exten,"_DGE.RHD.V","_cm.pr.vr123g",".csv", sep="")
  write.csv(GMMEM.RHD,file = GMMEM.filename)
  RHD.filename <- paste(filedir,trw.fn0exten,"_DGE.RHD.V","_cm.pr.vr",".csv", sep="")
  write.csv(GMMEM.RHD[,c(1,2,4,8,9,10)],file = RHD.filename)
} 



trwf <- list.files(path = filedir,pattern = ".rwl",full.names=TRUE)
for (i in 113){
  trw.filename<-basename(trwf[i])
  trw.fn0exten<-sub('\\.rwl$', '', trw.filename) 
  trw<-read.rwl(trwf[i])
  DGE(trw,trw.fn0exten)
}