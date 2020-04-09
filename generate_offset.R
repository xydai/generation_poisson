library(tidyverse)
library(plyr)
library(dplyr)
library(poisson)
source("try_pp_time.R")
library(random)
library(MASS)
library(survival)
library(lme4)
gen_Data_ag <- function(n, beta, r){
  age<-NULL
  treatmentA<-NULL
  gender<-NULL
  rate_run<-NULL
  censored<-NULL
  process<-NULL
  jump<-NULL
  y_c<-NULL
  y_run<-NULL
  id<-NULL
  event_ag<-NULL
  event<-NULL
  time_end<-NULL
  time_start<-NULL
  ag_covariates<-NULL
  age_ag<-NULL
  gender_ag<-NULL
  treatmentA_ag<-NULL
  status<-NULL
  os<-NULL
  
  age<-round(runif(n=n,min=2,max=12))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  rate_run<-r*exp(beta[1]+beta[2]*age+beta[3]*gender+beta[4]*treatmentA)
  censored<-gen_censor2(rate_run,4)
  
  process<-gen_timing(n=100,rate_run,censored)
  jump<-process$jump
  y_c<-process$y_c
  y1<-process$y
  y_run<-process$y_new
  os<-process$off
  status<-process$status
  
  id <- gen_id(y_c)
  event_ag <-gen_eventag(y_c)
  event <- gen_event(y_c)
  time_end <- gen_time_nh(jump=jump)
  time_start <-gen_t0(id,time_end,event)
  
  ag_covariates <-gen_covariate(age,gender,treatmentA,y_c)
  age_ag <-ag_covariates$age1
  gender_ag <-ag_covariates$gender1
  treatmentA_ag <-ag_covariates$treatment
  
  return(list(age=age, gender=gender, treatmentA=treatmentA, y1=y1,
              y_c=y_c,y_run=y_run,id=id, status=status, event=event, event_ag=event_ag,
              time_end=time_end, time_start=time_start, age_ag=age_ag,
              gender_ag=gender_ag, treatmentA_ag=treatmentA_ag,os=os))
}

gen_Data_mcar <- function(n, beta, r){
  age<-NULL
  treatmentA<-NULL
  gender<-NULL
  rate_run<-NULL
  censored<-NULL
  process<-NULL
  jump<-NULL
  y_c<-NULL
  y_run<-NULL
  id<-NULL
  event_ag<-NULL
  event<-NULL
  time_end<-NULL
  time_start<-NULL
  ag_covariates<-NULL
  age_ag<-NULL
  gender_ag<-NULL
  treatmentA_ag<-NULL
  status<-NULL
  os<-NULL
  
  age<-round(runif(n=n,min=2,max=12))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  rate_run<-r*exp(beta[1]+beta[2]*age+beta[3]*gender+beta[4]*treatmentA)
  censored<-gen_censor(rate_run,0.001)
  
  process<-gen_timing(n=100,rate_run,censored)
  jump<-process$jump
  y_c<-process$y_c
  y_run<-process$y_new
  os<-process$off
  status<-process$status
  
  id <- gen_id(y_c)
  event_ag <-gen_eventag(y_c)
  event <- gen_event(y_c)
  time_end <- gen_time_nh(jump=jump)
  time_start <-gen_t0(id,time_end,event)
  
  ag_covariates <-gen_covariate(age,gender,treatmentA,y_c)
  age_ag <-ag_covariates$age1
  gender_ag <-ag_covariates$gender1
  treatmentA_ag <-ag_covariates$treatment
  
  return(list(age=age, gender=gender, treatmentA=treatmentA,
              y_c=y_c,y_run=y_run,id=id, status=status, event=event, event_ag=event_ag,
              time_end=time_end, time_start=time_start, age_ag=age_ag,
              gender_ag=gender_ag, treatmentA_ag=treatmentA_ag,os=os))
}


gen_y <- function(x, z, m, beta){
  mu<-beta[1]+beta[2]*x+beta[3]*z+beta[4]*m
  y<-rpois(length(x),exp(mu))
  return(y)
}



gen_id<-function(y){
  n=length(y)
  max<-max(y)
  
  temp=matrix(0,n,max)
  ids=rep(0,1)
  for (i in 1:n){
    #if (y[i] != 0){
    for (j in 1:y[i]){
      temp[i,j]=i
    }
    #} else{
    #  temp[i,1]=i
    #}
    ids=c(ids,temp[i,])
  }
  
  ids=ids[ids!=0]%>%unlist
  return(ids)
}

gen_eventag<-function(y){
  n<-length(y)
  max<-max(y)
  
  temp=matrix(0,n,max)
  eventag=rep(0,1)
  for (i in 1:n){
    if (y[i] != 1 ){
      for (j in 1:max(y[i])){
        temp[i,j]=1
      }
    } else{
      temp[i,1]=2
    }
    eventag=c(eventag,temp[i,])
  }
  
  eventag=eventag[eventag!=0]%>%unlist
  eventag[which(eventag==2)]=0
  return(eventag)
}

gen_event<-function(y){
  n<-length(y)
  max<-max(y)
  
  temp=matrix(0,n,max)
  event=rep(0,1)
  for (i in 1:n){
    #if (y[i] != 0){
    for (j in 1:max(y[i])){
      temp[i,j]=j
    }
    #} else{
    #  temp[i,1]=max+1
    #}
    event=c(event,temp[i,])
  }
  
  event=event[event!=0]%>%unlist
  event[which(event==max+1)]=0
  return(event)
}

gen_pp<-function(rate,y){
  n<-length(y)
  max<-max(y)
  
  #poisson process generation
  jump<-matrix(0,nrow=1,ncol=max+1)
  jump<-data.frame(jump)
  for (i in 1:n){
    path<-round(pp.sim(rate,y[i]))
    a<-t(data.frame(path))
    a<-data.frame(a)
    jump<-rbind.fill(jump,a)
  }
  
  for (i in 1:(n+1)){
    if (is.na(jump[i,ncol(jump)]) == FALSE & is.na(jump[i,2]) == TRUE){
      jump[i,2]<-jump[i,ncol(jump)]
    }
  }
  jump<-jump[-1,-ncol(jump)][,-1]
  return(jump)
}

gen_pp_nh<-function(rate,y){
  n<-length(y)
  max<-max(y)
  
  #poisson process generation
  jump<-matrix(0,nrow=1,ncol=max+1)
  jump<-data.frame(jump)
  for (i in 1:n){
    path<-round(pp.sim(rate[i],y[i]))
    a<-t(data.frame(path))
    a<-data.frame(a)
    jump<-rbind.fill(jump,a)
  }
  
  for (i in 1:(n+1)){
    if (is.na(jump[i,ncol(jump)]) == FALSE & is.na(jump[i,2]) == TRUE){
      jump[i,2]<-jump[i,ncol(jump)]
    }
  }
  jump<-jump[-1,-ncol(jump)][,-1]
  return(jump)
}

#gen_time<-function(jump){
#  temp=jump


#time=rep(0,1)
#for (i in 1:100){
#  time=c(time,temp[i,])
#}
#time=time[is.na(time)==F]%>%unlist
#time=time[-1]
#time=data.frame(time)
#time=time$time
#return(time)
#}

gen_time_nh<-function(jump){
  temp=jump
  #time=rep(0,1)
  time=NULL
  for (i in 1:100){
    time=c(time,temp[i,])
  }
  #time=time[time!=500]%>%unlist
  time=time[!is.na(time)]%>%unlist
  #time=time[-1]
  time=data.frame(time)
  time=time$time
  return(time)
}

gen_t0<-function(id,time_end,event){
  data=data.frame(event,time_end,id)
  
  t0=rep(0,length(id))
  for (i in 2:length(t0)){
    if (data$id[i] == data$id[i-1]){
      t0[i]=data$time_end[i-1]
    }
  }
  
  return(t0)
}

gen_covariate<-function(age,gender,treatmentA,y){
  #arrange covariates
  age1=rep(0,1)
  treatment=rep(0,1)
  gender1=rep(0,1)
  n=length(y)
  agetemp=matrix(0,n,max(y))
  trttemp=matrix(0,n,max(y))
  gentemp=matrix(0,n,max(y))
  
  treatmentA[which(treatmentA==0)]=2
  gender[which(gender==0)]=2
  for (i in 1:n){
    if (y[i]!=0){
      for (j in 1:max(y[i])){
        agetemp[i,j]=age[i]
        trttemp[i,j]=treatmentA[i]
        gentemp[i,j]=gender[i]
      }
    } else{
      agetemp[i,1]=age[i]
      trttemp[i,1]=treatmentA[i]
      gentemp[i,1]=gender[i]
    }
    age1=c(age1,agetemp[i,])
    treatment=c(treatment,trttemp[i,])
    gender1=c(gender1,gentemp[i,])
  }
  
  age1=age1[age1!=0]%>%unlist
  treatment= treatment[treatment!=0]%>%unlist
  treatment[which(treatment==2)]=0
  gender1=gender1[gender1!=0]%>%unlist
  gender1[which(gender1==2)]=0
  return(data.frame(age1,gender1,treatment))
}

gen_covariate_pwp<-function(age,gender,treatmentA,y,time_end){
  #arrange covariates
  age1=rep(0,1)
  treatment=rep(0,1)
  gender1=rep(0,1)
  n=length(y)
  agetemp=matrix(0,n,max(y))
  trttemp=matrix(0,n,max(y))
  gentemp=matrix(0,n,max(y))
  
  treatmentA[which(treatmentA==0)]=2
  gender[which(gender==0)]=2
  for (i in 1:n){
    if (y[i]!=0){
      for (j in 1:max(y[i])){
        agetemp[i,j]=age[i]
        trttemp[i,j]=1
        gentemp[i,j]=gender[i]
      }
      trttemp[i,1]=treatmentA[i]
    } else{
      agetemp[i,1]=age[i]
      trttemp[i,1]=treatmentA[i]
      gentemp[i,1]=gender[i]
    }
    age1=c(age1,agetemp[i,])
    treatment=c(treatment,trttemp[i,])
    gender1=c(gender1,gentemp[i,])
  }
  
  age1=age1[age1!=0]%>%unlist
  age1=age1+time_end/360
  treatment= treatment[treatment!=0]%>%unlist
  treatment[which(treatment==2)]=0
  gender1=gender1[gender1!=0]%>%unlist
  gender1[which(gender1==2)]=0
  return(data.frame(age1,gender1,treatment))
}


test_one<-function(n,beta,r){
  dt<-NULL
  dt<-tryCatch({gen_Data_ag(n,beta,r)},
               error = function(e) {
                 NA
               })
  
  #dt<-gen_Data_ag(n,beta,r)
  if (!is.na(dt)){
    data<-data.frame(id=dt$id,time_start=dt$time_start,
                     time_end=dt$time_end,event=dt$event_ag,eventnums=dt$event,
                     age=dt$age_ag,gender=dt$gender_ag,
                     A=dt$treatmentA_ag)
    data1=data[!duplicated(data[,1]),]
    data1=data.frame(data1,status=dt$status)
    data_count<-data.frame(y=dt$y_run,OS=dt$os,Age=dt$age,
                           Gender=dt$gender,A=dt$treatmentA)
    ag<-ag_one(data)
    cox<-cox_one(data1)
    pwp<-pwp_one(data)
    poi<-poi.one(data_count)
    poi_off<-poi.offset(data_count)
    over<-over.one(data_count)
    over_off<-over.offset(data_count)
    nb<-nb.one(data_count)
    nb_off<-nb.offset(data_count)
  } else {
    ag<-NA
    cox<-NA
    pwp<-NA
    poi<-NA
    poi_off<-NA
    over<-NA
    over_off<-NA
    nb<-NA
    nb_off<-NA
  }     
  
  return(list(ag=ag,cox=cox,pwp=pwp,poi=poi,poi_off=poi_off,
              over=over,over_off=over_off,nb=nb,nb_off=nb_off))
}

get_power<-function(nrep,n,beta,r,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_one(n,beta,r))
  pow_ag<-mean(test1[1,]<level,na.rm=T)
  pow_cox<-mean(test1[2,]<level,na.rm=T)
  pow_pwp<-mean(test1[3,]<level,na.rm=T)
  pow_poi<-mean(test1[4,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[5,]<level,na.rm=T)
  pow_over<-mean(test1[6,]<level,na.rm=T)
  pow_over_offset<-mean(test1[7,]<level,na.rm=T)
  pow_nb<-mean(test1[8,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[9,]<level,na.rm=T)
  nanum<-sum(is.na(test1))
  return(list(pow_ag=pow_ag,pow_cox=pow_cox,pow_pwp=pow_pwp,
              pow_poi=pow_poi,pow_poi_offset=pow_poi_offset,
              pow_over=pow_over,pow_over_offset=pow_over_offset,
              pow_nb=pow_nb,pow_nb_offset=pow_nb_offset,nanum=nanum))
}

test_one_mcar<-function(n,beta,r){
  dt<-NULL
  dt<-tryCatch({gen_Data_mcar(n,beta,r)},
               error = function(e) {
                 NA
               })
  
  #dt<-gen_Data_ag(n,beta,r)
  if (!is.na(dt)){
    data<-data.frame(id=dt$id,time_start=dt$time_start,
                     time_end=dt$time_end,event=dt$event_ag,eventnums=dt$event,
                     age=dt$age_ag,gender=dt$gender_ag,
                     A=dt$treatmentA_ag)
    data1=data[!duplicated(data[,1]),]
    data1=data.frame(data1,status=dt$status)
    data_count<-data.frame(y=dt$y_run,OS=dt$os,Age=dt$age,
                           Gender=dt$gender,A=dt$treatmentA)
    ag<-ag_one(data)
    cox<-cox_one(data1)
    pwp<-pwp_one(data)
    poi<-poi.one(data_count)
    poi_off<-poi.offset(data_count)
    over<-over.one(data_count)
    over_off<-over.offset(data_count)
    nb<-nb.one(data_count)
    nb_off<-nb.offset(data_count)
  } else {
    ag<-NA
    cox<-NA
    pwp<-NA
    poi<-NA
    poi_off<-NA
    over<-NA
    over_off<-NA
    nb<-NA
    nb_off<-NA
  }     
  
  return(list(ag=ag,cox=cox,pwp=pwp,poi=poi,poi_off=poi_off,
              over=over,over_off=over_off,nb=nb,nb_off=nb_off))
}

get_power_mcar<-function(nrep,n,beta,r,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_one_mcar(n,beta,r))
  pow_ag<-mean(test1[1,]<level,na.rm=T)
  pow_cox<-mean(test1[2,]<level,na.rm=T)
  pow_pwp<-mean(test1[3,]<level,na.rm=T)
  pow_poi<-mean(test1[4,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[5,]<level,na.rm=T)
  pow_over<-mean(test1[6,]<level,na.rm=T)
  pow_over_offset<-mean(test1[7,]<level,na.rm=T)
  pow_nb<-mean(test1[8,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[9,]<level,na.rm=T)
  nanum<-sum(is.na(test1))
  return(list(pow_ag=pow_ag,pow_cox=pow_cox,pow_pwp=pow_pwp,
              pow_poi=pow_poi,pow_poi_offset=pow_poi_offset,
              pow_over=pow_over,pow_over_offset=pow_over_offset,
              pow_nb=pow_nb,pow_nb_offset=pow_nb_offset,nanum=nanum))
}

ag_one<-function(data){
  mod_ag <- tryCatch({coxph(Surv(time_start, time_end, event) ~ 
                              gender + age + A +cluster(id),
                            data=data)},
                 error=function(e){
                   NA
                   #glm.nb(y~Age+Gender+A, data=data_count)
                 })
  if (is.na(mod_ag)){
    p.val<-NA
  } else {
    p.val <- tryCatch({summary(mod_ag)$coefficient[3,6]},
                      error=function(e){NA})
  }
  return(p.val)
}

ag_simulation<-function(nrep,n,beta,r,level=0.05){
  test<- replicate(nrep, test_one(n,beta,r)$ag)
  pow=mean(test<level)
  return(pow)
}

cox_one<-function(data1){
  mod<-coxph(Surv(time_end, status) ~ 
               gender + age + A,
             data=data1)
  return(summary(mod)$coefficient[3,5])
}

cox_simulation<-function(nrep,n,beta,r,level=0.05){
  test<- replicate(nrep, test_one(n,beta,r)$cox)
  pow=mean(test<level)
  return(pow)
}

gen_e<-function(age_ag){
  n<-length(age_ag)
  e<-rep(1,n)
}

pwp_one<-function(data){
  mod_pwp <- tryCatch({coxph(Surv(time_start,time_end,event)~
                      age+gender+A+cluster(id)+strata(eventnums), data=data)},
                     error=function(e){
                       NA
                       #glm.nb(y~Age+Gender+A, data=data_count)
                     })
  if (is.na(mod_pwp)){
    p.val<-NA
  } else {
    p.val <- tryCatch({summary(mod_pwp)$coefficient[3,6]},
                      error=function(e){NA})
  }
  return(p.val)
 # mod_pwp<-coxph(Surv(time_start,time_end,event)~
  #                 age+gender+A+cluster(id)+strata(eventnums), data=data)
  #return(summary(mod_pwp)$coefficient[3,6])
}


pwp_simulation<-function(nrep,n,beta,r,level=0.05){
  test<- replicate(nrep, test_one(n,beta,r)$pwp)
  pow=mean(test<level)
  return(pow)
}

#pwp_simulation_p<-function(nrep,n,beta,r_pwp,level=0.05){
#  test<- replicate(nrep, pwp_one_p(n,beta,r_pwp))
#  pow=mean(test<level)
#  return(pow)
#}

gen_timing<-function(n,rate,censor){
  jump<-matrix(0,nrow=1,ncol=40)
  jump<-data.frame(jump)
  off=NULL
  y<-matrix(0,nrow=n,ncol=1)
  
  for (i in 1:n){
    event<-pp.sim(rate[i],40)
    y[i]<-max(which(event<=180))-1
    number<-y[i]+1
    b<-t(data.frame(event[1:number]))
    b<-data.frame(b)
    jump<-rbind.fill(jump,b)
  }
  
  for (i in 1:(n+1)){
    if (is.na(jump[i,ncol(jump)]) == FALSE & is.na(jump[i,2]) == TRUE){
      jump[i,2]<-jump[i,ncol(jump)]
    }
  }
  
  jump<-jump[-1,]
  jump<-subset(jump,select=-c(X1,b))
  jump[is.na(jump)]<-0
  
  y_c=y
  y_new=y
  
  status=rep(1,n)#1 uncensored, 0 censoreddata
  for (i in 1:n){
    if (jump[i,1]==0||is.na(jump[i,1])){
      jump[i,1]<-censor[i]
      y_c[i]<-1
      if(censor[i]!=180){
        status[i]<-0
      }
    } else if (mean(jump[i,]>censor[i])!=0 && !is.na(mean(jump[i,]>censor[i]))) {
      c<-which(jump[i,]>censor[i])
      status[i]<-0
      if (is.na(c) || c==0){
        jump[i,1]<-censor[i]
        y_c[i]<-1
        y_new[i]<-0
        #jump[i,2:ncol(jump)]<-0
      } else {
        jump[i,c[1]]<-censor[i]
        y_c[i]<-c[1]
        y_new[i]<-c[1]
        #jump[i,(y_c[i]+1):ncol(jump)]<-0
      }
      jump[i,(y_c[i]+1):ncol(jump)]<-0
    } 
  }
  for (i in 1:100){
    off[i]<-max(jump[i,])
  }
  
  jump[jump==0]<-NA
  
  return(list(y=y,y_c=y_c,y_new=y_new,jump=jump,status=status,off=off))
}


gen_pwp<-function(n,rate_pwp){
  jump<-matrix(0,nrow=1,ncol=11)
  jump<-data.frame(jump)
  y<-matrix(0,nrow=n,ncol=1)
  for (i in 1:n){
    event<-matrix(0,nrow=1,ncol=11)
    for (j in 2:11){
      temp<-pp.sim(rate_pwp[i,j-1],1)
      event[j]<-temp[2]
    }
    event=cumsum(event)
    y[i]<-max(which(event<=90))-1
    number<-y[i]+1
    b<-t(data.frame(event[1:number]))
    b<-data.frame(b)
    jump<-rbind.fill(jump,b)
  }
  max<-max(y)+2
  for (i in 1:(n+1)){
    if (is.na(jump[i,ncol(jump)]) == FALSE & is.na(jump[i,2]) == TRUE){
      jump[i,2]<-jump[i,ncol(jump)]
    }
  }
  jump<-jump[-1,c(-max:-ncol(jump))][,-1]
  return(list(y=y,jump=jump))
}

poi.offset = function(data_count){
  #dt<-gen_Data_ag(n,beta,r)
  #data_count<-data.frame(y=dt$y_c,age=dt$age,gender=dt$gender,A=dt$treatmentA)
  fit <- glm(y~Age+Gender+A,offset=log(OS), data=data_count,family=poisson)
  p.val <- summary(fit)$coefficients[4,4]
  return(p.val)
}

poi.glmm=function(data_count){
  rowid=row.names(data_count)
  data_count<-data.frame(data_count,rowid=rowid)
  fit <- tryCatch({glmer(y~Age+Gender+A+(1|rowid), offset=log(OS),family=poisson,data=data_count)},
                 error=function(e){
                   NA
                   #glm.nb(y~Age+Gender+A, data=data_count)
                 })
  if (is.na(fit)){
    p.val<-NA
  } else {
    p.val <- summary(fit)$coefficients[4,4]}
  
  return(p.val)
}

poi.one = function(data_count){
  #dt<-gen_Data_ag(n,beta,r)
  #data_count<-data.frame(y=dt$y_c,age=dt$age,gender=dt$gender,A=dt$treatmentA)
  fit <- glm(y~Age+Gender+A, data=data_count,family=poisson)
  p.val <- summary(fit)$coefficients[4,4]
  return(p.val)
}

over.offset<-function(data_count){
  #dt<-gen_Data_ag(n,beta,r)
  #data_count<-data.frame(y=dt$y_c,age=dt$age,gender=dt$gender,A=dt$treatmentA)
  quasi <- glm(y~Age+Gender+A,offset = log(OS), data=data_count,family=quasipoisson)
  p.val <- summary(quasi)$coefficients[4,4]
  return(p.val)
}

over.one<-function(data_count){
  #dt<-gen_Data_ag(n,beta,r)
  #data_count<-data.frame(y=dt$y_c,age=dt$age,gender=dt$gender,A=dt$treatmentA)
  quasi <- glm(y~Age+Gender+A, data=data_count,family=quasipoisson)
  p.val <- summary(quasi)$coefficients[4,4]
  return(p.val)
}

nb.offset<-function(data_count){
  #dt<-gen_Data_ag(n,beta,r)
  #data_count<-data.frame(y=dt$y_c,age=dt$age,gender=dt$gender,A=dt$treatmentA)
  #nb<-glm.nb(y~Age+Gender+A+offset(log(OS)), data=data_count)
  #p.val <- summary(nb)$coefficients[4,4]
  nb <- tryCatch({glm.nb(y~Age+Gender+A+offset(log(OS)), data=data_count)},
                 error=function(e){
                   NA
                   #glm.nb(y~Age+Gender+A, data=data_count)
                 })
  if (is.na(nb)){
    p.val<-NA
  } else {
    p.val <- summary(nb)$coefficients[4,4]
  }
  
  return(p.val)
}

nb.one<-function(data_count){
  #dt<-gen_Data_ag(n,beta,r)
  #data_count<-data.frame(y=dt$y_c,age=dt$age,gender=dt$gender,A=dt$treatmentA)
  #nb <- glm.nb(y~Age+Gender+A, data=data_count)
  #p.val <- summary(nb)$coefficients[4,4]
  #return(p.val)
  nb <- tryCatch({glm.nb(y~Age+Gender+A, data=data_count)},
                 error=function(e){
                   NA
                   #glm.nb(y~Age+Gender+A, data=data_count)
                 })
  if (is.na(nb)){
    p.val<-NA
  } else {
    p.val <- summary(nb)$coefficients[4,4]
  }
  
}

calc.pow.poi = function(nrep, n, beta, r, level = 0.05){
  trials <- replicate(nrep, poi.one(n, beta, r))
  pow <- mean(trials < level)
  return(pow)
}

calc.pow.over = function(nrep, n, beta, r, level = 0.05){
  trials <- replicate(nrep, over.one(n, beta, r))
  pow <- mean(trials < level)
  return(pow)
}

calc.pow.nb = function(nrep, n, beta, r, level = 0.05){
  trials <- replicate(nrep, nb.one(n, beta,r))
  pow <- mean(trials < level)
  return(pow)
}

gen_censor=function(x,rc){
  censor=NULL
  n=length(x)
  for (i in 1:n){
    censor[i]<-min(180,rexp(1,rc))
  }
  
  return(censor=censor)
} 

gen_censor2=function(x,s){
  censor=NULL
  n=length(x)
  for (i in 1:n){
    rc=x[i]
    c<-rc/s
    censor[i]<-min(180,rexp(1,c))
  }
  censor[is.na(censor)]<-180
  return(censor=censor)
} 

gen_e<-function(age_ag){
  n<-length(age_ag)
  e<-rep(1,n)
}

test_censor<-function(x,s){
  censor<-gen_censor2(x,s)
  n<-length(x)
  process<-gen_timing(n,x,censor)
  mean<-1-sum(process$status)/100
  return(cen=mean)
}

try_censor<-function(t,x,s){
  trials <- replicate(t, test_censor(x,s))
  get<-mean(trials)
  return(get)
}

run_test<-function(n,r,c){
  process<-NULL
  process<-gen_timing(n,r,c)
  test_event<-gen_event(process$y_c)
  test_id<-gen_id(process$y_c)
  test_end<-gen_time_nh(process$jump)
  len1<-length(test_event)
  len2<-length(test_end)
  len3<-length(test_id)
  return(list(test_id=test_id,y=process$y_c,te=test_event,td=test_end,
              jump=process$jump,len1=len1,len2=len2,len3=len3))
}

run_mcar<-function(nrep,x){
  result<-NULL
  i<-0
  for (rc in seq(0.001,0.01,0.001)){
    i=i+1
    result[i]<-rep_mcar(nrep,x,rc)
  }
  return(result)
}

rate_mcar<-function(x,rc){
  ce<-gen_censor(x,rc)
  process<-gen_timing(100,x,ce)
  rmc<-1-sum(process$status)/100
  return(rmc)
}

rep_mcar<-function(nrep,x,rc){
  repp<-NULL
  repp<-replicate(nrep,rate_mcar(x,rc))
  return(mean(repp,na.rm=T))
}

rate_mar<-function(x,rc){
  ce<-gen_censor(x,rc)
  process<-gen_timing(100,x,ce)
  rmc<-1-sum(process$status)/100
  return(rmc)
}

rep_mcar<-function(nrep,x,rc){
  repp<-NULL
  repp<-replicate(nrep,rate_mcar(x,rc))
  return(mean(repp,na.rm=T))
}

gen_timing_nc<-function(n,rate){
  jump<-matrix(0,nrow=1,ncol=40)
  jump<-data.frame(jump)
  off=NULL
  y<-matrix(0,nrow=n,ncol=1)
  
  for (i in 1:n){
    event<-pp.sim(rate[i],40)
    y[i]<-max(which(event<=180))-1
    #y[i]<-max(which(event>0))
    number<-y[i]+1
    b<-t(data.frame(event[1:number]))
    b<-data.frame(b)
    jump<-rbind.fill(jump,b)
  }
  
  for (i in 1:(n+1)){
    if (is.na(jump[i,ncol(jump)]) == FALSE & is.na(jump[i,2]) == TRUE){
      jump[i,2]<-jump[i,ncol(jump)]
    }
  }
  
  jump<-jump[-1,]
  jump<-subset(jump,select=-c(X1,b))
  #jump<-subset(jump,select=-c(X12))
  jump[is.na(jump)]<-0
  y_c=y
  y_c[y_c==0]<-1
  status=y#1 uncensored, 0 censorddata
  status[status!=0]<-1
  
  for(i in 1:n){
    if (jump[i,1]==0){
      jump[i,1]<-180
    }
  }
  for (i in 1:100){
    off[i]<-max(jump[i,])
  }
  jump[jump==0]<-NA
  return(list(y=y,y_c=y_c,jump=jump,status=status,off=off))
}

gen_Data_nc <- function(n, beta, r){
  age<-NULL
  treatmentA<-NULL
  gender<-NULL
  rate_run<-NULL
  censored<-NULL
  process<-NULL
  jump<-NULL
  y_c<-NULL
  y_run<-NULL
  id<-NULL
  event_ag<-NULL
  event<-NULL
  time_end<-NULL
  time_start<-NULL
  ag_covariates<-NULL
  age_ag<-NULL
  gender_ag<-NULL
  treatmentA_ag<-NULL
  status<-NULL
  os<-NULL
  
  age<-round(runif(n=n,min=2,max=12))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  rate_run<-r*exp(beta[1]+beta[2]*age+beta[3]*gender+beta[4]*treatmentA)
  
  process<-gen_timing_nc(n=100,rate_run)
  jump<-process$jump
  y_c<-process$y_c
  y_run<-process$y
  os<-process$off
  status<-process$status
  
  id <- gen_id(y_c)
  event_ag <-gen_eventag(y_c)
  event <- gen_event(y_c)
  time_end <- gen_time_nh(jump=jump)
  time_start <-gen_t0(id,time_end,event)
  
  ag_covariates <-gen_covariate(age,gender,treatmentA,y_c)
  age_ag <-ag_covariates$age1
  gender_ag <-ag_covariates$gender1
  treatmentA_ag <-ag_covariates$treatment
  
  return(list(age=age, gender=gender, treatmentA=treatmentA,
              y_c=y_c,y_run=y_run,id=id, status=status, event=event, event_ag=event_ag,
              time_end=time_end, time_start=time_start, age_ag=age_ag,
              gender_ag=gender_ag, treatmentA_ag=treatmentA_ag,os=os))
}

test_one_nc<-function(n,beta,r){
  dt<-NULL
  dt<-gen_Data_nc(n,beta,r)
  
  #if (!is.na(dt)){
  data<-data.frame(id=dt$id,time_start=dt$time_start,
                   time_end=dt$time_end,event=dt$event_ag,eventnums=dt$event,
                   age=dt$age_ag,gender=dt$gender_ag,
                   A=dt$treatmentA_ag)
  data1=data[!duplicated(data[,1]),]
  data1=data.frame(data1,status=dt$status)
  data_count<-data.frame(y=dt$y_run,OS=dt$os,Age=dt$age,
                         Gender=dt$gender,A=dt$treatmentA)
  ag<-ag_one(data)
  cox<-cox_one(data1)
  pwp<-pwp_one(data)
  poi<-poi.one(data_count)
  poi_off<-poi.offset(data_count)
  over<-over.one(data_count)
  over_off<-over.offset(data_count)
  nb<-nb.one(data_count)
  nb_off<-nb.offset(data_count)
  #} else {
  
  return(list(ag=ag,cox=cox,pwp=pwp,poi=poi,poi_off=poi_off,
              over=over,over_off=over_off,nb=nb,nb_off=nb_off))
}

get_power_nc<-function(nrep,n,beta,r,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_one_nc(n,beta,r))
  pow_ag<-mean(test1[1,]<level,na.rm=T)
  pow_cox<-mean(test1[2,]<level,na.rm=T)
  pow_pwp<-mean(test1[3,]<level,na.rm=T)
  pow_poi<-mean(test1[4,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[5,]<level,na.rm=T)
  pow_over<-mean(test1[6,]<level,na.rm=T)
  pow_over_offset<-mean(test1[7,]<level,na.rm=T)
  pow_nb<-mean(test1[8,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[9,]<level,na.rm=T)
  nanum<-sum(is.na(test1))
  return(list(pow_ag=pow_ag,pow_cox=pow_cox,pow_pwp=pow_pwp,
              pow_poi=pow_poi,pow_poi_offset=pow_poi_offset,
              pow_over=pow_over,pow_over_offset=pow_over_offset,
              pow_nb=pow_nb,pow_nb_offset=pow_nb_offset,nanum=nanum))
}

gen_pp_of<-function(n,beta,r){
  age<-round(runif(n=n,min=2,max=12))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  rate_run<-r*exp(beta[1]+beta[2]*age+beta[3]*gender+beta[4]*treatmentA)
  
  total_time<-get_total_time(n,rate_run)
  count_ex<-get_ppoffset(total_time,rate_run)
  
  return(list(age=age,treatmentA=treatmentA,gender=gender,rate_run=rate_run,
              total_time=total_time,count_ex=count_ex))
}

gen_pp_of2<-function(n,beta,r){
  age<-round(runif(n=n,min=2,max=12))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  rate_run<-r*exp(beta[1]+beta[2]*age+beta[3]*gender+beta[4]*treatmentA)
  count_ex<-get_ppoffset(180,rate_run)
  return(list(age=age,treatmentA=treatmentA,gender=gender,rate_run=rate_run,
              total_time=total_time,count_ex=count_ex))
}

get_total_time<-function(n,rate){
  total_t<-rep(0,n)
  for (i in 1:n){
    total_t[i]<-min(180,rexp(1,rate[i]))
  }
  return(total_t=total_t)
}

#get_ppoffset<-function(tt,rate){
#  cou_exacer<-rep(0,n)
#  total_ti<-rep(0,n)
#  for(i in 1:n){
#    run_pp<-NA
#    run_pp<-sim_pp1(tt,rate[i])
#    cou_exacer[i]<-length(unlist(run_pp[2]))
#    total_ti[i]<-1
#  }
#  return(cou_exacer=cou_exacer)
#}

get_ppoffset<-function(tt,rate){
  cou_exacer<-rep(0,n)
  for(i in 1:n){
    run_pp<-NA
    run_pp<-sim_pp2(tt[i],rate[i])
    cou_exacer[i]<-length(unlist(run_pp[2]))
  }
  return(cou_exacer=cou_exacer)
}

test_one_off<-function(n,beta,r){
  dt<-NULL
  dt<-gen_pp_of(n,beta,r)
  
  data_count<-data.frame(y=dt$count_ex,OS=dt$total_time,Age=dt$age,
                         Gender=dt$gender,A=dt$treatmentA)
  poi<-poi.one(data_count)
  poi_off<-poi.offset(data_count)
  over<-over.one(data_count)
  over_off<-over.offset(data_count)
  nb<-nb.one(data_count)
  nb_off<-nb.offset(data_count)
  #} else {
  
  return(list(poi=poi,poi_off=poi_off,
              over=over,over_off=over_off,nb=nb,nb_off=nb_off))
}

get_power_off<-function(nrep,n,beta,r,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_one_off(n,beta,r))
  pow_poi<-mean(test1[1,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[2,]<level,na.rm=T)
  pow_over<-mean(test1[3,]<level,na.rm=T)
  pow_over_offset<-mean(test1[4,]<level,na.rm=T)
  pow_nb<-mean(test1[5,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[6,]<level,na.rm=T)
  nanum<-sum(is.na(test1))
  return(list( pow_poi=pow_poi,pow_poi_offset=pow_poi_offset,
              pow_over=pow_over,pow_over_offset=pow_over_offset,
              pow_nb=pow_nb,pow_nb_offset=pow_nb_offset,nanum=nanum))
}

gen_timing<-function(n,rate,censor){
  jump<-matrix(0,nrow=1,ncol=20)
  jump<-data.frame(jump)
  off=NULL
  y<-matrix(0,nrow=n,ncol=1)
  
  for (i in 1:n){
    event<-pp.sim(rate[i],20)
    y[i]<-max(which(event<=180))-1
    number<-y[i]+1
    b<-t(data.frame(event[1:number]))
    b<-data.frame(b)
    jump<-rbind.fill(jump,b)
  }
  
  for (i in 1:(n+1)){
    if (is.na(jump[i,ncol(jump)]) == FALSE & is.na(jump[i,2]) == TRUE){
      jump[i,2]<-jump[i,ncol(jump)]
    }
  }
  
  jump<-jump[-1,]
  jump<-subset(jump,select=-c(X1,b))
  jump[is.na(jump)]<-0
  
  y_c=y
  y_new=y
  
  status=rep(1,n)#1 uncensored, 0 censoreddata
  for (i in 1:n){
    if (jump[i,1]==0||is.na(jump[i,1])){
      jump[i,1]<-censor[i]
      y_c[i]<-1
      if(censor[i]!=180){
        status[i]<-0
      }
    } else if (mean(jump[i,]>censor[i])!=0 && !is.na(mean(jump[i,]>censor[i]))) {
      c<-which(jump[i,]>censor[i])
      status[i]<-0
      if (is.na(c) || c==0){
        jump[i,1]<-censor[i]
        y_c[i]<-1
        y_new[i]<-0
        #jump[i,2:ncol(jump)]<-0
      } else {
        jump[i,c[1]]<-censor[i]
        y_c[i]<-c[1]
        y_new[i]<-c[1]
        #jump[i,(y_c[i]+1):ncol(jump)]<-0
      }
      jump[i,(y_c[i]+1):ncol(jump)]<-0
    } 
  }
  for (i in 1:100){
    off[i]<-max(jump[i,])
  }
  
  jump[jump==0]<-NA
  
  return(list(y=y,y_c=y_c,y_new=y_new,jump=jump,status=status,off=off))
}


get_ppoffset_time<-function(rate){
  cou_exacer<-rep(0,n)
  time_offset<-rep(0,n)
  for(i in 1:n){
    run_pp<-NA
    run_pp<-generate_pp(180,rate[i])
    cou_exacer[i]<-unlist(run_pp[2])
    time_offset[i]<-max(unlist(run_pp[3]))
  }
  return(list(cou_exacer=cou_exacer,time_offset=time_offset))
}

gen_pp_of_time<-function(n,beta,r){
  age<-round(runif(n=n,min=2,max=12))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  rate_run<-r*exp(beta[1]+beta[2]*age+beta[3]*gender+beta[4]*treatmentA)
  
  result<-get_ppoffset_time(rate_run)
  total_time<-result$time_offset
  count_ex<-result$cou_exacer
  return(list(age=age,treatmentA=treatmentA,gender=gender,rate_run=rate_run,
              total_time=total_time,count_ex=count_ex))
}

test_one_off_time<-function(n,beta,r){
  dt<-NULL
  dt<-gen_pp_of_time(n,beta,r)
  
  data_count<-data.frame(y=dt$count_ex,OS=dt$total_time,Age=dt$age,
                         Gender=dt$gender,A=dt$treatmentA)
  poi<-poi.one(data_count)
  poi_off<-poi.offset(data_count)
  over<-over.one(data_count)
  over_off<-over.offset(data_count)
  nb<-nb.one(data_count)
  nb_off<-nb.offset(data_count)
  #} else {
  
  return(list(poi=poi,poi_off=poi_off,
              over=over,over_off=over_off,nb=nb,nb_off=nb_off))
}

get_power_off_time<-function(nrep,n,beta,r,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_one_off_time(n,beta,r))
  pow_poi<-mean(test1[1,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[2,]<level,na.rm=T)
  pow_over<-mean(test1[3,]<level,na.rm=T)
  pow_over_offset<-mean(test1[4,]<level,na.rm=T)
  pow_nb<-mean(test1[5,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[6,]<level,na.rm=T)
  nanum<-sum(is.na(test1))
  return(list( pow_poi=pow_poi,pow_poi_offset=pow_poi_offset,
               pow_over=pow_over,pow_over_offset=pow_over_offset,
               pow_nb=pow_nb,pow_nb_offset=pow_nb_offset,nanum=nanum))
}

test_20<-function(n,beta,r){
  dt<-NULL
  dt<-gen_Data_t(n,beta,r)
  
  data_count<-data.frame(y=dt$y_run,OS=dt$os,Age=dt$age,
                         Gender=dt$gender,A=dt$treatmentA)
  poi<-poi.one(data_count)
  poi_off<-poi.offset(data_count)
  over<-over.one(data_count)
  over_off<-over.offset(data_count)
  nb<-nb.one(data_count)
  nb_off<-nb.offset(data_count)
  #} else {
  
  return(list(poi=poi,poi_off=poi_off,
              over=over,over_off=over_off,nb=nb,nb_off=nb_off))
}

get_power_20<-function(nrep,n,beta,r,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_20(n,beta,r))
  pow_poi<-mean(test1[1,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[2,]<level,na.rm=T)
  pow_over<-mean(test1[3,]<level,na.rm=T)
  pow_over_offset<-mean(test1[4,]<level,na.rm=T)
  pow_nb<-mean(test1[5,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[6,]<level,na.rm=T)
  nanum<-sum(is.na(test1))
  return(list( pow_poi=pow_poi,pow_poi_offset=pow_poi_offset,
               pow_over=pow_over,pow_over_offset=pow_over_offset,
               pow_nb=pow_nb,pow_nb_offset=pow_nb_offset,nanum=nanum))
}

gen_timing_t<-function(n,t,rate,censor){
  jump<-matrix(0,nrow=1,ncol=40)
  jump<-data.frame(jump)
  off=NULL
  y<-matrix(0,nrow=n,ncol=1)
  
  for (i in 1:n){
    event<-pp.sim(rate[i],40)
    y[i]<-max(which(event<=t[i]))-1
    number<-y[i]+1
    b<-t(data.frame(event[1:number]))
    b<-data.frame(b)
    jump<-rbind.fill(jump,b)
  }
  
  for (i in 1:(n+1)){
    if (is.na(jump[i,ncol(jump)]) == FALSE & is.na(jump[i,2]) == TRUE){
      jump[i,2]<-jump[i,ncol(jump)]
    }
  }
  
  jump<-jump[-1,]
  jump<-subset(jump,select=-c(X1,b))
  jump[is.na(jump)]<-0
  
  y_c=y
  y_new=y
  
  status=rep(1,n)#1 uncensored, 0 censoreddata
  for (i in 1:n){
    if (jump[i,1]==0||is.na(jump[i,1])){
      jump[i,1]<-censor[i]
      y_c[i]<-1
      if(censor[i]!=180){
        status[i]<-0
      }
    } else if (mean(jump[i,]>censor[i])!=0 && !is.na(mean(jump[i,]>censor[i]))) {
      c<-which(jump[i,]>censor[i])
      status[i]<-0
      if (is.na(c) || c==0){
        jump[i,1]<-censor[i]
        y_c[i]<-1
        y_new[i]<-0
        #jump[i,2:ncol(jump)]<-0
      } else {
        jump[i,c[1]]<-censor[i]
        y_c[i]<-c[1]
        y_new[i]<-c[1]
        #jump[i,(y_c[i]+1):ncol(jump)]<-0
      }
      jump[i,(y_c[i]+1):ncol(jump)]<-0
    } 
  }
  for (i in 1:100){
    off[i]<-max(jump[i,])
  }
  
  jump[jump==0]<-NA
  
  return(list(y=y,y_c=y_c,y_new=y_new,jump=jump,status=status,off=off))
}

gen_Data_t <- function(n, beta, r){
  age<-NULL
  treatmentA<-NULL
  gender<-NULL
  rate_run<-NULL
  censored<-NULL
  process<-NULL
  jump<-NULL
  y_c<-NULL
  y_run<-NULL
  id<-NULL
  event_ag<-NULL
  event<-NULL
  time_end<-NULL
  time_start<-NULL
  ag_covariates<-NULL
  age_ag<-NULL
  gender_ag<-NULL
  treatmentA_ag<-NULL
  status<-NULL
  os<-NULL
  
  age<-round(runif(n=n,min=2,max=12))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  rate_run<-r*exp(beta[1]+beta[2]*age+beta[3]*gender+beta[4]*treatmentA)
  censored<-gen_censor2(rate_run,4)
  total_time<-get_total_time(n,rate_run)
  
  process<-gen_timing_t(n=100,total_time,rate_run,censored)
  jump<-process$jump
  y_c<-process$y_c
  y1<-process$y
  y_run<-process$y_new
  os<-process$off
  status<-process$status
  
  id <- gen_id(y_c)
  event_ag <-gen_eventag(y_c)
  event <- gen_event(y_c)
  time_end <- gen_time_nh(jump=jump)
  time_start <-gen_t0(id,time_end,event)
  
  ag_covariates <-gen_covariate(age,gender,treatmentA,y_c)
  age_ag <-ag_covariates$age1
  gender_ag <-ag_covariates$gender1
  treatmentA_ag <-ag_covariates$treatment
  
  return(list(age=age, gender=gender, treatmentA=treatmentA, y1=y1,
              y_c=y_c,y_run=y_run,id=id, status=status, event=event, event_ag=event_ag,
              time_end=time_end, time_start=time_start, age_ag=age_ag,
              gender_ag=gender_ag, treatmentA_ag=treatmentA_ag,os=os))
}

get_simple<-function(simple_time,rate){
  cou_exacer<-rep(0,n)
  #time_offset<-rep(0,n)
  for(i in 1:n){
    run_pp<-NA
    run_pp<-simple_pp(simple_time[i],rate[i])
    cou_exacer[i]<-run_pp
    #time_offset[i]<-max(unlist(run_pp[2]))
  }
  return(cou_exacer=cou_exacer)
}

gen_simple<-function(n,beta,r){
  age<-round(runif(n=n,min=2,max=12))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  rate_run<-r*exp(beta[1]+beta[2]*age+beta[3]*gender+beta[4]*treatmentA)
  
  total_time<-get_total_time(n,rate_run)
  count_ex<-get_simple(total_time,rate_run)
  
  return(list(age=age,treatmentA=treatmentA,gender=gender,rate_run=rate_run,
              total_time=total_time,count_ex=count_ex))
}

test_simple<-function(n,beta,r){
  dt<-NULL
  dt<-gen_simple(n,beta,r)
  
  data_count<-data.frame(y=dt$count_ex,OS=dt$total_time,Age=dt$age,
                         Gender=dt$gender,A=dt$treatmentA)
  poi<-poi.one(data_count)
  poi_off<-poi.offset(data_count)
  over<-over.one(data_count)
  over_off<-over.offset(data_count)
  nb<-nb.one(data_count)
  nb_off<-nb.offset(data_count)
  #} else {
  
  return(list(poi=poi,poi_off=poi_off,
              over=over,over_off=over_off,nb=nb,nb_off=nb_off))
}

get_power_simple<-function(nrep,n,beta,r,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_simple(n,beta,r))
  pow_poi<-mean(test1[1,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[2,]<level,na.rm=T)
  pow_over<-mean(test1[3,]<level,na.rm=T)
  pow_over_offset<-mean(test1[4,]<level,na.rm=T)
  pow_nb<-mean(test1[5,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[6,]<level,na.rm=T)
  nanum<-sum(is.na(test1))
  return(list( pow_poi=pow_poi,pow_poi_offset=pow_poi_offset,
               pow_over=pow_over,pow_over_offset=pow_over_offset,
               pow_nb=pow_nb,pow_nb_offset=pow_nb_offset,nanum=nanum))
}

get_new<-function(simple_time,rate,n){
  cou_exacer<-rep(0,n)
  time_offset<-matrix(0,nrow=n,ncol=40)
  status<-rep(1,n)
  for(i in 1:n){
    run_pp<-NA
    run_pp<-new_pp2(simple_time[i],rate[i])
    cou_exacer[i]<-unlist(run_pp[1])
    
    if (cou_exacer[i]==0){
      status[i]<-0
    }
    
    len_time<-unlist(run_pp[2])
    for (j in 1:length(len_time)){
      time_offset[i,j]<-len_time[j]
    }
  }
  
  check_end<-time_offset
  
  for (i in 1:n){
    if (max(time_offset[i,])!=simple_time[i]){
      time_offset[i,which.max(time_offset[i,])+1]<-simple_time[i]
    }
  }
  
  #event_re<-rep(1,length(cou_exacer))
  #for (i in 1: length(event_re)){
  #  if (count_ex[i]==1||count_ex[i]==0){
  #    event_re[i]<-0
  #  }
  #}
  count_tte<-rowSums(time_offset!=0)
  newid<-matrix(0,nrow=n,ncol=40)
  t_eventstart<-matrix(0,nrow=n,ncol=40)
  t_someend<-matrix(0,nrow=n,ncol=40)
  event_num<-matrix(0,nrow=n,ncol=40)
  event_recurrent<-matrix(0,nrow=n,ncol=40)
  for (i in 1:n){
     # if (cou_exacer[i]>1){
    #    newid[i,1:cou_exacer[i]]<-i
    #    for (j in 1:cou_exacer[i]){
    #      event_num[i,j]<-j
    #    }
     if (count_tte[i]>1){
        newid[i,1:count_tte[i]]<-i
        for (j in 1:count_tte[i]){
          event_num[i,j]<-j
        }
      } else {
        newid[i,1]<-i
        event_num[i,1]<-1
      }
  }
  
  
  for (i in 1:n){
    #for (j in 1:40){
    if (count_tte[i]>2){
      for (j in 2:count_tte[i]){
        t_eventstart[i,j]<-time_offset[i,j-1]
      }
      t_eventstart[i,(count_tte[i]+1):40]<-NA
    } else if (count_tte[i]==2){
      t_eventstart[i,2]<-time_offset[i,1]
      t_eventstart[i,3:40]<-NA
    } else {
      t_eventstart[i,2:40]<-NA
    }
  }
  for (i in 1:n){
    #for (j in 1:40){
    if (count_tte[i]>=2){
      for (j in 1:count_tte[i]){
        t_someend[i,j]<-time_offset[i,j]
      }
    }  else {
      t_someend[i,1]<-time_offset[i,1]
    }
  }
  
  for(i in 1:n){
    if (count_tte[i]<=1){
      event_recurrent[i,1]<-0
      event_recurrent[i,2:40]<-NA
    } else {
      event_recurrent[i,1:count_tte[i]]<-1
      if (max(check_end[i,])!=simple_time[i]){
        event_recurrent[i,count_tte[i]]<-0
      } 
      event_recurrent[i,(count_tte[i]+1):40]<-NA
    }
  }
  
  new_id=0
  t_eventend=0
  t_start=0
  event_count=0
  event_rec=0
  for(i in 1:n){
    new_id<-c(new_id,newid[i,])
    t_eventend<-c(t_eventend,t_someend[i,])
    t_start<-c(t_start,t_eventstart[i,])
    event_count<-c(event_count,event_num[i,])
    event_rec<-c(event_rec,event_recurrent[i,])
  }
  new_id<-new_id[new_id!=0]%>%unlist
  t_eventend<-t_eventend[t_eventend!=0]%>%unlist
  t_start<-t_start[!is.na(t_start)]%>%unlist
  t_start<-t_start[-1]
  event_rec<-event_rec[!is.na(event_rec)]%>%unlist
  event_rec<-event_rec[-1]
  event_count<-event_count[event_count!=0]%>%unlist
  data.frame(t_start,t_eventend)
  return(list(cou_exacer=cou_exacer,time_offset=time_offset,status=status,
              new_id=new_id,t_eventend=t_eventend,count_tte=count_tte,
              t_start=t_start,event_count=event_count,event_rec=event_rec))
}


gen_new<-function(n,beta,r){
  age<-round(runif(n=n,min=2,max=12))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  #error<-runif(n,-0.001,0.001)
  #error<-rnorm(n,0,0.0005)
  rate_run<-r*exp(beta[1]+beta[2]*age+beta[3]*gender+beta[4]*treatmentA)#+error
  
  total_time<-get_total_time(n,rate_run)
  off_cou<-get_new(total_time,rate_run,n)
  count_ex<-off_cou$cou_exacer
  #offset_run<-off_cou$time_offset
  count_ttevent<-off_cou$count_tte
  events<-off_cou$event_count
  t_to_start<-off_cou$t_start
  t_to_end<-off_cou$t_eventend
  id_new<-off_cou$new_id
  event_reccur<-off_cou$event_rec
  status<-off_cou$status
  ag_covariates <-gen_covariate(age,gender,treatmentA,count_ttevent)
  age_ag <-ag_covariates$age1
  gender_ag <-ag_covariates$gender1
  treatmentA_ag <-ag_covariates$treatment
  
  return(list(age=age,treatmentA=treatmentA,gender=gender,rate_run=rate_run,
              total_time=total_time,count_ex=count_ex,events=events,
              age_ag=age_ag,gender_ag=gender_ag,treatmentA_ag=treatmentA_ag,
              t_to_start=t_to_start,t_to_end=t_to_end,id_new=id_new,
              event_reccur=event_reccur,status=status,count_ttevent=count_ttevent))
}

test_one_offset<-function(n,beta,r){
  dt<-NULL
 # dt<-tryCatch({gen_new(n,beta,r)},
  #               error = function(e) {
  #               NA
  #             })
  dt<-gen_new(n,beta,r)
  #dt<-gen_Data_ag(n,beta,r)
 # if (!is.na(dt)){
    data<-data.frame(id=dt$id_new,time_start=dt$t_to_start,
                     time_end=dt$t_to_end,event=dt$event_reccur,
                     eventnums=dt$events,
                     age=dt$age_ag,gender=dt$gender_ag,
                     A=dt$treatmentA_ag)
    data1=data[!duplicated(data[,1]),]
    data1=data.frame(data1,status=dt$status)
    data_count<-data.frame(y=dt$count_ex,OS=dt$total_time,Age=dt$age,
                           Gender=dt$gender,A=dt$treatmentA)
    ag<-ag_one(data)
    cox<-cox_one(data1)
    pwp<-pwp_one(data)
    
    #poi<-poi.one(data_count)
    poi_off<-poi.offset(data_count)
    #over<-over.one(data_count)
    over_off<-over.offset(data_count)
    #nb<-nb.one(data_count)
    nb_off<-nb.offset(data_count)
 # } else {
#    ag<-NA
#    cox<-NA
#    pwp<-NA
#    poi_off<-NA
#    over_off<-NA
#    nb_off<-NA
#  }     
  
  return(list(ag=ag,pwp=pwp,cox=cox,poi_off=poi_off,
              over_off=over_off,nb_off=nb_off))
}

get_power_offset<-function(nrep,n,beta,r,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_one_offset(n,beta,r))
  pow_ag<-mean(test1[1,]<level,na.rm=T)
  pow_pwp<-mean(test1[2,]<level,na.rm=T)
  pow_cox<-mean(test1[3,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[4,]<level,na.rm=T)
  #pow_over<-mean(test1[5,]<level,na.rm=T)
  pow_over_offset<-mean(test1[5,]<level,na.rm=T)
  #pow_nb<-mean(test1[7,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[6,]<level,na.rm=T)
  nanum<-sum(is.na(test1))
  return(list(pow_ag=pow_ag,pow_pwp=pow_pwp,
              pow_cox=pow_cox,pow_poi_offset=pow_poi_offset,
              pow_over_offset=pow_over_offset,
              pow_nb_offset=pow_nb_offset,nanum=nanum))
}

#start
gen_new_mcar<-function(n,beta,r){
  age<-round(runif(n=n,min=2,max=12))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  error<-rgamma(n,1,300)
  rate_run<-r*exp(beta[1]+beta[2]*age+beta[3]*gender+beta[4]*treatmentA)+error
  
  total_time<-get_total_time_mcar(n,0.001)
  off_cou<-get_new(total_time,rate_run,n)
  count_ex<-off_cou$cou_exacer
  #offset_run<-off_cou$time_offset
  count_ttevent<-off_cou$count_tte
  events<-off_cou$event_count
  t_to_start<-off_cou$t_start
  t_to_end<-off_cou$t_eventend
  id_new<-off_cou$new_id
  event_reccur<-off_cou$event_rec
  status<-off_cou$status
  ag_covariates <-gen_covariate(age,gender,treatmentA,count_ttevent)
  age_ag <-ag_covariates$age1
  gender_ag <-ag_covariates$gender1
  treatmentA_ag <-ag_covariates$treatment
  
  return(list(age=age,treatmentA=treatmentA,gender=gender,rate_run=rate_run,
              total_time=total_time,count_ex=count_ex,events=events,
              age_ag=age_ag,gender_ag=gender_ag,treatmentA_ag=treatmentA_ag,
              t_to_start=t_to_start,t_to_end=t_to_end,id_new=id_new,
              event_reccur=event_reccur,status=status,count_ttevent=count_ttevent))
}

test_one_offset_mcar<-function(n,beta,r){
  dt<-NULL
  # dt<-tryCatch({gen_new(n,beta,r)},
  #               error = function(e) {
  #               NA
  #             })
  dt<-gen_new_mcar(n,beta,r)
  #dt<-gen_Data_ag(n,beta,r)
  # if (!is.na(dt)){
  data<-data.frame(id=dt$id_new,time_start=dt$t_to_start,
                   time_end=dt$t_to_end,event=dt$event_reccur,
                   eventnums=dt$events,
                   age=dt$age_ag,gender=dt$gender_ag,
                   A=dt$treatmentA_ag)
  data1=data[!duplicated(data[,1]),]
  data1=data.frame(data1,status=dt$status)
  data_count<-data.frame(y=dt$count_ex,OS=dt$total_time,Age=dt$age,
                         Gender=dt$gender,A=dt$treatmentA)
  ag<-ag_one(data)
  cox<-cox_one(data1)
  pwp<-pwp_one(data)
  
  #poi<-poi.one(data_count)
  poi_off<-poi.offset(data_count)
  #over<-over.one(data_count)
  over_off<-over.offset(data_count)
  #nb<-nb.one(data_count)
  nb_off<-nb.offset(data_count)
  # } else {
  #    ag<-NA
  #    cox<-NA
  #    pwp<-NA
  #    poi_off<-NA
  #    over_off<-NA
  #    nb_off<-NA
  #  }     
  
  return(list(ag=ag,pwp=pwp,cox=cox,poi_off=poi_off,
              over_off=over_off,nb_off=nb_off))
}

get_power_offset_mcar<-function(nrep,n,beta,r,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_one_offset_mcar(n,beta,r))
  pow_ag<-mean(test1[1,]<level,na.rm=T)
  pow_pwp<-mean(test1[2,]<level,na.rm=T)
  pow_cox<-mean(test1[3,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[4,]<level,na.rm=T)
  #pow_over<-mean(test1[5,]<level,na.rm=T)
  pow_over_offset<-mean(test1[5,]<level,na.rm=T)
  #pow_nb<-mean(test1[7,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[6,]<level,na.rm=T)
  nanum<-sum(is.na(test1))
  return(list(pow_ag=pow_ag,pow_pwp=pow_pwp,
              pow_cox=pow_cox,pow_poi_offset=pow_poi_offset,
              pow_over_offset=pow_over_offset,
              pow_nb_offset=pow_nb_offset,nanum=nanum))
}

get_new_mcar<-function(simple_time,rate){
  cou_exacer<-rep(0,n)
  time_offset<-matrix(0,nrow=n,ncol=40)
  status<-rep(1,n)
  for(i in 1:n){
    run_pp<-NA
    run_pp<-new_pp2(simple_time[i],rate[i])
    cou_exacer[i]<-unlist(run_pp[1])
    
    if (cou_exacer[i]==0){
      status[i]<-0
    }
    
    len_time<-unlist(run_pp[2])
    for (j in 1:length(len_time)){
      time_offset[i,j]<-len_time[j]
    }
  }
  
  for (i in 1:n){
    if (max(time_offset[i,])!=simple_time[i]){
      time_offset[i,which.max(time_offset[i,])+1]<-simple_time[i]
    }
  }
  
  #event_re<-rep(1,length(cou_exacer))
  #for (i in 1: length(event_re)){
  #  if (count_ex[i]==1||count_ex[i]==0){
  #    event_re[i]<-0
  #  }
  #}
  count_tte<-rowSums(time_offset!=0)
  newid<-matrix(0,nrow=n,ncol=40)
  t_eventstart<-matrix(0,nrow=n,ncol=40)
  t_someend<-matrix(0,nrow=n,ncol=40)
  event_num<-matrix(0,nrow=n,ncol=40)
  event_recurrent<-matrix(0,nrow=n,ncol=40)
  for (i in 1:n){
    # if (cou_exacer[i]>1){
    #    newid[i,1:cou_exacer[i]]<-i
    #    for (j in 1:cou_exacer[i]){
    #      event_num[i,j]<-j
    #    }
    if (count_tte[i]>1){
      newid[i,1:count_tte[i]]<-i
      for (j in 1:count_tte[i]){
        event_num[i,j]<-j
      }
    } else {
      newid[i,1]<-i
      event_num[i,1]<-1
    }
  }
  
  for(i in 1:n){
    if (count_tte[i]<=1){
      event_recurrent[i,1]<-0
      event_recurrent[i,2:40]<-NA
    } else {
      event_recurrent[i,1:count_tte[i]]<-1
      event_recurrent[i,(count_tte[i]+1):40]<-NA
    }
  }
  for (i in 1:n){
    #for (j in 1:40){
    if (count_tte[i]>2){
      for (j in 2:count_tte[i]){
        t_eventstart[i,j]<-time_offset[i,j-1]
      }
      t_eventstart[i,(count_tte[i]+1):40]<-NA
    } else if (count_tte[i]==2){
      t_eventstart[i,2]<-time_offset[i,1]
      t_eventstart[i,3:40]<-NA
    } else {
      t_eventstart[i,2:40]<-NA
    }
  }
  for (i in 1:n){
    #for (j in 1:40){
    if (count_tte[i]>=2){
      for (j in 1:count_tte[i]){
        t_someend[i,j]<-time_offset[i,j]
      }
    }  else {
      t_someend[i,1]<-time_offset[i,1]
    }
  }
  new_id=0
  t_eventend=0
  t_start=0
  event_count=0
  event_rec=0
  for(i in 1:n){
    new_id<-c(new_id,newid[i,])
    t_eventend<-c(t_eventend,t_someend[i,])
    t_start<-c(t_start,t_eventstart[i,])
    event_count<-c(event_count,event_num[i,])
    event_rec<-c(event_rec,event_recurrent[i,])
  }
  new_id<-new_id[new_id!=0]%>%unlist
  t_eventend<-t_eventend[t_eventend!=0]%>%unlist
  t_start<-t_start[!is.na(t_start)]%>%unlist
  t_start<-t_start[-1]
  event_rec<-event_rec[!is.na(event_rec)]%>%unlist
  event_rec<-event_rec[-1]
  event_count<-event_count[event_count!=0]%>%unlist
  data.frame(t_start,t_eventend)
  return(list(cou_exacer=cou_exacer,time_offset=time_offset,status=status,
              new_id=new_id,t_eventend=t_eventend,count_tte=count_tte,
              t_start=t_start,event_count=event_count,event_rec=event_rec))
}

get_total_time_mcar<-function(n,rc){
  total_t<-rep(0,n)
  #censor_t<-rep(0,n)
  for (i in 1:n){
    #total_t[i]<-min(180,rexp(1,rate[i]))
    total_t[i]<-min(180,rexp(1,rc))
    #total_t[i]<-min(total_t[i],censor_t[i])
  }
  return(total_t=total_t)
}

#MAR
gen_new_mar<-function(n,beta,r){
  age<-round(runif(n=n,min=2,max=12))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  error<-rgamma(n,1,300)
  rate_run<-r*exp(beta[1]+beta[2]*age+beta[3]*gender+beta[4]*treatmentA)+error
  
  total_time<-get_total_time_mar(n,rate_run,4)
  off_cou<-get_new(total_time,rate_run,n)
  count_ex<-off_cou$cou_exacer
  #offset_run<-off_cou$time_offset
  count_ttevent<-off_cou$count_tte
  events<-off_cou$event_count
  t_to_start<-off_cou$t_start
  t_to_end<-off_cou$t_eventend
  id_new<-off_cou$new_id
  event_reccur<-off_cou$event_rec
  status<-off_cou$status
  ag_covariates <-gen_covariate(age,gender,treatmentA,count_ttevent)
  age_ag <-ag_covariates$age1
  gender_ag <-ag_covariates$gender1
  treatmentA_ag <-ag_covariates$treatment
  
  return(list(age=age,treatmentA=treatmentA,gender=gender,rate_run=rate_run,
              total_time=total_time,count_ex=count_ex,events=events,
              age_ag=age_ag,gender_ag=gender_ag,treatmentA_ag=treatmentA_ag,
              t_to_start=t_to_start,t_to_end=t_to_end,id_new=id_new,
              event_reccur=event_reccur,status=status,count_ttevent=count_ttevent))
}

test_one_offset_mar<-function(n,beta,r){
  dt<-NULL
  # dt<-tryCatch({gen_new(n,beta,r)},
  #               error = function(e) {
  #               NA
  #             })
  dt<-gen_new_mar(n,beta,r)
  #dt<-gen_Data_ag(n,beta,r)
  # if (!is.na(dt)){
  data<-data.frame(id=dt$id_new,time_start=dt$t_to_start,
                   time_end=dt$t_to_end,event=dt$event_reccur,
                   eventnums=dt$events,
                   age=dt$age_ag,gender=dt$gender_ag,
                   A=dt$treatmentA_ag)
  data1=data[!duplicated(data[,1]),]
  data1=data.frame(data1,status=dt$status)
  data_count<-data.frame(y=dt$count_ex,OS=dt$total_time,Age=dt$age,
                         Gender=dt$gender,A=dt$treatmentA)
  ag<-ag_one(data)
  cox<-cox_one(data1)
  pwp<-pwp_one(data)
  
  #poi<-poi.one(data_count)
  poi_off<-poi.offset(data_count)
  #over<-over.one(data_count)
  over_off<-over.offset(data_count)
  #nb<-nb.one(data_count)
  nb_off<-nb.offset(data_count)
  # } else {
  #    ag<-NA
  #    cox<-NA
  #    pwp<-NA
  #    poi_off<-NA
  #    over_off<-NA
  #    nb_off<-NA
  #  }     
  
  return(list(ag=ag,pwp=pwp,cox=cox,poi_off=poi_off,
              over_off=over_off,nb_off=nb_off))
}

get_power_offset_mar<-function(nrep,n,beta,r,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_one_offset_mar(n,beta,r))
  pow_ag<-mean(test1[1,]<level,na.rm=T)
  pow_pwp<-mean(test1[2,]<level,na.rm=T)
  pow_cox<-mean(test1[3,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[4,]<level,na.rm=T)
  #pow_over<-mean(test1[5,]<level,na.rm=T)
  pow_over_offset<-mean(test1[5,]<level,na.rm=T)
  #pow_nb<-mean(test1[7,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[6,]<level,na.rm=T)
  nanum<-sum(is.na(test1))
  return(list(pow_ag=pow_ag,pow_pwp=pow_pwp,
              pow_cox=pow_cox,pow_poi_offset=pow_poi_offset,
              pow_over_offset=pow_over_offset,
              pow_nb_offset=pow_nb_offset,nanum=nanum))
}

get_total_time_mar<-function(n,rate,rc){
  total_t<-rep(0,n)
  #censor_t<-rep(0,n)
  for (i in 1:n){
    #total_t[i]<-min(180,rexp(1,rate[i]))
    total_t[i]<-min(180,rexp(1,rate[i]/rc))
    #total_t[i]<-min(total_t[i],censor_t[i])
  }
  return(total_t=total_t)
}

#non-missing
gen_new_nc<-function(n,beta,r){
  age<-round(runif(n=n,min=14,max=24))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  #error<-rnorm(n,mean=0,sd=0.001)
  error<-rgamma(n,1,300)
  #error<-runif(n,-0.002,0.002)
  rate_run<-r*exp(beta[1]+beta[2]*age_group+beta[3]*gender+beta[4]*treatmentA)+error
  
  total_time<-rep(180,n)
  off_cou<-get_new_nc(total_time,rate_run,n)
  count_ex<-off_cou$cou_exacer
  #offset_run<-off_cou$time_offset
  count_ttevent<-off_cou$count_tte
  events<-off_cou$event_count
  t_to_start<-off_cou$t_start
  t_to_end<-off_cou$t_eventend
  id_new<-off_cou$new_id
  event_reccur<-off_cou$event_rec
  status<-off_cou$status
  ag_covariates <-gen_covariate(age,gender,treatmentA,count_ttevent)
  age_ag <-ag_covariates$age1
  gender_ag <-ag_covariates$gender1
  treatmentA_ag <-ag_covariates$treatment
  
  return(list(age=age,treatmentA=treatmentA,gender=gender,rate_run=rate_run,
              total_time=total_time,count_ex=count_ex,events=events,
              age_ag=age_ag,gender_ag=gender_ag,treatmentA_ag=treatmentA_ag,
              t_to_start=t_to_start,t_to_end=t_to_end,id_new=id_new,
              event_reccur=event_reccur,status=status,count_ttevent=count_ttevent))
}

test_one_offset_nc<-function(n,beta,r){
  dt<-NULL
  # dt<-tryCatch({gen_new(n,beta,r)},
  #               error = function(e) {
  #               NA
  #             })
  dt<-gen_new_nc(n,beta,r)
  #dt<-gen_Data_ag(n,beta,r)
  # if (!is.na(dt)){
  data<-data.frame(id=dt$id_new,time_start=dt$t_to_start,
                   time_end=dt$t_to_end,event=dt$event_reccur,
                   eventnums=dt$events,
                   age=dt$age_ag,gender=dt$gender_ag,
                   A=dt$treatmentA_ag)
  data1=data[!duplicated(data[,1]),]
  data1=data.frame(data1,status=dt$status)
  data_count<-data.frame(y=dt$count_ex,OS=dt$total_time,Age=dt$age,
                         Gender=dt$gender,A=dt$treatmentA)
  ag<-ag_one(data)
  cox<-cox_one(data1)
  pwp<-pwp_one(data)
  
  #poi<-poi.one(data_count)
  poi_off<-poi.offset(data_count)
  #over<-over.one(data_count)
  over_off<-over.offset(data_count)
  #nb<-nb.one(data_count)
  nb_off<-nb.offset(data_count)
  glmm_off<-poi.glmm(data_count)
  # } else {
  #    ag<-NA
  #    cox<-NA
  #    pwp<-NA
  #    poi_off<-NA
  #    over_off<-NA
  #    nb_off<-NA
  #  }     
  
  return(list(ag=ag,pwp=pwp,cox=cox,poi_off=poi_off,
              over_off=over_off,nb_off=nb_off,glmm_off=glmm_off))
}

get_power_offset_nc<-function(nrep,n,beta,r,thres=0.05,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_one_offset_nc(n,beta,r))
  pow_ag<-mean(test1[1,]<level,na.rm=T)
  pow_pwp<-mean(test1[2,]<level,na.rm=T)
  pow_cox<-mean(test1[3,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[4,]<level,na.rm=T)
  pow_poi_offset_adjusted<-mean(test1[4,]<thres,na.rm=T)
  #pow_over<-mean(test1[5,]<level,na.rm=T)
  pow_over_offset<-mean(test1[5,]<level,na.rm=T)
  #pow_nb<-mean(test1[7,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[6,]<level,na.rm=T)
  pow_glmm<-mean(test1[7,]<level,na.rm=T)
  nanum<-sum(is.na(test1))
  return(list(pow_ag=pow_ag,pow_pwp=pow_pwp,
              pow_cox=pow_cox,pow_poi_offset=pow_poi_offset,
              pow_over_offset=pow_over_offset,
              pow_nb_offset=pow_nb_offset,
              pow_poi_offset_adjusted=pow_poi_offset_adjusted,
              pow_glmm=pow_glmm,
              nanum=nanum,
              test_poi=test1[4,]))
}


get_new_nc<-function(simple_time,rate,n){
  cou_exacer<-rep(0,n)
  time_offset<-matrix(0,nrow=n,ncol=40)
  status<-rep(1,n)
  for(i in 1:n){
    run_pp<-NA
    run_pp<-new_pp2(simple_time[i],rate[i])
    cou_exacer[i]<-unlist(run_pp[1])
    
    if (cou_exacer[i]==0){
      status[i]<-0
    }
    
    len_time<-unlist(run_pp[2])
    for (j in 1:length(len_time)){
      time_offset[i,j]<-len_time[j]
    }
  }
  check_end<-time_offset
  
  for (i in 1:n){
    if (max(time_offset[i,])!=simple_time[i]){
      time_offset[i,which.max(time_offset[i,])+1]<-simple_time[i]
    }
  }
  
  #event_re<-rep(1,length(cou_exacer))
  #for (i in 1: length(event_re)){
  #  if (count_ex[i]==1||count_ex[i]==0){
  #    event_re[i]<-0
  #  }
  #}
  count_tte<-rowSums(time_offset!=0)
  newid<-matrix(0,nrow=n,ncol=40)
  t_eventstart<-matrix(0,nrow=n,ncol=40)
  t_someend<-matrix(0,nrow=n,ncol=40)
  event_num<-matrix(0,nrow=n,ncol=40)
  event_recurrent<-matrix(0,nrow=n,ncol=40)
  for (i in 1:n){
    # if (cou_exacer[i]>1){
    #    newid[i,1:cou_exacer[i]]<-i
    #    for (j in 1:cou_exacer[i]){
    #      event_num[i,j]<-j
    #    }
    if (count_tte[i]>1){
      newid[i,1:count_tte[i]]<-i
      for (j in 1:count_tte[i]){
        event_num[i,j]<-j
      }
    } else {
      newid[i,1]<-i
      event_num[i,1]<-1
    }
  }
  
 
  for (i in 1:n){
    #for (j in 1:40){
    if (count_tte[i]>2){
      for (j in 2:count_tte[i]){
        t_eventstart[i,j]<-time_offset[i,j-1]
      }
      t_eventstart[i,(count_tte[i]+1):40]<-NA
    } else if (count_tte[i]==2){
      t_eventstart[i,2]<-time_offset[i,1]
      t_eventstart[i,3:40]<-NA
    } else {
      t_eventstart[i,2:40]<-NA
    }
  }
  for (i in 1:n){
    #for (j in 1:40){
    if (count_tte[i]>=2){
      for (j in 1:count_tte[i]){
        t_someend[i,j]<-time_offset[i,j]
      }
    }  else {
      t_someend[i,1]<-time_offset[i,1]
    }
  }
  
  for(i in 1:n){
    if (count_tte[i]<=1){
      event_recurrent[i,1]<-0
      event_recurrent[i,2:40]<-NA
    } else {
      event_recurrent[i,1:count_tte[i]]<-1
      if (max(check_end[i,])!=simple_time[i]){
        event_recurrent[i,count_tte[i]]<-0
      } 
      event_recurrent[i,(count_tte[i]+1):40]<-NA
    }
  }
  
  new_id=0
  t_eventend=0
  t_start=0
  event_count=0
  event_rec=0
  for(i in 1:n){
    new_id<-c(new_id,newid[i,])
    t_eventend<-c(t_eventend,t_someend[i,])
    t_start<-c(t_start,t_eventstart[i,])
    event_count<-c(event_count,event_num[i,])
    event_rec<-c(event_rec,event_recurrent[i,])
  }
  new_id<-new_id[new_id!=0]%>%unlist
  t_eventend<-t_eventend[t_eventend!=0]%>%unlist
  t_start<-t_start[!is.na(t_start)]%>%unlist
  t_start<-t_start[-1]
  event_rec<-event_rec[!is.na(event_rec)]%>%unlist
  event_rec<-event_rec[-1]
  event_count<-event_count[event_count!=0]%>%unlist
  data.frame(t_start,t_eventend)
  return(list(cou_exacer=cou_exacer,time_offset=time_offset,status=status,
              new_id=new_id,t_eventend=t_eventend,count_tte=count_tte,
              t_start=t_start,event_count=event_count,event_rec=event_rec))
}

#pwp
gen_new_nc_pwp<-function(n,beta,r){
  age<-round(runif(n=n,min=2,max=12))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  #error<-rnorm(n,mean=0,sd=0.001)
  error<-rgamma(n,1,300)
  #error<-runif(n,-0.002,0.002)
  #error<-rgamma(100,1,100)
  rate_run<-exp(beta[1]+beta[2]*age+beta[3]*gender+beta[4]*treatmentA)+error
  
  total_time<-rep(180,n)
  off_cou<-get_new_nc_pwp(total_time,rate_run,error)
  count_ex<-off_cou$cou_exacer
  #offset_run<-off_cou$time_offset
  count_ttevent<-off_cou$count_tte
  events<-off_cou$event_count
  t_to_start<-off_cou$t_start
  t_to_end<-off_cou$t_eventend
  id_new<-off_cou$new_id
  event_reccur<-off_cou$event_rec
  status<-off_cou$status
  ag_covariates <-gen_covariate(age,gender,treatmentA,count_ttevent)
  age_ag <-ag_covariates$age1
  gender_ag <-ag_covariates$gender1
  treatmentA_ag <-ag_covariates$treatment
  
  return(list(age=age,treatmentA=treatmentA,gender=gender,rate_run=rate_run,
              total_time=total_time,count_ex=count_ex,events=events,
              age_ag=age_ag,gender_ag=gender_ag,treatmentA_ag=treatmentA_ag,
              t_to_start=t_to_start,t_to_end=t_to_end,id_new=id_new,
              event_reccur=event_reccur,status=status,count_ttevent=count_ttevent))
}

test_one_offset_nc_pwp<-function(n,beta,r){
  dt<-NULL
  # dt<-tryCatch({gen_new(n,beta,r)},
  #               error = function(e) {
  #               NA
  #             })
  dt<-gen_new_nc_pwp(n,beta,r)
  #dt<-gen_Data_ag(n,beta,r)
  # if (!is.na(dt)){
  data<-data.frame(id=dt$id_new,time_start=dt$t_to_start,
                   time_end=dt$t_to_end,event=dt$event_reccur,
                   eventnums=dt$events,
                   age=dt$age_ag,gender=dt$gender_ag,
                   A=dt$treatmentA_ag)
  data1=data[!duplicated(data[,1]),]
  data1=data.frame(data1,status=dt$status)
  data_count<-data.frame(y=dt$count_ex,OS=dt$total_time,Age=dt$age,
                         Gender=dt$gender,A=dt$treatmentA)
  ag<-ag_one(data)
  cox<-cox_one(data1)
  pwp<-pwp_one(data)
  
  #poi<-poi.one(data_count)
  poi_off<-poi.offset(data_count)
  #over<-over.one(data_count)
  over_off<-over.offset(data_count)
  #nb<-nb.one(data_count)
  nb_off<-nb.offset(data_count)
  # } else {
  #    ag<-NA
  #    cox<-NA
  #    pwp<-NA
  #    poi_off<-NA
  #    over_off<-NA
  #    nb_off<-NA
  #  }     
  
  return(list(ag=ag,pwp=pwp,cox=cox,poi_off=poi_off,
              over_off=over_off,nb_off=nb_off))
}

get_power_offset_nc_pwp<-function(nrep,n,beta,r,threshold,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_one_offset_nc_pwp(n,beta,r))
  pow_ag<-mean(test1[1,]<level,na.rm=T)
  pow_pwp<-mean(test1[2,]<level,na.rm=T)
  pow_cox<-mean(test1[3,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[4,]<threshold,na.rm=T)
  #pow_over<-mean(test1[5,]<level,na.rm=T)
  pow_over_offset<-mean(test1[5,]<level,na.rm=T)
  #pow_nb<-mean(test1[7,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[6,]<level,na.rm=T)
  nanum<-sum(is.na(test1))
  return(list(pow_ag=pow_ag,pow_pwp=pow_pwp,
              pow_cox=pow_cox,pow_poi_offset=pow_poi_offset,
              pow_over_offset=pow_over_offset,
              pow_nb_offset=pow_nb_offset,nanum=nanum,test_poi=test1[4,]))
}


get_new_nc_pwp<-function(simple_time,rate,error){
  cou_exacer<-rep(0,n)
  time_offset<-matrix(0,nrow=n,ncol=40)
  status<-rep(1,n)
  for(i in 1:n){
    run_pp<-NA
    run_pp<-new_pp3(simple_time[i],rate[i],error)
    cou_exacer[i]<-unlist(run_pp[1])
    
    if (cou_exacer[i]==0){
      status[i]<-0
    }
    
    len_time<-unlist(run_pp[2])
    for (j in 1:length(len_time)){
      time_offset[i,j]<-len_time[j]
    }
  }
  
  for (i in 1:n){
    if (max(time_offset[i,])!=simple_time[i]){
      time_offset[i,which.max(time_offset[i,])+1]<-simple_time[i]
    }
  }
  
  #event_re<-rep(1,length(cou_exacer))
  #for (i in 1: length(event_re)){
  #  if (count_ex[i]==1||count_ex[i]==0){
  #    event_re[i]<-0
  #  }
  #}
  count_tte<-rowSums(time_offset!=0)
  newid<-matrix(0,nrow=n,ncol=40)
  t_eventstart<-matrix(0,nrow=n,ncol=40)
  t_someend<-matrix(0,nrow=n,ncol=40)
  event_num<-matrix(0,nrow=n,ncol=40)
  event_recurrent<-matrix(0,nrow=n,ncol=40)
  for (i in 1:n){
    # if (cou_exacer[i]>1){
    #    newid[i,1:cou_exacer[i]]<-i
    #    for (j in 1:cou_exacer[i]){
    #      event_num[i,j]<-j
    #    }
    if (count_tte[i]>1){
      newid[i,1:count_tte[i]]<-i
      for (j in 1:count_tte[i]){
        event_num[i,j]<-j
      }
    } else {
      newid[i,1]<-i
      event_num[i,1]<-1
    }
  }
  
  for(i in 1:n){
    if (count_tte[i]<=1){
      event_recurrent[i,1]<-0
      event_recurrent[i,2:40]<-NA
    } else {
      event_recurrent[i,1:count_tte[i]]<-1
      event_recurrent[i,(count_tte[i]+1):40]<-NA
    }
  }
  for (i in 1:n){
    #for (j in 1:40){
    if (count_tte[i]>2){
      for (j in 2:count_tte[i]){
        t_eventstart[i,j]<-time_offset[i,j-1]
      }
      t_eventstart[i,(count_tte[i]+1):40]<-NA
    } else if (count_tte[i]==2){
      t_eventstart[i,2]<-time_offset[i,1]
      t_eventstart[i,3:40]<-NA
    } else {
      t_eventstart[i,2:40]<-NA
    }
  }
  for (i in 1:n){
    #for (j in 1:40){
    if (count_tte[i]>=2){
      for (j in 1:count_tte[i]){
        t_someend[i,j]<-time_offset[i,j]
      }
    }  else {
      t_someend[i,1]<-time_offset[i,1]
    }
  }
  new_id=0
  t_eventend=0
  t_start=0
  event_count=0
  event_rec=0
  for(i in 1:n){
    new_id<-c(new_id,newid[i,])
    t_eventend<-c(t_eventend,t_someend[i,])
    t_start<-c(t_start,t_eventstart[i,])
    event_count<-c(event_count,event_num[i,])
    event_rec<-c(event_rec,event_recurrent[i,])
  }
  new_id<-new_id[new_id!=0]%>%unlist
  t_eventend<-t_eventend[t_eventend!=0]%>%unlist
  t_start<-t_start[!is.na(t_start)]%>%unlist
  t_start<-t_start[-1]
  event_rec<-event_rec[!is.na(event_rec)]%>%unlist
  event_rec<-event_rec[-1]
  event_count<-event_count[event_count!=0]%>%unlist
  data.frame(t_start,t_eventend)
  return(list(cou_exacer=cou_exacer,time_offset=time_offset,status=status,
              new_id=new_id,t_eventend=t_eventend,count_tte=count_tte,
              t_start=t_start,event_count=event_count,event_rec=event_rec))
}

run_type1_poi<-function(nrep,type1,alpha=0.05){
  type1<-sort(type1)
  threshold<-type1[nrep*alpha+1]
  return(threshold=threshold)
}

#testing_for inflation
gen_new_nc_test<-function(n,beta,r){
  age<-round(runif(n=n,min=2,max=12))
  treatmentA<-round(runif(n=n,min=0,max=1))
  gender<-round(runif(n=n,min=0,max=1))
  #error<-rnorm(n,mean=0,sd=0.001)
  #error<-rgamma(n,1,300)
  #error<-runif(n,-0.002,0.002)
  rate_run<-r*exp(beta[1]+beta[2]*age+beta[3]*gender+beta[4]*treatmentA)#+error
  
  total_time<-rep(180,n)
  off_cou<-get_new_nc_test(total_time,rate_run)
  count_ex<-off_cou$cou_exacer
  #offset_run<-off_cou$time_offset
  count_ttevent<-off_cou$count_tte
  events<-off_cou$event_count
  t_to_start<-off_cou$t_start
  t_to_end<-off_cou$t_eventend
  id_new<-off_cou$new_id
  event_reccur<-off_cou$event_rec
  status<-off_cou$status
  ag_covariates <-gen_covariate(age,gender,treatmentA,count_ttevent)
  age_ag <-ag_covariates$age1
  gender_ag <-ag_covariates$gender1
  treatmentA_ag <-ag_covariates$treatment
  
  return(list(age=age,treatmentA=treatmentA,gender=gender,rate_run=rate_run,
              total_time=total_time,count_ex=count_ex,events=events,
              age_ag=age_ag,gender_ag=gender_ag,treatmentA_ag=treatmentA_ag,
              t_to_start=t_to_start,t_to_end=t_to_end,id_new=id_new,
              event_reccur=event_reccur,status=status,count_ttevent=count_ttevent))
}

test_one_offset_nc_test<-function(n,beta,r){
  dt<-NULL
  # dt<-tryCatch({gen_new(n,beta,r)},
  #               error = function(e) {
  #               NA
  #             })
  dt<-gen_new_nc_test(n,beta,r)
  #dt<-gen_Data_ag(n,beta,r)
  # if (!is.na(dt)){
  data<-data.frame(id=dt$id_new,time_start=dt$t_to_start,
                   time_end=dt$t_to_end,event=dt$event_reccur,
                   eventnums=dt$events,
                   age=dt$age_ag,gender=dt$gender_ag,
                   A=dt$treatmentA_ag)
  data1=data[!duplicated(data[,1]),]
  data1=data.frame(data1,status=dt$status)
  data_count<-data.frame(y=dt$count_ex,OS=dt$total_time,Age=dt$age,
                         Gender=dt$gender,A=dt$treatmentA)
  ag<-ag_one(data)
  cox<-cox_one(data1)
  pwp<-pwp_one(data)
  
  #poi<-poi.one(data_count)
  poi_off<-poi.offset(data_count)
  #over<-over.one(data_count)
  over_off<-over.offset(data_count)
  #nb<-nb.one(data_count)
  nb_off<-nb.offset(data_count)
  # } else {
  #    ag<-NA
  #    cox<-NA
  #    pwp<-NA
  #    poi_off<-NA
  #    over_off<-NA
  #    nb_off<-NA
  #  }     
  
  return(list(ag=ag,pwp=pwp,cox=cox,poi_off=poi_off,
              over_off=over_off,nb_off=nb_off))
}

get_power_offset_nc_test<-function(nrep,n,beta,r,thres=0.05,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_one_offset_nc_test(n,beta,r))
  pow_ag<-mean(test1[1,]<level,na.rm=T)
  pow_pwp<-mean(test1[2,]<level,na.rm=T)
  pow_cox<-mean(test1[3,]<level,na.rm=T)
  pow_poi_offset_adjusted<-mean(test1[4,]<thres,na.rm=T)
  #pow_over<-mean(test1[5,]<level,na.rm=T)
  pow_over_offset<-mean(test1[5,]<level,na.rm=T)
  #pow_nb<-mean(test1[7,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[6,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[4,]<level,na.rm=T)
  nanum<-sum(is.na(test1))
  return(list(pow_ag=pow_ag,pow_pwp=pow_pwp,
              pow_cox=pow_cox,pow_poi_offset=pow_poi_offset,
              pow_over_offset=pow_over_offset,
              pow_nb_offset=pow_nb_offset,
              pow_poi_offset_adjusted=pow_poi_offset_adjusted,
              nanum=nanum,
              test_poi=test1[4,]))
}


get_new_nc_test<-function(simple_time,rate){
  cou_exacer<-rep(0,n)
  time_offset<-matrix(0,nrow=n,ncol=40)
  status<-rep(1,n)
  for(i in 1:n){
    run_pp<-NA
    run_pp<-new_pp4(simple_time[i],rate[i])
    cou_exacer[i]<-unlist(run_pp[1])
    
    if (cou_exacer[i]==0){
      status[i]<-0
    }
    
    len_time<-unlist(run_pp[2])
    for (j in 1:length(len_time)){
      time_offset[i,j]<-len_time[j]
    }
  }
  
  for (i in 1:n){
    if (max(time_offset[i,])!=simple_time[i]){
      time_offset[i,which.max(time_offset[i,])+1]<-simple_time[i]
    }
  }
  
  #event_re<-rep(1,length(cou_exacer))
  #for (i in 1: length(event_re)){
  #  if (count_ex[i]==1||count_ex[i]==0){
  #    event_re[i]<-0
  #  }
  #}
  count_tte<-rowSums(time_offset!=0)
  newid<-matrix(0,nrow=n,ncol=40)
  t_eventstart<-matrix(0,nrow=n,ncol=40)
  t_someend<-matrix(0,nrow=n,ncol=40)
  event_num<-matrix(0,nrow=n,ncol=40)
  event_recurrent<-matrix(0,nrow=n,ncol=40)
  for (i in 1:n){
    # if (cou_exacer[i]>1){
    #    newid[i,1:cou_exacer[i]]<-i
    #    for (j in 1:cou_exacer[i]){
    #      event_num[i,j]<-j
    #    }
    if (count_tte[i]>1){
      newid[i,1:count_tte[i]]<-i
      for (j in 1:count_tte[i]){
        event_num[i,j]<-j
      }
    } else {
      newid[i,1]<-i
      event_num[i,1]<-1
    }
  }
  
  for(i in 1:n){
    if (count_tte[i]<=1){
      event_recurrent[i,1]<-0
      event_recurrent[i,2:40]<-NA
    } else {
      event_recurrent[i,1:count_tte[i]]<-1
      event_recurrent[i,(count_tte[i]+1):40]<-NA
    }
  }
  for (i in 1:n){
    #for (j in 1:40){
    if (count_tte[i]>2){
      for (j in 2:count_tte[i]){
        t_eventstart[i,j]<-time_offset[i,j-1]
      }
      t_eventstart[i,(count_tte[i]+1):40]<-NA
    } else if (count_tte[i]==2){
      t_eventstart[i,2]<-time_offset[i,1]
      t_eventstart[i,3:40]<-NA
    } else {
      t_eventstart[i,2:40]<-NA
    }
  }
  for (i in 1:n){
    #for (j in 1:40){
    if (count_tte[i]>=2){
      for (j in 1:count_tte[i]){
        t_someend[i,j]<-time_offset[i,j]
      }
    }  else {
      t_someend[i,1]<-time_offset[i,1]
    }
  }
  new_id=0
  t_eventend=0
  t_start=0
  event_count=0
  event_rec=0
  for(i in 1:n){
    new_id<-c(new_id,newid[i,])
    t_eventend<-c(t_eventend,t_someend[i,])
    t_start<-c(t_start,t_eventstart[i,])
    event_count<-c(event_count,event_num[i,])
    event_rec<-c(event_rec,event_recurrent[i,])
  }
  new_id<-new_id[new_id!=0]%>%unlist
  t_eventend<-t_eventend[t_eventend!=0]%>%unlist
  t_start<-t_start[!is.na(t_start)]%>%unlist
  t_start<-t_start[-1]
  event_rec<-event_rec[!is.na(event_rec)]%>%unlist
  event_rec<-event_rec[-1]
  event_count<-event_count[event_count!=0]%>%unlist
  data.frame(t_start,t_eventend)
  return(list(cou_exacer=cou_exacer,time_offset=time_offset,status=status,
              new_id=new_id,t_eventend=t_eventend,count_tte=count_tte,
              t_start=t_start,event_count=event_count,event_rec=event_rec))
}
