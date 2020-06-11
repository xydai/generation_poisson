setwd("/home/students/daix4/Thesis/")
gen_data_offset<-function(n,beta,r,method="nc"){
  age<-round(runif(n=n,min=1,max=18))
  treatmentA<-rbinom(n=n,1,0.5)
  gender<-rbinom(n=n,1,0.5)
  
  ploteff=rnorm(n,0,0.6)
  error<-rgamma(n,1,1)
  age_group<-rep(0,n)
  age_group[age<=3]=1
  age_group[age>3 & age<=6]=2
  age_group[age>6 & age<=12]=3
  age_group[age>12]=4

  if (beta[5]==1){
    gam_err<-error
  } else if (beta[5]==0){
    gam_err<-1
  }
  
  if (beta[6]==1){
    ran<-read.csv("random_effect.csv")
    ran_1<-ran$x
    raneff<-gen_raneff(n,ran_1)
  } else if (beta[6]==0){
    raneff<-0
  }
  
  rate_run<-r*exp(beta[1]*age_group+beta[2]*gender+beta[3]*treatmentA+beta[4]*ploteff
                  +beta[6]*raneff)*gam_err
  
  if(method=="nc"){
    total_time<-rep(180,n)
    off_cou<-get_new_nc(total_time,rate_run,n)
  } else if (method=="mcar"){
    total_time<-get_total_time_mcar(n,0.001)
    off_cou<-get_new(total_time,rate_run,n)
  } else if (method=="mar"){
    total_time<-get_total_time_mar(n,rate_run,6)
    off_cou<-get_new(total_time,rate_run,n)
  }
  
  
  count_ex<-off_cou$cou_exacer
  #offset_run<-off_cou$time_offset
  count_ttevent<-off_cou$count_tte
  events<-off_cou$event_count
  t_to_start<-off_cou$t_start
  t_to_end<-off_cou$t_eventend
  id_new<-off_cou$new_id
  event_reccur<-off_cou$event_rec
  status<-off_cou$status
  ag_covariates <-gen_covariate(age_group,gender,treatmentA,count_ttevent)
  age_ag <-ag_covariates$age1
  gender_ag <-ag_covariates$gender1
  treatmentA_ag <-ag_covariates$treatment
  
  return(list(age=age_group,treatmentA=treatmentA,gender=gender,rate_run=rate_run,
              total_time=total_time,count_ex=count_ex,events=events,
              age_ag=age_ag,gender_ag=gender_ag,treatmentA_ag=treatmentA_ag,
              t_to_start=t_to_start,t_to_end=t_to_end,id_new=id_new,
              event_reccur=event_reccur,status=status,count_ttevent=count_ttevent,
              age_notgroupgroup=age))
}

test_offset<-function(n,beta,r,method="nc"){
  dt<-NULL
  dt<-tryCatch({gen_data_offset(n,beta,r,method)},
               error = function(e) {
                 NA
               })

  if (!is.na(dt)){
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
    poi_off<-poi.offset(data_count)
    over_off<-over.offset(data_count)
    nb_off<-nb.offset(data_count)
    glmm_off<-poi.glmm(data_count)
  } else {
    ag<-NA
    cox<-NA
    pwp<-NA
    poi_off<-NA
    over_off<-NA
    nb_off<-NA
    glmm_off<-NA
  }    
  
  return(list(ag=ag,pwp=pwp,cox=cox,poi_off=poi_off,
              over_off=over_off,nb_off=nb_off,glmm_off=glmm_off))
}

power_offset<-function(nrep,n,beta,r,method="nc",thres=0.05,level=0.05){
  test1<-NULL
  test1<- replicate(nrep, test_offset(n,beta,r,method))
  pow_ag<-mean(test1[1,]<level,na.rm=T)
  pow_pwp<-mean(test1[2,]<level,na.rm=T)
  pow_cox<-mean(test1[3,]<level,na.rm=T)
  pow_poi_offset<-mean(test1[4,]<level,na.rm=T)
  pow_poi_adjusted<-mean(test1[4,]<thres,na.rm=T)
  pow_over_offset<-mean(test1[5,]<level,na.rm=T)
  pow_nb_offset<-mean(test1[6,]<level,na.rm=T)
  pow_glmm<-mean(test1[7,]<level,na.rm=T)
  
  nanum<-sum(is.na(test1))
  return(list(pow_ag=pow_ag,pow_pwp=pow_pwp,
              pow_cox=pow_cox,pow_poi_offset=pow_poi_offset,
              pow_over_offset=pow_over_offset,
              pow_poi_adjusted=pow_poi_adjusted,
              pow_nb_offset=pow_nb_offset,pow_glmm=pow_glmm,
              nanum=nanum,test_poi=test1[4,]))
}

get_total_time_mcar<-function(n,rc){
  total_t<-rep(0,n)
  for (i in 1:n){
    total_t[i]<-min(180,rexp(1,rc))
  }
  return(total_t=total_t)
}

get_total_time_mar<-function(n,rate,rc){
  total_t<-rep(0,n)
  for (i in 1:n){
    total_t[i]<-min(180,rexp(1,rate[i]/rc))
  }
  return(total_t=total_t)
}

gen_raneff<-function(n,dis){
  ran_eff<-rep(0,n)
  for (i in 1:n){
    ran_eff[i]<-sample(dis,1)
  }
  return(ran_eff=ran_eff)
}
