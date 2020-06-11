setwd("/Users/xydai/Desktop/Thesis/Thesis_data")
study1<-read.csv("AZ4_03132020.csv")
library(lme4)
library(poisson)
library(MASS)
library(AER)
library(survival)
library(survMisc)
library(dplyr)
library(lattice)
library(survminer)
library(ggplot2)
#Poisson family
study1$pex_yn[is.na(study1$pex_yn)]<-0
study1$pex1_days[is.na(study1$pex1_days)]<-0
study1
study1$trt<-as.numeric(study1$trt)
study1$trt[study1$trt==2]<-0
select_study1<-study1[!duplicated(study1[,1]),]
study1_id<-select_study1$tdnccid
length(study1_id)
study1_num<-select_study1$pex_num
poi_study1<-select_study1[1,]

#for(i in 1:nrow(select_study1)){
#  test1<-NULL
#  order1<-NULL
#  test1<-filter(study1,tdnccid==study1_id[i])
#  test1<-test1[!duplicated(test1[,13]),]
#  order1<-test1[order(test1$fu_days),]
#  if (study1_num[i]>0){
#    poi_study1<-rbind(poi_study1,order1[study1_num[i],])
#  } else {
#    poi_study1<-rbind(poi_study1,order1[1,])
#  }
#}

for(i in 1:nrow(select_study1)){
  test1<-NULL
  order1<-NULL
  test1<-filter(study1,tdnccid==study1_id[i])
  len_test1<-nrow(test1)
  order1<-test1[order(test1$pex_days),]
  order1<-test1[order(test1$fu_days),]
  #if (study1_num[i]>0){
  #  poi_study1<-rbind(poi_study1,order1[study1_num[i],])
  #} else {
  #  poi_study1<-rbind(poi_study1,order1[1,])
  #}
  add<-order1[len_test1,]
  if(add$fu_days==0 & add$pex_days>add$fu_days){
    add$fu_days<-add$pex_days
  }
  poi_study1<-rbind(poi_study1,add)
}

poi_study1<-poi_study1[-1,]
nrow(poi_study1)
nrow(select_study1)

poi_study1$trt<-as.numeric(poi_study1$trt)
poi_study1$trt[poi_study1$trt==2]<-0
poi_study1$sex<-poi_study1$sex-1
poi_study1$fu_days

#Poisson regression with offset
mod_poi<-glm(pex_num~as.factor(age_group)+as.factor(sex)+trt, offset=log(fu_days),data=poi_study1,family=poisson)
summary(mod_poi)
pchisq(mod_poi$deviance, df=mod_poi$df.residual, lower.tail=FALSE)
#p-value=0.65 fit pretty well

mod_glmm<-glmer(pex_num~age_group+sex+trt+(1|tdnccid),offset=log(fu_days),family=poisson,data=poi_study1)
summary(mod_glmm)
ran_az<-ranef(mod_glmm,condVar = TRUE)
ran_a<-ran_az$tdnccid
ran_a<-ran_a[,1]
random_az=data.frame(random_az=exp(ran_a))
az<-exp(ran_a)
ggplot(random_az,aes(x=random_az))+
  geom_histogram(fill="ivory4",color="ivory4",binwidth = 0.1,aes(y=..density..))+
  #geom_density()+
  stat_function(fun=dlnorm, args=list(mean=0.02736077, sd=0.26341334),color="mediumpurple4")+
  stat_function(fun=dgamma, args=list(shape=12.57665, rate=11.75399),color="salmon4")+
  labs(x = "random effect")+
  scale_colour_manual("distribution", values = c("gaussian"="mediumpurple4","gamma"="salmon4"))+
  ggtitle("Histogram for random effect")
fit.gamma <- fitdist(az, "gamma", lower = c(0, 0))
fit.norm<-fitdist(ran_a,"norm")
fit.norm
fit.gamma
mod_quasi<-glm(pex_num~age_group+sex+trt, offset=log(fu_days),data=poi_study1,family=quasipoisson)
summary(mod_quasi)
pchisq(mod_quasi$deviance, df=mod_quasi$df.residual, lower.tail=FALSE)
#p-value=0.65 fit pretty well
#overdispersion 1.1, basically no overdispersion

mod_nb<-glm.nb(pex_num~age_group+sex+trt+offset(log(fu_days)),data=poi_study1)
summary(mod_nb)
pchisq(mod_nb$deviance, df=mod_nb$df.residual, lower.tail=FALSE)
#p-value=0.96 fit pretty well

hist_study1<-hist(poi_study1$pex_num)
test_disperse1<-dispersiontest(mod_poi)
test_disperse1

#cox
poi_study1
cox_study1<-poi_study1
#cox_study1<-cox_study1[cox_study1$pex_days!=0,]
cox_study1
nrow(cox_study1)

for (i in 1:nrow(cox_study1)){
  if (cox_study1$pex_days[i]!=cox_study1$fu_days[i]){
    cox_study1$pex_days[i]<-cox_study1$fu_days[i]
  }
}

for (i in 1:nrow(cox_study1)){
  if (cox_study1$pex_days[i]!=cox_study1$pex1_days[i] && cox_study1$pex1_days[i]!=0){
    cox_study1$pex_days[i]<-cox_study1$pex1_days[i]
  }
}

cox_study1$status<-as.numeric(cox_study1$pex_num>0)
cox_study1$status
cox_study1
surv_study1<-Surv(cox_study1$pex_days, cox_study1$status)
surv_study1
mod_cox<-coxph(surv_study1 ~trt+sex+strata(age_group),
                    data=cox_study1)
mod_cox
plot(az4.ft,xlab="Observed follow-up time (days)",ylab="Survival") 
title("Kaplan-Meier Survival Estimate (Overall)")
az4.ft<-survfit(surv_study1~1,conf.type="log-log",conf.int=0.95)

cox_study_f_1<-filter(cox_study1,sex==1,age_group==1)
cox_study_f_2<-filter(cox_study1,sex==1,age_group==2)

cox_study_m_1<-filter(cox_study1,sex==0,age_group==1)
cox_study_m_2<-filter(cox_study1,sex==0,age_group==2)

surv_studyf1<-with(cox_study_f_1,Surv(pex_days, status))
surv_studyf2<-with(cox_study_f_2,Surv(pex_days, status))
surv_studym1<-with(cox_study_m_1,Surv(pex_days, status))
surv_studym2<-with(cox_study_m_2,Surv(pex_days, status))

library(ggplot2)
plot1.1<-ggsurvplot(trt_fit1, #conf.int = TRUE,
                    #conf.int.style = "step",
                    data = cox_study_f_1,
                    legend.title = "Treatment Arm",
                    legend.labs = 
                      c("Placebo", "Azithromycin"),
                    ylab="PEx free probability",
                    xlab = "Time in days")
plot1.1+ggtitle("Kaplan-Meier Estimate (6-12 yr, Female)")
trt_fit1 <- survfit(surv_studyf1 ~ trt, data=cox_study_f_1,conf.type="log-log",conf.int=0.95)
plot(trt_fit1,xlab="Observed follow-up time (days)",ylab="Survival") 
title("Kaplan-Meier Survival Estimate (trt)")
gof(mod_cox)
mod_coxm1<-coxph(surv_studym2 ~trt,
               data=cox_study_m_2)
test.ph<-cox.zph(mod_coxm1)
test.ph

#pwp and AG
study1$time_start<-0
study1$event_pwp<-1
study1$event_ag<-0
ag_study1<-study1[1,]
for(i in 1:length(study1_id)){
  test1<-NULL
  order1<-NULL
  
  test1<-filter(study1,tdnccid==study1_id[i])
  len_test1<-nrow(test1)
  maxi_test1<-max(test1$pex_num)
  maxi_fu1<-max(test1$fu_days)
  
  if(maxi_test1==0){
    order1<-test1[order(test1$pex_days),]
    order1<-test1[order(test1$fu_days),]
    order1<-order1[len_test1,]
    order1$pex_days<-order1$fu_days
  } else{
    #pex_date<-max(order1$pex1_days)
    largest_days<-test1[test1$fu_days==maxi_fu1,][1,]
    largest_days$pex_days<-largest_days$fu_days
    order1<-test1[1:maxi_test1,]
    max_pex=max(order1$pex_days)
    if (largest_days$pex_days!=max_pex){
      order1=rbind(order1,largest_days)
    }
    if(nrow(order1)>1){
      for (j in 2:nrow(order1)){
        #if (order1[j,9] < pex_date){
        order1[j,15]<-order1[j-1,9]
        order1[j,16]<-j
        #} else {
        #  order1[j,13]<-pex_date
        #  order1[j,15]<-order1[j-1,13]
      }
    }
    
    if (maxi_test1>1){
      order1[,17]<-1
      if(!identical(order1[maxi_test1+1,13],max_pex)){
        order1[maxi_test1+1,17]<-0
      }
    } else if (maxi_test1==1){
      order1[1,17]<-1
      if(!identical(order1[2,13],max_pex)){
        order1[2,17]<-0
      }
    }
    
    }
    
  #}
  ag_study1<-rbind(ag_study1,order1)
}

ag_study1<-ag_study1[-1,]
ag_study1
ag_study1$trt<-as.numeric(ag_study1$trt)
ag_study1$trt[ag_study1$trt==2]<-0
mod_ag<-coxph(Surv(pwp_0, pex_days, event_ag) ~ 
               strata(age_group) +sex+trt +cluster(tdnccid),
              data=ag_study1)
mod_ag
ag_study1$pwp_0<-0
summary(mod_ag)
mod_pwp<-coxph(Surv(time_start,pex_days,event_ag)~
                 strata(age_group)+sex+trt+cluster(tdnccid)+strata(event_pwp), data=ag_study1)
mod_pwp
mod_pwp<-coxph(Surv(pwp_0,pex_days,event_ag)~
                 strata(age_group)+sex+trt+cluster(tdnccid)+strata(event_pwp), data=ag_study1)
mod_pwp
summary(mod_pwp)
which(ag_study1$pex_days<=ag_study1$time_start)
ag_study1[40,]
ag_study1[75,]
ag_study1[302,]

length(which(poi_study1$pex_num>3))

library(tableone)
vars=c("age_rand","age_group","sex","pex_num","fu_days")
factorvars=c("age_group","sex","pex_num")
tableOneAZE<-CreateTableOne(vars=vars,strata="trt",
                            data=poi_study1,factorVars = factorvars)
print(tableOneAZE)
str(poi_study1)

#second and gap
cox2_study1<-study1[1,]
for(i in 1:nrow(select_study1)){
  test2<-NULL
  order2<-NULL
  test2<-filter(study1,tdnccid==study1_id[i])
  #len_test2<-nrow(test2)
  num_pex<-test2$pex_num[1]
  order2<-test2[order(-test2$fu_days),]
  order2<-order2[order(order2$pex_days),]
  order2<-order2[!duplicated(order2$pex_days),]
  #if (study2_num[i]>0){
  #  poi_study2<-rbind(poi_study2,order2[study2_num[i],])
  #} else {
  #  poi_study2<-rbind(poi_study2,order2[1,])
  #}
  if(num_pex>=2){
    event_some<-2
  } else {
    event_some<-1
    if (order2[event_some,]$fu_days>order2[event_some,]$pex_days){
      order2[event_some,]$pex_days<-order2[event_some,]$fu_days
    }
    
  }
  cox2_study1<-rbind(cox2_study1,order2[event_some,])
}
cox2_study1<-cox2_study1[-1,]
new_cox2<-cox2_study1[cox2_study1$pex_num>0,]
nrow(new_cox2)
nrow(gap_study1)
new_cox2$status<-as.numeric(new_cox2$pex_num>=2)
new_cox2$trt[new_cox2$trt==2]<-0
which(cox2_study1$pex_days==0)
a<-as.numeric(new_cox2$trt)
a
new_cox2
surv1_second<-Surv(new_cox2$pex_days, new_cox2$status)
surv1_second
new_cox2
length(which(cox2_study1$status==1))
nrow(cox2_study1)
trt_fit1 <- survfit(surv_gap1 ~ trt, data=gao_study1,conf.type="log-log",conf.int=0.95)
plot1.1<-ggsurvplot(trt_fit1, #conf.int = TRUE,
                    #conf.int.style = "step",
                    data = gao_study1,
                    legend.title="Treatment Arm",
                    legend.labs = 
                      c("Placebo","Azithromycin"),
                    ylab=">=2nd PE free probability",
                    xlab = "Time in days")
plot1.1+ggtitle("Kaplan-Meier Estimate (second event since first)")
nrow(gap_study2)
#plot(opt.ft,xlab="Observed follow-up time (days)",ylab="Survival") 
#title("Kaplan-Meier Survival Estimate (0-second)")
#opt.ft<-survfit(surv_second~1,conf.type="log-log",conf.int=0.95)
##model cox
mod_second1<-coxph(surv1_second~trt+strata(age_group)+sex,
                  data=new_cox2)
summary(mod_second1)
new_cox2
#gap_cox_first to second

study1$start_first<-0
gap_study1<-study1[1,]
for(i in 1:nrow(select_study1)){
  test2<-NULL
  order2<-NULL
  test2<-filter(study1,tdnccid==study1_id[i])
  #len_test2<-nrow(test2)
  num_pex<-test2$pex_num[1]
  order2<-test2[order(-test2$fu_days),]
  order2<-order2[order(order2$pex_days),]
  order2<-order2[!duplicated(order2$pex_days),]
  #if (study2_num[i]>0){
  #  poi_study2<-rbind(poi_study2,order2[study2_num[i],])
  #} else {
  #  poi_study2<-rbind(poi_study2,order2[1,])
  #}
  if(num_pex>=2){
    event_some<-2
    needed<-order2[event_some,]
    needed$start_first<-order2[event_some-1,]$pex_days
  } else if (num_pex==1){
    event_some<-1
    needed<-order2[event_some,]
    needed$start_first<-order2[event_some,]$pex_days
    needed$pex_days<-order2[event_some,]$fu_days
  } else {
    event_some<-1
    needed<-order2[event_some,]
  }
  
  gap_study1<-rbind(gap_study1,needed)
}
gap_study1$trt<-as.numeric(gap_study1$trt)
gap_study1$trt[gap_study1$trt==2]<-0
gap_study1$trt
gap_study1<-gap_study1[-1,]
nrow(gap_study1)
which(gap_study2$pex_days==0)
gap_study2[203,]
gap_study1$status<-as.numeric(gap_study1$pex_num>=2)
gap_study1
which(gap_study2$pex_days<gap_study2$start_first)
gap_study2[37,]
#gap_study1<-gap_study1[(gap_study1$pex1_days>0|gap_study1$pex_num>0),]
surv_gap1<-with(gao_study1,Surv(start_first,pex_days,status))
surv_gap1

kmfit<- survfit(surv_gap1~1,data=gao_study1,type="kaplan-meier")
kmfit
a<-gap_study2[gap_study2$start_first>=13,]
a$start_first
min(a$start_first)
a$pex_num
#plot(opt.ft,xlab="Observed follow-up time (days)",ylab="Survival") 
#title("Kaplan-Meier Survival Estimate (gap)")
#opt.ft<-survfit(surv_gap~1,conf.type="log-log",conf.int=0.95)
##model cox
mod_gap<-coxph(surv_gap1~trt+strata(age_group)+sex,
               data=gap_study1)
summary(mod_gap)
gao_study1<-gap_study1[cox2_study1$pex_num>0,]
nrow(gao_study1)
#plot(opt.ft,xlab="Observed follow-up time (days)",ylab="Survival") 
#title("Kaplan-Meier Survival Estimate (gap)")
#opt.ft<-survfit(surv_gap~1,conf.type="log-log",conf.int=0.95)
##model cox
mod_gap1<-coxph(surv_gap1~trt+strata(age_group)+sex,
               data=gao_study1)
mod_gap1
summary(mod_gap1)

leftcen1<-gap_study1$pex1_days
leftcen1<-leftcen1[!is.na(leftcen1)]
leftcen1<-leftcen1[leftcen1>0]
mincen1<-min(leftcen1)
length(leftcen1)
mincen1
leftcen1

for (i in (1:nrow(gap_study1))){
  if (gap_study1[i,]$pex_num<1){
    gap_study1[i,]$start_first<-0
  }
}

gap_study1

#cox3
cox3_study1<-study1[1,]
for(i in 1:nrow(select_study1)){
  test2<-NULL
  order2<-NULL
  test2<-filter(study1,tdnccid==study1_id[i])
  #len_test2<-nrow(test2)
  num_pex<-test2$pex_num[1]
  order2<-test2[order(-test2$fu_days),]
  order2<-order2[order(order2$pex_days),]
  order2<-order2[!duplicated(order2$pex_days),]
  #if (study2_num[i]>0){
  #  poi_study2<-rbind(poi_study2,order2[study2_num[i],])
  #} else {
  #  poi_study2<-rbind(poi_study2,order2[1,])
  #}
  if(num_pex>=3){
    event_some<-3
  } else if (num_pex==2){
    event_some<-2
    if (order2[event_some,]$fu_days>order2[event_some,]$pex_days){
      order2[event_some,]$pex_days<-order2[event_some,]$fu_days
    }
  } else {
    event_some<-1
  }
    cox3_study1<-rbind(cox3_study1,order2[event_some,])

}
cox3_study1<-cox3_study1[-1,]
new_cox3<-cox3_study1[cox3_study1$pex_num>1,]
nrow(new_cox3)
nrow(gap_study1)
new_cox3$status<-as.numeric(new_cox3$pex_num>2)
new_cox2$trt[new_cox2$trt==2]<-0
which(cox2_study1$pex_days==0)
a<-as.numeric(new_cox2$trt)
a
new_cox2
surv1_3rd<-Surv(new_cox3$pex_days, new_cox3$status)
surv1_3rd
new_cox2
length(which(cox2_study1$status==1))
nrow(cox2_study1)
trt_fit1 <- survfit(surv1_3rd ~ trt, data=new_cox3,conf.type="log-log",conf.int=0.95)
plot1.1<-ggsurvplot(trt_fit1, #conf.int = TRUE,
                    #conf.int.style = "step",
                    data = new_cox3,
                    legend.title="Treatment Arm",
                    legend.labs = 
                      c("Placebo","Azithromycin"),
                    ylab=">=2nd PE free probability",
                    xlab = "Time in days")
plot1.1+ggtitle("Kaplan-Meier Estimate (second event since first)")
nrow(gap_study2)
#plot(opt.ft,xlab="Observed follow-up time (days)",ylab="Survival") 
#title("Kaplan-Meier Survival Estimate (0-second)")
#opt.ft<-survfit(surv_second~1,conf.type="log-log",conf.int=0.95)
##model cox
mod_cox3<-coxph(surv1_3rd~trt+strata(age_group)+sex,
                   data=new_cox3)
summary(mod_cox3)
new_cox3
