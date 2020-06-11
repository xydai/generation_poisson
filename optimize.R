setwd("/Users/xydai/Desktop/Thesis")
study2<-read.csv("OPTIMIZE_03132020.csv")
library(lme4)
library(poisson)
library(MASS)
library(AER)
library(survival)
library(survMisc)
library(dplyr)
library(lattice)
library(fitdistrplus)
library(ggplot2)
#Poisson family
study2
study2$pex_yn[is.na(study2$pex_yn)]<-0
study2$pex1_days[is.na(study2$pex1_days)]<-0
study2
select_study2<-study2[!duplicated(study2[,1]),]
study2_id<-select_study2$tdnccid
study2_num<-select_study2$pex_num
poi_study2<-select_study2[1,]
length(study2_id)


for(i in 1:nrow(select_study2)){
  test2<-NULL
  order2<-NULL
  test2<-filter(study2,tdnccid==study2_id[i])
  len_test2<-nrow(test2)
  order2<-test2[order(test2$fu_days),]
  #if (study2_num[i]>0){
  #  poi_study2<-rbind(poi_study2,order2[study2_num[i],])
  #} else {
  #  poi_study2<-rbind(poi_study2,order2[1,])
  #}
  add<-order2[len_test2,]
  #if(add$fu_days==0 & add$pex_days>add$fu_days){
  #  add$fu_days<-add$pex_days
  #}
  poi_study2<-rbind(poi_study2,add)
}

poi_study2<-poi_study2[-1,]
nrow(poi_study2)
nrow(select_study2)

#poi_study2$trt<-as.numeric(poi_study2$trt)
poi_study2$trt[poi_study2$trt==2]<-0
poi_study2$sex<-poi_study2$sex-1
which(poi_study2$fu_days==0)
poi_study2[203,]
#one person do not have valid follow up days
poi_study2<-poi_study2[poi_study2$fu_days>0,]
nrow(poi_study2)
poi_study2$sex
#fit poi family models
mod_poi2<-glm(pex_num~as.factor(age_group)+trt, offset=log(fu_days),data=poi_study2,family=poisson)
summary(mod_poi2)
pchisq(mod_poi2$deviance, df=mod_poi2$df.residual, lower.tail=FALSE)
#p_value=0.0001, fit pretty bad
#sex as confounder (not significant), with interaction, significant


mod_glmm2<-glmer(pex_num~1+as.factor(age_group)+trt+(1|tdnccid), offset=log(fu_days), family=poisson,data=poi_study2)
summary(mod_glmm2)

#str(ran_glmm2<-ranef(mod_glmm2))
ran<-ranef(mod_glmm2,condVar = TRUE)
ran

#random_opt<-qqmath(ranef(mod_glmm2, condVar=TRUE)) 
#dotplot(ranef(mod_glmm2, condVar=TRUE))
#randoms<-ranef(mod_glmm2, condVar = TRUE)
#qq <- attr(ranef(mod_glmm2, condVar = TRUE)[[1]],"postVar")
#rand.interc<-randoms$tdnccid
#df<-data.frame(Intercepts=rand.interc[,1],
#               sd.interc=2*sqrt(qq[,,1:length(qq)]),
#               lev.names=rownames(rand.interc))
#df$lev.names<-factor(df$lev.names,levels=df$lev.names[order(df$Intercepts)])
#p <- ggplot(df,aes(lev.names,Intercepts,shape=lev.names))
#qqnorm(rand.interc[,1])
#qqline(rand.interc[,1])
ran1<-ran$tdnccid
ran1<-ran1[,1]
random_effect=data.frame(random_effect=exp(ran1))
random_effect<-random_effect$random_effect
effect<-exp(ran1)
ran1
write.csv(ran1,paste0("random_effect.csv"), row.names = T)

fit.gamma <- fitdist(effect, "gamma", lower = c(0, 0))
fit.norm<-fitdist(ran1,"norm")
fit.gamma
fit.norm

ploy_ran_2<-ggplot(random_effect,aes(x=random_effect))+
  geom_histogram(fill="ivory4",color="ivory4",binwidth = 0.1,aes(y=..density..))+
  stat_function(fun=dgamma, args=list(shape=5.416021, rate=4.682933),color="salmon4")+
  stat_function(fun=dlnorm, args=list(mean=0.05031582, sd=0.40193182),color="mediumpurple4")+
  #geom_line(aes(x=exp(ran1), y=dgamma(exp(ran1),fit.gamma$estimate["shape"], fit.gamma$estimate["rate"])), color="salmon4")+
  labs(x = "random effect")+
  #geom_density()+
  scale_colour_manual("distribution", values = c("gaussian"="mediumpurple4","gamma"="salmon4"))+
  ggtitle("Histogram for random effect")
  
f=data.frame(f)
plots<-ggplot(f,aes(x=f))+
  geom_histogram(fill="ivory4",color="ivory4",binwidth = 0.1,aes(y=..density..))+
  #stat_function(fun=dgamma, args=list(shape=5.182835, rate=2.013795),color="salmon4")+
  stat_function(fun=dnorm, args=list(mean=-0.04861307, sd=0.6),color="mediumpurple4")+
  #geom_line(aes(x=exp(ran1), y=dgamma(exp(ran1),fit.gamma$estimate["shape"], fit.gamma$estimate["rate"])), color="salmon4")+
  labs(x = "random effect")+
  scale_colour_manual("distribution", values = c("gaussian"="mediumpurple4","gamma"="salmon4"))+
  ggtitle("Histogram for random effect")
plots
ploy_ran_2
#fit.params <- fitdistr(exp(ran1), "gamma")
#fit.params
den<-density(effect)
summary(den)
shapiro.test(ran1)
mod_quasi2<-glm(pex_num~as.factor(age_group)+trt, offset=log(fu_days),data=poi_study2,family=quasipoisson)
summary(mod_quasi2)
pchisq(mod_quasi2$deviance, df=mod_quasi2$df.residual, lower.tail=FALSE)
#sex as confounder (not significant), with interaction, significant
#but it is for 0.1 not 0.05 (p=0.0646)
#p value for goodness of fit is 0.0001, fit pretty bad
#overdispersion 1.7
mod_nb2<-glm.nb(pex_num~as.factor(age_group)+trt+offset(log(fu_days)),data=poi_study2)
summary(mod_nb2)
pchisq(mod_nb2$deviance, df=mod_nb2$df.residual, lower.tail=FALSE)
#p value for goodness of fit is 0.602, fit pretty well
#sex as confounder (not significant), with interaction, significant
#overdispersion=1.48
PEx_for_OPTIMIZE=poi_study2$pex_num
hist_study2<-hist(PEx_for_OPTIMIZE)

test_disperse<-dispersiontest(mod_poi2)
test_disperse

#cox
cox_study2<-poi_study2
#cox_study2<-cox_study2[cox_study2$pex_days!=0,]
cox_study2
nrow(cox_study2)

for (i in 1:nrow(cox_study2)){
  if (cox_study2$pex_days[i]!=cox_study2$fu_days[i]&&cox_study2$pex_days[i]<=cox_study2$fu_days[i]){
    cox_study2$pex_days[i]<-cox_study2$fu_days[i]
  }
}

for (i in 1:nrow(cox_study2)){
  if (cox_study2$pex_days[i]!=cox_study2$pex1_days[i] && cox_study2$pex1_days[i]!=0){
    cox_study2$pex_days[i]<-cox_study2$pex1_days[i]
  }
}

cox_study2$status<-as.numeric(cox_study2$pex_num>0)
cox_study2$status
#cox_study2$month=cox_study2$pex_days/30
surv_study2<-Surv(cox_study2$pex_days, cox_study2$status)
surv_study2

##model cox
mod_cox2<-coxph(surv_study2 ~trt+strata(age_group),
               data=cox_study2)
summary(mod_cox2)

##by trt
summary(cox_study2)
cox_study2$trt[cox_study2$trt==2]<-0
cox_1<-filter(cox_study2,age_group==1)
cox_2<-filter(cox_study2,age_group==2)
cox_3<-filter(cox_study2,age_group==3)
cox_4<-filter(cox_study2,age_group==4)
surv_study_age4<-with(cox_1,Surv(pex_days, status))
surv_study_age1
cox_1

trt_fit <- survfit(surv_second ~ trt, data=cox2_study2,conf.type="log-log",conf.int=0.95)
plot1.1<-ggsurvplot(trt_fit, #conf.int = TRUE,
                    #conf.int.style = "step",
                    data = cox2_study2,
                    legend.title = "Treatment Arm",
                    legend.labs = 
                      c("Placebo","Azithromycin"),
                    ylab="PEx free probability",
                    xlab = "Time in days")
plot1.1+ggtitle("Kaplan-Meier Estimate (second event from time 0)")
plot(trt_fit,col=c("red","blue"),xlab="Observed follow-up time (days)",ylab="Survival") 
title("Kaplan-Meier Survival Estimate by Treatment(gap)")
par(mfrow=c(3,1))
legend("topright",c("Placebo","Azithromycin"),
       col=c("red","blue"),lwd=rep(3,4))
cox_study2

mod_cox_age<-coxph(surv_study_age3 ~trt,
                   data=cox_3)
gof(mod_cox2)
test.ph2<-cox.zph(mod_cox_age)
test.ph2


mod_cox_age4<-coxph(surv_study_age4 ~trt,
                data=cox_4)
mod_cox_age3
test1<-cox.zph(mod_cox_age4)
test1
plot(test1)
print(test1)

#by trt
study2_p<-cox_study2[cox_study2$trt==0,]
study2_a<-cox_study2[cox_study2$trt==1,]
study2_p$fu_mo<-study2_p$fu_days/30
summary(study2_p)
#pwp and AG
study2$time_start<-0
study2$event_pwp<-1
study2$event_ag<-0
#study2$pwp_0<-0
ag_study2<-study2[1,]
for(i in 1:length(study2_id)){
  test2<-NULL
  order2<-NULL
  
  test2<-filter(study2,tdnccid==study2_id[i])
  len_test2<-nrow(test2)
  maxi_test2<-max(test2$pex_num)
  maxi_fu2<-max(test2$fu_days)
  
  if(maxi_test2==0){
    order2<-test2[order(test2$pex_days),]
    order2<-test2[order(test2$fu_days),]
    order2<-order2[len_test2,]
    order2$pex_days<-order2$fu_days
  } else{
    #pex_date<-max(order2$pex1_days)
    largest_days<-test2[test2$fu_days==maxi_fu2,][1,]
    largest_days$pex_days<-largest_days$fu_days
    order2<-test2[1:maxi_test2,]
    max_pex=max(order2$pex_days)
    if (largest_days$pex_days!=max_pex){
      order2=rbind(order2,largest_days)
    }
    if(nrow(order2)>1){
      for (j in 2:nrow(order2)){
        #if (order2[j,9] < pex_date){
        order2[j,15]<-order2[j-1,9]
        order2[j,16]<-j
        #} else {
        #  order2[j,13]<-pex_date
        #  order2[j,15]<-order2[j-1,13]
      }
    }
    if (maxi_test2>1){
      order2[,17]<-1
      if(!identical(order2[maxi_test2+1,13],max_pex)){
        order2[maxi_test2+1,17]<-0
      }
    } else if (maxi_test2==1){
      order2[1,17]<-1
      if(!identical(order2[2,13],max_pex)){
        order2[2,17]<-0
      }
    }
  }
  
  
  ag_study2<-rbind(ag_study2,order2)
}

ag_study2<-ag_study2[-1,]
ag_study2<-ag_study2[!is.na(ag_study2$tdnccid),]
ag_study2$trt[ag_study2$trt==2]<-0

for (i in 1:nrow(ag_study2)){
  a<-ag_study2[i,]$pex_days
  b<-ag_study2[i,]$time_start
  if (a<b){
    x<-ag_study2[i,]$time_start
    ag_study2[i,]$time_start<-ag_study2[i,]$pex_days
    ag_study2[i,]$pex_days<-x
  }
}

which(ag_study2$pex_days==ag_study2$time_start)
ag_study2[57,]

mod_ag2<-coxph(Surv(time_start, pex_days, event_ag) ~ 
                strata(age_group) + trt +cluster(tdnccid),
              data=ag_study2)
summary(mod_ag2)

pwp_study2<-ag_study2
pwp_study2<-pwp_study2[pwp_study2$event_pwp<=3,]
#ag_study2<-ag_study2[pwp_study2$event_pwp<3,]
pwp_study2
nrow(ag_study2)
nrow(pwp_study2)
mod_pwp2<-coxph(Surv(time_start,pex_days,event_ag)~
                 strata(age_group)+trt+cluster(tdnccid)+strata(event_pwp), data=pwp_study2)
summary(mod_pwp2)
mod_pwp2<-coxph(Surv(pwp_0,pex_days,event_ag)~
                  strata(age_group)+trt+cluster(tdnccid)+strata(event_pwp), data=pwp_study2)
summary(mod_pwp2)
#mod_pwp_full<-coxph(Surv(time_start,pex_days,event_ag)~
#                      strata(age_group)+trt+cluster(tdnccid)+strata(event_pwp), data=ag_study2)
#mod_pwp_full
which(ag_study2$pex_days<=ag_study2$time_start)

#logistiv/linear
poi_study2$recurrent<-0
poi_study2$recurrent[poi_study2$pex_num>=1]<-1
mod_log<-glm(recurrent ~ trt+age_group, data = poi_study2, family = binomial)
mod_log
summary(mod_log)
pchisq(mod_log$deviance, df=mod_log$df.residual, lower.tail=FALSE)

mod_linear<-glm(pex_num ~ trt+age_group, data = poi_study2)
summary(mod_linear)
#significant in trt
pchisq(mod_linear$deviance, df=mod_linear$df.residual, lower.tail=FALSE)
mod_linear$df.residual
mod_linear$deviance

dat_1<-filter(poi_study2,age_group==1)
mod_log1<-glm(recurrent ~ trt, data = dat_1, family = binomial)
#pchisq(mod_linear$deviance, df=mod_linear$df.residual, lower.tail=FALSE)
#plot(recurrent~trt,data=poi_study2)
newdat <- data.frame(trt=seq(min(dat_1$trt), max(dat_1$trt),len=100))
newdat$recu = predict(mod_log1, newdata=newdat, type="response")
plot(recurrent~trt, data=dat_1, col="red4")
lines(recu~trt, newdat, col="green4", lwd=2)

#risk set problme for pwp
length(which(poi_study2$event_pwp==8))

#cox_second
cox2_study2<-study2[1,]
for(i in 1:nrow(select_study2)){
  test2<-NULL
  order2<-NULL
  test2<-filter(study2,tdnccid==study2_id[i])
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
  cox2_study2<-rbind(cox2_study2,order2[event_some,])
}
cox2_study2<-cox2_study2[-1,]
nrow(cox2_study2)
cox2_study2$status<-as.numeric(cox2_study2$pex_num>=2)
cox2_study2$trt[cox2_study2$trt==2]<-0
which(cox2_study2$pex_days==0)
nrow(cox2_study2)
cox2_study2<-cox2_study2[cox2_study2$pex_num>0,]
cox2_study2[203,]
nrow(cox2_study2)
surv_second<-Surv(cox2_study2$pex_days, cox2_study2$status)
surv_second

trt_fit <- survfit(surv_second ~ trt, data=cox2_study2,
                   conf.type="log-log",conf.int=0.95)
library(survminer)
plot1.1<-ggsurvplot(trt_fit, #conf.int = TRUE,
                    #conf.int.style = "step",
                    data = cox2_study2,
                    legend.title="Treatment Arm",
                    legend.labs = 
                      c("Placebo","Azithromycin"),
                    ylab=">=2nd PE free probability",
                    xlab = "Time in days")
plot1.1+ggtitle("Kaplan-Meier Estimate (second event since time 0)")

#plot(opt.ft,xlab="Observed follow-up time (days)",ylab="Survival") 
#title("Kaplan-Meier Survival Estimate (0-second)")
#opt.ft<-survfit(surv_second~1,conf.type="log-log",conf.int=0.95)
##model cox
mod_second<-coxph(surv_second~trt+strata(age_group),
                data=cox2_study2)
summary(mod_second)

#gap_cox_first to second

study2$start_first<-0
gap_study2<-study2[1,]
for(i in 1:nrow(select_study2)){
  test2<-NULL
  order2<-NULL
  test2<-filter(study2,tdnccid==study2_id[i])
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
    if (needed$fu_days>order2[event_some,]$pex_days){
      needed$pex_days<-order2[event_some,]$fu_days
    }
    
  } else {
    event_some<-1
    needed<-order2[event_some,]
  }
  
  gap_study2<-rbind(gap_study2,needed)
}
gap_study2<-gap_study2[-1,]
nrow(gap_study2)
which(gap_study2$pex_days==0)
gap_study2[197,]
gap_study2$status<-as.numeric(gap_study2$pex_num>=2)
gap_study2
which(gap_study2$pex_days==gap_study2$start_first)
gap_study2[37,]
#gap_study2
leftcen2<-gap_study2$pex1_days
leftcen2<-leftcen2[!is.na(leftcen2)]
leftcen2<-leftcen2[leftcen2>0]
mincen2<-min(leftcen2)
length(leftcen2)
mincen2
leftcen2

#for (i in (1:nrow(gap_study2))){
#  if (gap_study2[i,]$pex_num<1){
#    gap_study2[i,]$start_first<-13
#  }
#}

gap_study2$start_first
nrow(gap_study2)
gap_study2<-gap_study2[gap_study2$start_first>=13,]
surv_gap<-with(gap_study2,Surv(start_first,pex_days,status))
surv_gap
mod_gap2<-coxph(surv_gap~trt+strata(age_group),
                data=gap_study2)
mod_gap2
kmfit<- survfit(surv_gap~1,data=gap_study2,subset=(start_first>=13),
                  type="kaplan-meier")
kmfit
a<-gap_study2[gap_study2$start_first>=13,]
a$start_first
min(a$start_first)
a$pex_num
#plot(opt.ft,xlab="Observed follow-up time (days)",ylab="Survival") 
#title("Kaplan-Meier Survival Estimate (gap)")
#opt.ft<-survfit(surv_gap~1,conf.type="log-log",conf.int=0.95)
##model cox
mod_gap<-coxph(surv_gap~trt+strata(age_group),
                  data=gap_study2,subset=(start_first>=13))
summary(mod_gap)
#stratum
gap1<-filter(gap_study2,age_group==1)
surv_gap1<-with(gap1,Surv(start_first,pex_days,status))
surv_gap1
mod_gap1<-coxph(surv_gap1~trt,
               data=gap1)
mod_gap1
cox.zph(mod_gap1)

gap2<-filter(gap_study2,age_group==2)
surv_gap2<-with(gap2,Surv(start_first,pex_days,status))
surv_gap2
mod_gap2<-coxph(surv_gap2~trt,
                data=gap2)
mod_gap2
cox.zph(mod_gap2)

gap3<-filter(gap_study2,age_group==3)
surv_gap3<-with(gap3,Surv(start_first,pex_days,status))
surv_gap3
mod_gap3<-coxph(surv_gap3~trt,
                data=gap3)
mod_gap3
cox.zph(mod_gap3)

gap4<-filter(gap_study2,age_group==4)
surv_gap4<-with(gap4,Surv(start_first,pex_days,status))
surv_gap4
mod_gap4<-coxph(surv_gap4~trt,
                data=gap4)
mod_gap4
cox.zph(mod_gap4)

#create table one
poi_study2$sex=poi_study2$sex+1
vars=c("age_group","sex","pex_num")
factorvars=c("age_group","sex")
library(tableone)
tableOneOpt<-CreateTableOne(vars=vars,strata="trt",
                            data=poi_study2,factorVars = factorvars)
print(tableOneOpt)
poi_study2$sex
poi_study2$trt
poi_study2$Treatment<-as.factor(poi_study2$trt)
#"Azithromycin"
ggplot(poi_study2, aes(x=pex_num, color=Treatment,fill=Treatment)) +
  geom_histogram(position="identity",
                 panel.labs = c("Azithromycin","Placebo"),
                binwidth = 1,alpha=0.5)+
  ggtitle("PEx numbers by treatment group")
  #stat_function(fun=dgamma, args=list(shape=0.5, rate=1))

#ranef without large events
poi_study_small<-poi_study2#[poi_study2$pex_num<4,]
poi_study_small$pex_num[poi_study_small$pex_num>3]<-3

nrow(poi_study_small)
max(poi_study_small$pex_num)
mod_poi2s<-glm(pex_num~age_group+trt, offset=log(fu_days),data=poi_study_small,family=poisson)
summary(mod_poi2s)
pchisq(mod_poi2s$deviance, df=mod_poi2s$df.residual, lower.tail=FALSE)
#p_value=0.0001, fit pretty bad
#sex as confounder (not significant), with interaction, significant


mod_glmm2s<-glmer(pex_num~1+age_group+trt+(1|tdnccid), offset=log(fu_days), family=poisson,data=poi_study_small)
summary(mod_glmm2s)
#sex as confounder (not significant), with interaction, significant
rans<-ranef(mod_glmm2s,condVar = TRUE)
ran1s<-rans$tdnccid
ran1s<-ran1s[,1]
random_effect_small=data.frame(random_effect=exp(ran1s))

effect<-exp(ran1s)

fit.gamma <- fitdist(effect, "gamma", lower = c(0, 0))
fit.norm<-fitdist(ran1s,"norm")
fit.gamma
fit.norm

library(ggplot2)
plot_ran_s<-ggplot(random_effect_small,aes(x=random_effect))+
  geom_histogram(fill="ivory4",color="ivory4",binwidth = 0.1,aes(y=..density..))+
  stat_function(fun=dgamma, args=list(shape=67.11786, rate=66.20331),color="salmon4")+
  stat_function(fun=dlnorm, args=list(mean=0.006258826, sd=0.120162542),color="mediumpurple4")+
  #geom_line(aes(x=exp(ran1), y=dgamma(exp(ran1),fit.gamma$estimate["shape"], fit.gamma$estimate["rate"])), color="salmon4")+
  labs(x = "random effect")+
  dist()
  scale_colour_manual("distribution", values = c("gaussian"="mediumpurple4","gamma"="salmon4"))+
  ggtitle("Histogram for random effect (larger events assigned 3)")

plot_ran_s

fit.params <- fitdistr(exp(ran1), "gamma")
#geom_density() 
mod_quasi2s<-glm(pex_num~age_group+trt, offset=log(fu_days),data=poi_study_small,family=quasipoisson)
summary(mod_quasi2s)
pchisq(mod_quasi2s$deviance, df=mod_quasi2s$df.residual, lower.tail=FALSE)
#sex as confounder (not significant), with interaction, significant
#but it is for 0.1 not 0.05 (p=0.0646)
#p value for goodness of fit is 0.0001, fit pretty bad
#overdispersion 1.7
mod_nb2s<-glm.nb(pex_num~age_group+trt+offset(log(fu_days)),data=poi_study_small)
summary(mod_nb2s)
pchisq(mod_nb2s$deviance, df=mod_nb2s$df.residual, lower.tail=FALSE)
#p value for goodness of fit is 0.602, fit pretty well
#sex as confounder (not significant), with interaction, significant
#overdispersion=1.48
PEx_for_OPTIMIZE_small=poi_study_small$pex_num
hist_study2<-hist(PEx_for_OPTIMIZE_small)

test_disperse<-dispersiontest(mod_poi2s)
test_disperse
