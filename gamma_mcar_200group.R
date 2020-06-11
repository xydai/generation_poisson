#$ -S /usr/local/bin/Rscript
setwd("/home/students/daix4/Thesis/")
library(tidyverse) 
library(plyr)
library(dplyr)
library(poisson)
library(lme4)
library(MASS)
library(survival)
library(random)
source("try_pp_time.R")
source("generate_offset.R")
source("data_simulation_group.R")
#n_run=c(100,200,400,600,800)
#beta0<-0.05
beta1<-0.2
beta2<-0.4
beta3<- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8)
r<-0.002
beta4=0
beta5=1
beta6=0
n=200

  powersag <- rep(0,length(beta3))
  powerspwp <- rep(0,length(beta3))
  powerscox <- rep(0,length(beta3))
  powerspoi_off <- rep(0,length(beta3))
  powersover_off <- rep(0,length(beta3))
  powersnb_off <- rep(0,length(beta3))
  pow_poi_adjusted<-rep(0,length(beta3))
  pow_glmm<-rep(0,length(beta3))
  NUM_off<-rep(0,length(beta3))
  
  #for (j in 1:length(method_run)){
    method="mcar"
    #for(h in 1:length(frail)){
      term="gamma"
      
      beta=NULL
      power=NULL
      beta=c(beta1,beta2,beta3[1],beta4,beta5,beta6)
      parameter=list(5000,n,beta,r,method,0.05)
      power=do.call(power_offset,parameter)
      
      powersag[1]<-power$pow_ag
      powerspwp[1]<-power$pow_pwp
      powerscox[1]<-power$pow_cox
      powerspoi_off[1]<-power$pow_poi_offset
      powersover_off[1]<-power$pow_over_offset
      powersnb_off[1]<-power$pow_nb_offset
      pow_poi_adjusted[1]<-power$pow_poi_adjusted
      pow_glmm[1]<-power$pow_glmm
      
      find_thres<-power$test_poi
      threshold1=run_type1_poi(5000,unlist(find_thres),0.05)
      pow_poi_adjusted[1]<-mean(find_thres<threshold1,na.rm=T)
      
      for(i in 2:length(beta3)){
        beta=NULL
        power=NULL
        
        beta=c(beta1,beta2,beta3[i],beta4,beta5,beta6)
        parameter=list(5000,n,beta,r,method,threshold1)
        power=do.call(power_offset,parameter)
        
        powersag[i]<-power$pow_ag
        powerspwp[i]<-power$pow_pwp
        powerscox[i]<-power$pow_cox
        powerspoi_off[i]<-power$pow_poi_offset
        powersover_off[i]<-power$pow_over_offset
        powersnb_off[i]<-power$pow_nb_offset
        pow_poi_adjusted[i]<-power$pow_poi_adjusted
        pow_glmm[i]<-power$pow_glmm
        NUM_off[i]<-power$nanum
      }
      method2<-c(rep('Anderson Gill',length(beta3)),
                rep('PWP',length(beta3)),
                rep('Cox-PH',length(beta3)),
                rep('Poisson with offset',length(beta3)),
                rep('Overdispersion with offset',length(beta3)),
                rep('Negative Binomial with offset',length(beta3)),
                rep('Poisson GLMM',length(beta3)))
      power<-c(powersag,powerspwp,powerscox,
               powerspoi_off,powersover_off,powersnb_off,pow_glmm)
      
      power_adj<-c(powersag,powerspwp,powerscox,
                   pow_poi_adjusted,powersover_off,powersnb_off,pow_glmm)
      
      powers=data.frame(beta3,power=power,power_adj=power_adj,method_used=method2)
      print(powers)
      setwd("/home/students/daix4/Thesis/")
      write.csv(powers,paste0(term,"_",method,"_",n,".csv"), row.names = T)

