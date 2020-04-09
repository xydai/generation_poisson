library(tidyverse)
library(plyr)
library(dplyr)
sim_pp2 <- function(t, rate) {
  
  path <- matrix(0, nrow = 1, ncol = 2)
  
  jumps_number <- rpois(1, lambda = rate * t)
  jumps_time <- runif(n = jumps_number, min = 0, max = t) %>% sort()
  
  for(j in seq_along(jumps_time)) {
    jump <- matrix(c(jumps_time[j], path[nrow(path), 2],
                     jumps_time[j], path[nrow(path), 2]  + 1),
                   nrow = 2, ncol = 2, byrow = TRUE)
    path <- rbind(path, jump)
  }
  
  path <- rbind(path,
                c(t, path[nrow(path), 2]))
  
  list(path, jumps_time,jumps_number)
  
}

sim_pp1 <- function(t, rate) { 
  
  path <- matrix(0, nrow = 1, ncol = 2)
  
  jumps_time <- rexp(1, rate)
  
  while(jumps_time[length(jumps_time)] < t) {
    
    jump <- matrix(c(jumps_time[length(jumps_time)], path[nrow(path), 2],
                     jumps_time[length(jumps_time)], path[nrow(path), 2]  + 1),
                   nrow = 2, ncol = 2, byrow = TRUE)
    
    path <- rbind(path, jump)
    
    jumps_time <- c(jumps_time,
                    jumps_time[length(jumps_time)] + rexp(1, rate))
  }
  
  path <- rbind(path,
                c(t, path[nrow(path), 2]))
  
  list(path, jumps_time)
}

generate_pp<-function(tt,rate){
  jumps_number<-rep(0,tt)
  time_record<-rep(0,tt)
  for (t in 1:tt){
    jumps_number[t] <- rpois(1, lambda = rate )
    if (jumps_number[t]>0){
      time_record[t]<-t
    }
  }
  
  sum_number<-sum(jumps_number)
  time_records<-time_record[time_record!=0]
  if (length(time_records)==0){
    time_records<-tt
  } 
  return(list(jumps_number=jumps_number,sum_number=sum_number,time_record=time_records))
}

pp.sim <- function (rate, num.events, t0 = 0) 
{
  x = t0 + cumsum(rexp(n = num.events, rate = rate))
  return(c(t0,x))
  
}

pp.sim(0.1,5)
pp.sim.t <- function (rate, num.events, t0 = 0) 
{
  x = t0 + cumsum(rexp(n = num.events, rate = rate))
  return(c(t0,x))
  
}

simple_pp<-function(tt,rate){
  jumps_number<-rep(0,tt)
  #time_record<-rep(0,tt)
  for (t in 1:tt){
    jumps_number[t] <- rpois(1, lambda = rate)
    #if (jumps_number[t]>0){
    #  time_record[t]<-t
    #}
  }
  
  sum_number<-sum(jumps_number)
  if (is.na(sum_number)){
    sum_number<-0
  }
  #time_records<-time_record[time_record!=0]
  #if (length(time_records)==0){
  #  time_records<-tt
  #} 
  #jumps_number<-jumps_number[jumps_number!=0]
  #return(list(jumps_number=jumps_number,sum_number=sum_number,time_record=time_records))
  return(sum_number=sum_number)
}

new_pp<-function(tt,rate){
  jumps_number<-rep(0,tt)
  time_record<-rep(0,tt)
  for (t in 1:tt){
    jumps_number[t] <- rpois(1, lambda = rate)
    if (jumps_number[t]>0){
      time_record[t]<-t
    }
  }
  #for (t in 1:tt){
  #  jumps_number[t] <- rpois(1, lambda = rate)
  #  if (jumps_number[t]==1){
  #    time_record<-c(time_record,t)
  #  } else if (!is.na(jumps_number[t]) && jumps_number[t]>1){
  #    interval<-jumps_number[t]
  #    for (i in 1:jumps_number[t]){
  #      choice[i]<-i/jumps_number[t]
  #    }
  #    time_record<-c(time_record,t-1+choice)
  #  }
  #}
  
  sum_number<-sum(jumps_number)
  if (is.na(sum_number)){
    sum_number<-0
  }
  
  time_record<-time_record[time_record!=0]
  if (length(time_record)==0){
    time_record<-tt
  } 
  
  #jumps_number<-jumps_number[jumps_number!=0]
  return(list(jumps_number=jumps_number,sum_number=sum_number,time_record=time_record))
  #return(sum_number=sum_number)
}

new_pp2<-function(tt,rate){
  jumps_number<-rep(0,tt)
  time_record<-0
  for (t in 1:tt){
    jumps_number[t] <- rpois(1, lambda = rate)#+rgamma(1,1,50)))
    if (identical(jumps_number[t],1)){
      time_record<-c(time_record,t)
    } else if (jumps_number[t]>=2 && !is.na(jumps_number[t])){
      interval<-jumps_number[t]
      choice<-rep(1,interval)
      for (i in 1:interval){
        choice[i]<-i/interval
      }
       time_record<-c(time_record,t-1+choice)
      }
  }
  
  sum_number<-sum(jumps_number)
  if (is.na(sum_number)){
    sum_number<-0
  }
  
  
  time_record<-time_record[time_record!=0]
  if (length(time_record)==0){
    time_record<-tt
  } 
  
  #jumps_number<-jumps_number[jumps_number!=0]
  return(list(sum_number=sum_number,time_record=time_record))
  #return(sum_number=sum_number)
}

new_pp3<-function(tt,rate,error){
  jumps_number<-rep(0,tt)
  time_record<-0
  r_initial<-0.002
  #error<-rgamma(100,1,100)
  jumps_number[1]<-rpois(1, lambda = rate*r_initial)
  for (t in 2:tt){
      size<-sum(jumps_number[1]:jumps_number[t-1])
      jumps_number[t] <- rpois(1, lambda = rate*error[size+1])
      if (identical(jumps_number[t],1)){
        time_record<-c(time_record,t)
      } else if (jumps_number[t]>=2 && !is.na(jumps_number[t])){
        interval<-jumps_number[t]
        choice<-rep(1,interval)
        for (i in 1:interval){
          choice[i]<-i/interval
        }
        time_record<-c(time_record,t-1+choice)
      }
  }
  
  sum_number<-sum(jumps_number)
  if (is.na(sum_number)){
    sum_number<-0
  }
  
  
  time_record<-time_record[time_record!=0]
  if (length(time_record)==0){
    time_record<-tt
  } 
  
  #jumps_number<-jumps_number[jumps_number!=0]
  return(list(sum_number=sum_number,time_record=time_record))
  #return(sum_number=sum_number)
}

new_pp4<-function(tt,rate){
  jumps_number<-rep(0,tt)
  time_record<-0
  for (t in 1:tt){
    jumps_number[t] <- rpois(1, lambda = rate+rgamma(1,1,100))#+rgamma(1,1,50)))
    if (identical(jumps_number[t],1)){
      time_record<-c(time_record,t)
    } else if (jumps_number[t]>=2 && !is.na(jumps_number[t])){
      interval<-jumps_number[t]
      choice<-rep(1,interval)
      for (i in 1:interval){
        choice[i]<-i/interval
      }
      time_record<-c(time_record,t-1+choice)
    }
  }
  
  sum_number<-sum(jumps_number)
  if (is.na(sum_number)){
    sum_number<-0
  }
  
  
  time_record<-time_record[time_record!=0]
  if (length(time_record)==0){
    time_record<-tt
  } 
  
  #jumps_number<-jumps_number[jumps_number!=0]
  return(list(sum_number=sum_number,time_record=time_record))
  #return(sum_number=sum_number)
}
