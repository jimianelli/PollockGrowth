# wt.R
library(ggplot2)
library(tidyr)
library(dplyr)
source("../R/readadmb.R")
d <- read_dat("../empirical/wt.dat")
d <- read_dat("wt2.dat")
write_dat(d,"../empirical/wtst.dat")
setwd("~/_mymods/PollockGrowth/empirical")
getwd()

i=2015
mod = 0
res <- list(list())
for (i in 2014:2000){
	mod = mod +1
  d$cur_yr <- i
  print(i)
  write_dat(d,"wt.dat")
	system("./wt2 -nox >t ")
	system(paste0("cp wt2.rep arc/wt2_",i,".rep"))
	res[[mod]] <- read_admb("wt2")
}
(res[[1]]$mnwt)
# Compile retrospectives
i=3
df <- data.frame()
yrs <- 1991:2014
for (i in 1:15)
{
  lastyr <- 2017-i-2
  yrsfit <- 1991:lastyr
  yrsprj <- (lastyr+1):(lastyr+3)
  nyrs   <- length(yrsfit) #+ length(yrsprj)
  #t <- unlist(res[[i]]$data[1:nyrs,])
  #t <- data.frame(t,yr=yrs[1:nyrs],             run=i,src="obs")
  t <- unlist(res[[i]]$data)
  t <- data.frame(t,yr=yrs,             run=i,src="obs")
  names(t) <-c(3:15,names(t[14:16]))
  r <- unlist(res[[i]]$W)
  r <- data.frame(r,yr=c(yrsfit,yrsprj),run=i,src=c(rep("est",length(yrsfit)),rep("projected",length(yrsprj)) ) )
  r <- bind_rows(r,data.frame(t(replicate(3,res[[i]]$mnwt)),yr=yrsprj,run=rep(i,3),src=rep("mean",3)) )
  names(r) <-c(3:15,names(r[14:16]))
  names(r)
  df <- rbind(df,gather(t,age,wt,1:13), gather(r,age,wt,1:13))
}
df$age <- as.numeric(df$age)+2
odf <- df %>% filter(src=="obs",age>3,age<9) # ,run<14) 
tdf <- df %>% filter(src=="projected",age>3,age<9) # ,run<14) 
edf <- df %>% filter(run==1,src=="est",age>3,age<9) # ,run<14) 
ggplot(odf, aes(x=yr,y=wt)) + geom_point() + geom_line(data=edf,aes(x=yr,y=wt)) + geom_line(data=tdf,aes(x=yr,y=wt,colour=as.factor(run)))  + facet_grid(age ~ .,scale="free_y") + ylab("Weight (kg)") + xlab("Year") + scale_color_discrete(guide=FALSE) #+ geom_hline(data = df, yintercept=(mnwt))
length(res)

ggplot(tdf, aes(x=yr,y=wt,colour=as.factor(run))) + geom_line() + geom_point(odf,aes(x=yr,y=wt,colour=as.factor(run))) + facet_grid(. ~ age )


