# wt.R
library(ggplot2)
library(tidyr)
library(dplyr)
source("../R/readadmb.R")
d <- read_dat("../empirical/wt.dat")
write_dat(d,"../empirical/wtst.dat")
setwd("../empirical")

mod = 0
res <- list(list())
for (i in 2014:2000){
	mod = mod +1
  d$retyr <- i
  print(i)
  write_dat(d,"wt.dat")
	system("./wt -nox >t ")
	system(paste0("cp wt.rep arc/wt_",i,".rep"))
	res[[mod]] <- read_admb("wt")
}
length(res)
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
  names(r) <-c(3:15,names(r[14:16]))
  df <- rbind(df,gather(t,age,wt,1:13), gather(r,age,wt,1:13))
}
df$age <- as.numeric(df$age)+2
odf <- df %>% filter(src=="obs",age>3,age<9) # ,run<14) 
tdf <- df %>% filter(src=="projected",age>3,age<9) # ,run<14) 
edf <- df %>% filter(run==1,src=="est",age>3,age<9) # ,run<14) 
ggplot(odf, aes(x=yr,y=wt)) + geom_point() + facet_grid(age ~ .,scale="free_y") + ylab("Weight (kg)") + xlab("Year") + geom_line(data=edf,aes(x=yr,y=wt)) + geom_line(data=tdf,aes(x=yr,y=wt,colour=as.factor(run)))  + scale_color_discrete(guide=FALSE)

ggplot(tdf, aes(x=yr,y=wt,colour=as.factor(run))) + geom_line() + geom_point(odf,aes(x=yr,y=wt,colour=as.factor(run))) + facet_grid(. ~ age )


