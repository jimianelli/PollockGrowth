# wt.R
library(ggplot2)
library(tidyverse)
library(data.table)
source("../R/readadmb.R")
mytheme <- theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank())# element_line(colour="grey60", linetype="dashed"))
mytheme <- mytheme + theme(text=element_text(size=18)) + theme(axis.title.x=element_text(size=22) ,axis.title.y=element_text(size=22))
mytheme <- mytheme + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank() )
mytheme <- mytheme + theme( panel.background = element_rect(fill="white"), panel.border = element_rect(colour="black", fill=NA, size=1))
setwd("~/_mymods/PollockGrowth/empirical")
getwd()

i=2015
mod = 0
# RE (random effects) results
reres <- list(list())
# First 15 w/o survey
# Second 15 w/ survey
#system(paste0("cp arc/wt2_no_srv.dat wt_",mod_opt,".dat") )
#system(paste0("cp arc/wt2_no_srv.pin wt_",mod_opt,".pin") )
# system(paste0("cp arc/wt2_with_srv.pin wt_",mod_opt,".pin") )
do_est = FALSE
for (ij in 1:3)
{
  mod_opt <- ifelse(ij==1,"nosrv",ifelse(ij==2,"both",ifelse(ij==3,"yreff","coheff")))
  if (ij>1)
  {
    system(paste0("cp arc/wt2_with_srv.dat wt_",mod_opt,".dat") )
    system(paste0("cp wt_both.pin wt_",mod_opt,".pin") )
  }
  else
  {
    mod_opt="both"
    system(paste0("cp arc/wt2_no_srv.dat wt_",mod_opt,".dat") )
    system(paste0("cp arc/wt2_no_srv.pin wt_",mod_opt,".pin") )
  }
  d <- read_dat(paste0("wt_",mod_opt,".dat"))
  print(mod); print(mod_opt)
  for (i in 2014:2000){
  	mod = mod +1
    d$cur_yr <- i
    print(i)
    write_dat(d,paste0("wt_",mod_opt,".dat"))
    fname <- paste0("arc/wt_",mod_opt,"_",ij,"_",i)
    if (do_est)
    {
  	  system(paste0("./wt_",mod_opt," -nox >t "))
  	  system(paste0("cp wt_",mod_opt,".rep ", fname,".rep"))
      system(paste0("cp wt_",mod_opt,".std ", fname,".std"))
      system(paste0("cp wt_",mod_opt,".par ", fname,".par"))
    }
  	reres[[mod]] <- read_admb(fname)
  }
}
# 1st 15 is w/o survey, second w/ and both, 3rd is w/ and yreff, 4th is w/ and coheff
getwd()
length(reres)

#--------------------------
# Prediction evaluations from means
#--------------------------
# Read in weights by age 3-15
stwt <- as.data.frame(read.table("~/OneDrive/ebswp/data/sampler/cases/ebswp/statwts.dat",as.is=TRUE))
stwt[11:13] <- 0
stwt
mnwt <- data.table(read.table("~/OneDrive/ebswp/data/sampler/cases/ebswp/wtage2016.dat",as.is=TRUE))
names(mnwt) <- 3:15;
mnwt$yr <- 1991:2015; 
mnwt
mnwt.m <- melt(data.table(mnwt), measure.vars = 1:13, variable.name = "age", value.name = "wt")
mnwt.m$age <- as.numeric(mnwt.m$age) +2
setkey(mnwt.m,yr,age) # NOTE the order of this matters...
mnwt.m$mygam <- as.numeric(stwt)
mnwt.m
reres[[45]]
length(reres)
repred <- data.table(subn=integer(),yr=double(),age=double(),pred=double(),obs=double(),pyr=integer(), ayr=double(),stat=double())[]
for (i in 1:length(reres))
{
  mysub <- ifelse(i<=15,"No survey",ifelse(i<=30,"With survey","Year effects only"))
  if(i!=15|i!=30|i!=45)
  {
    iyr <- reres[[i]]$cur_yr
    pred <- data.table(age=rep(3:15,2),yr=rep(c(iyr,iyr+1),each=13))
    dtmp <- data.table(yr=reres[[i]]$yr,reres[[i]]$wt_pre)
    dtmp
    names(dtmp)
    names(dtmp) <- c("yr",3:15)
    pred <- melt( dtmp ,measure.vars = 2:14, variable.name="age",value.name="pred")
    pred$age <- as.numeric(pred$age)+2
    setkey(pred,yr,age)
    pred <- pred[yr %between% c(iyr,iyr+1)] #,.(wt=mean(wt)),age][,wt]
    pred
    setkey(pred,yr,age)
    t     <- merge(pred,mnwt.m )[,.(pred=pred,obs=wt,stat=sum(mygam*(wt-pred)^2) ),.(yr,age)]
    t$ayr <- iyr; t$pyr <- t$yr - iyr; 
    if (i<=15) 
      t$subn <- mysub #"No survey" # i #paste0("Model_",1*integer(i/15),i)
    else
      t$subn <- mysub #"With survey"   # i #paste0("Model_",1*integer(i/15),i)
    repred   <- rbind(repred,t)
  }
}
repred$pyr <- as.factor(repred$pyr)
#repred$subn <- as.factor(rep(c("with","without"),each=15))
setkey(repred,age,pyr)
repred
mnpred <- data.table(subn=integer(),yr=double(),age=double(),pyr=integer(), ayr=double(),stat=double(),pred=double(),obs=double())[]
for (subn in c(1,3,5,10))
{
  for (i in 2001:2014){
    pred <- data.table(age=rep(3:15,2),yr=c(rep(i,13),rep(i+1,13)))
    #if (subn==0)
    #  pred$pred <- mnwt.m[yr == i,.(wt=mean(wt)),age][,wt]
    #else
    pred$pred <- mnwt.m[yr %between% c(i-subn,i-1),.(wt=mean(wt)),age][,wt]
    setkey(pred,yr,age)
    t     <- merge(pred,mnwt.m )[,.(pred=pred,obs=wt,stat=sum(mygam*(wt-pred)^2) ),.(yr,age)]
    t$ayr <- i; t$pyr <- t$yr - i; t$subn <- subn
    t
    mnpred
    mnpred   <- rbind(mnpred,t)
  }
}
mnpred$pyr <- as.factor(mnpred$pyr)
mnpred$subn <- as.factor(mnpred$subn)
mnpred
mnpred[,.(score=sum(stat)),.(subn,pyr)] %>% ggplot(aes(x=subn,y=score,fill=pyr)) + geom_bar(stat="identity",position="dodge") + mytheme #stat_identity())
setkey(mnpred,yr)
mnpred

#repred[subn=="With survey"&age<18&age>10] %>%  ggplot(aes(x=yr,y=pred,colour=as.factor(age))) + labs(y="Body weight (kg)",x="Year") + 
p <- repred[subn=="With survey"&age<8&age>3] %>%  ggplot(aes(x=yr,y=pred,colour=as.factor(age))) + labs(y="Body weight (kg)",x="Year") + 
                  mytheme + theme(panel.grid.major.x = element_line(colour="grey",linetype="dashed")) +
                  scale_x_continuous(breaks=2001:2015) + annotate("text", x=2001, y=0.49,colour="red",         label="Age 4",size=9) + annotate("text", x=2003, y=0.68,colour="limegreen", label="Age 5",size=9) + annotate("text", x=2002, y=0.82,colour="darkcyan",   label="Age 6",size=9) + annotate("text", x=2001, y=0.95,colour="purple",      label="Age 7",size=9) +
                  geom_point(aes(x=yr-.2*(as.numeric(pyr)-1.5),shape=pyr,size=1.2)) + 
                  geom_point(aes(x=yr,y=obs),shape=8,size=4) + guides(fill=FALSE,shape=FALSE,size=FALSE,colour=FALSE) #+
                  geom_line(data=mnpred[subn=="1"&age<8&age>3],aes(x=yr,y=pred,colour=as.factor(age)),size=1) 
print(p)
t <- repred[subn=="With survey"&age<8&age>3]
p <- ggplot(t, aes(x=yr,y=pred,colour=as.factor(age))) + labs(y="Body weight (kg)",x="Year") + 
                  mytheme + theme(panel.grid.major.x = element_line(colour="grey",linetype="dashed")) +
                  scale_x_continuous(breaks=2001:2015) + annotate("text", x=2001, y=0.49,colour="red", label="Age 4",size=9) + annotate("text", x=2003, y=0.68,colour="limegreen", label="Age 5",size=9) + annotate("text", x=2002, y=0.82,colour="darkcyan",   label="Age 6",size=9) + annotate("text", x=2001, y=0.95,colour="purple",      label="Age 7",size=9) +
                  geom_point(aes(x=yr-.2*(as.numeric(pyr)-1.5),shape=pyr,size=1.2)) + 
                  geom_point(aes(x=yr,y=obs),shape=8,size=4) + guides(fill=FALSE,shape=FALSE,size=FALSE,colour=FALSE) +
                  geom_line(data=mnpred[subn=="1"&age<8&age>3], aes(x=yr,y=pred,colour=as.factor(age)),size=1) + 
                  geom_point(data=mnpred[subn=="1"&age<8&age>3], aes(x=yr-.2*(as.numeric(pyr)-1.5),shape=pyr,size=1.2))  
                  #geom_line(data=mnpred[subn=="1"&age<18&age>10],aes(x=yr,y=pred,colour=as.factor(age),size=1.1)) 

print(p)

t <- rbind(mnpred[subn=="1"],repred[subn=="No survey"])
ggplot(t,aes(x=yr,y=stat,col=pyr,as.factor(age))) + geom_point() + mytheme #stat_identity())
ggplot(t[pyr=="0"],aes(x=yr,y=stat,col=subn)) + geom_point() +facet_grid(age~.) +  mytheme #stat_identity())

#--------------------------
# Prediction evaluations all models...
#--------------------------
tt <- rbind(mnpred[,.(score=sum(stat)),.(subn,pyr)],repred[,.(score=sum(stat)),.(subn,pyr)] )
setkey(tt,pyr,subn)
tt
ggplot(tt,aes(x=subn,y=score,fill=pyr)) +labs(x="Model",y="Weighted score",fill="Projection \n year") + geom_bar(stat="identity",position="dodge") + mytheme #stat_identity())
#--------------------------
stwt

tt
repred

stwt <- data.table(read.table("statwts.dat",as.is=TRUE))
mnwt <- data.table(read.table("~/OneDrive/ebswp/data/sampler/cases/ebswp/wtage2016.dat",as.is=TRUE))
names(mnwt) <- 3:15;mnwt$yr <- 1991:2015; 
stwt
mnwt
mnwt.m <- melt(data.table(mnwt), measure.vars = 1:13, variable.name = "age", value.name = "wt")
mnwt.m$age <- as.numeric(mnwt.m$age) +2
mnwt.m$mygam <- as.numeric(stwt)
setkey(mnwt.m,yr,age)
mnwt.m
str(mnwt.m)
i=2001
subn=1
mnpred[subn=="1"&age<8&age>3] %>%  ggplot(aes(x=yr,y=pred,colour=as.factor(age))) + geom_line(aes(size=1.0)) + geom_point(aes(shape=pyr,size=1.2)) + geom_point(aes(x=yr,y=obs),shape=3,size=4) +mytheme

mnpred[subn=="1"&age<8&age>3] %>%  ggplot(aes(x=yr,y=pred,colour=as.factor(age))) + geom_line() + labs(y="Body weight (kg)",x="Year") + 
                         annotate("text", x=2001, y=0.54,colour="red", label="Age 4",size=9) +
                         annotate("text", x=2003, y=0.70,colour="forestgreen", label="Age 5",size=9) +
                         annotate("text", x=2001, y=0.95,colour="purple", label="Age 7",size=9) +
                         geom_point(aes(shape=pyr,size=1.2)) + geom_point(aes(x=yr,y=obs),shape=8,size=4) +mytheme + guides(fill=FALSE,shape=FALSE,size=FALSE,colour=FALSE) 

mnpred[subn=="1"&age<8&age>3] %>%  ggplot(aes(x=yr,y=pred,shape=pyr,colour=as.factor(age))) + geom_point(aes(size=4)) + geom_point(aes(x=yr,y=obs),shape=3,size=4) +mytheme
mnpred 
str(repred)
repred
repred[subn<=15,subn:="No survey"]
/new

#-----------------------------------------------------------
 
(res[[1]]$yr)
dim(res[[8]]$wt_pre)
names(res[[1]] )
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
  t <- unlist(res[[i]]$data)
  t <- data.frame(t,yr=yrs,             run=i,src="obs")
  names(t) <-c(3:15,names(t[14:16]))
  r <- unlist(res[[i]]$W[c(yrsfit,yrsprj)-1969,])
  r
  r <- data.frame(r,yr=c(yrsfit,yrsprj),run=i,src=c(rep("est",length(yrsfit)),rep("projected",length(yrsprj)) ) )
  #r <- bind_rows(r,data.frame(t(replicate(3,res[[i]]$mnwt)),yr=yrsprj,run=rep(i,3),src=rep("mean",3)) )
  names(r) <-c(3:15,names(r[14:16]))
  names(r)
  df <- rbind(df,gather(t,age,wt,1:13), gather(r,age,wt,1:13))
}
odf <- df %>% filter(src=="obs",age>3,age<9) # ,run<14) 
tdf <- df %>% filter(src=="projected",age>3,age<9) # ,run<14) 
edf <- df %>% filter(run==1,src=="est",age>3,age<9) # ,run<14) 
tdf
ggplot(odf, aes(x=yr,y=wt)) + geom_point() + geom_line(data=edf,aes(x=yr,y=wt)) + geom_line(data=tdf,aes(x=yr,y=wt,colour=as.factor(run)))  + facet_grid(age ~ .,scale="free_y") + ylab("Weight (kg)") + xlab("Year") + scale_color_discrete(guide=FALSE) #+ geom_hline(data = df, yintercept=(mnwt))
length(res)

ggplot(tdf, aes(x=yr,y=wt,colour=as.factor(run))) + geom_line() + geom_point(data=odf,aes(x=yr,y=wt,colour=as.factor(run))) + facet_grid(. ~ age )


