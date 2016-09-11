source("../R/jimtools.R")
#recompile() #compile("wtmb.cpp","-O0 -g")
compile("wtmb.cpp")
# Load compiled model
dyn.load(dynlib("wtmb"))
#data <- list( cur_yr=2014L, styr=1970L, endyr=2017L, f_age=3L, l_age=15L) #, ages=as.numeric(3:15),
data <- list( fyr_fsh=1991L, lyr_fsh=2014L,cur_yr=2015L, styr=1970L, endyr=2017L, f_age=3L, l_age=15L, ages=as.numeric(3:15),
	wt_obs=as.matrix(read.table("wtagefsh.dat",as.is=TRUE)),wt_sd=as.matrix(read.table("wtagefshsd.dat",as.is=TRUE)))
nyrs <- data$endyr-data$styr+1 -(data$endyr-data$cur_yr-2)
parameters <- list( L1=35., L2=55., log_K=-.12, log_alpha=-11.5, log_sigma_coh=-2.2, log_sigma_yr=-1.2,coh_eff=rep(0.,nyrs),yr_eff=rep(0.,nyrs))
map = list( log_alpha = factor(c(NA)))
obj  <- MakeADFun(data,parameters=parameters,map=map,random=c("coh_eff","yr_eff")); 
obj$hessian <- TRUE
opt  <- do.call("optim",obj)
get_cor(obj)
dtmp <- get_dt(obj)

#######################################
res1 <- list()
res  <- list()
for (i in 1:15)
{
  iyr <- 2015L - i
  data$cur_yr <- iyr
  nyrs <- data$endyr-data$styr+1 -(data$endyr-data$cur_yr-2)
  parameters <- list( L1=35., L2=55., log_K=-.12, log_alpha=-11.5, log_sigma_coh=-2.2, log_sigma_yr=-1.2,coh_eff=rep(0.,nyrs),yr_eff=rep(0.,nyrs))
  obj  <- MakeADFun(data,parameters=parameters,map=map,random=c("coh_eff","yr_eff")); 
  opt  <- do.call("optim",obj)
  dtmp <- get_dt(obj)
  res1[[i]] <- dtmp 
  res[[i]] <- obj$rep()
}
df <- data.table()
for (i in 1:15)
{
	yrs    <- 1970:(2017-i)
  cur_yr <- 2017-i-2
  yrsfit <- 1991:(cur_yr-1)
  yrsprj <- cur_yr:(cur_yr+2)
  nyrs   <- length(yrsfit) #+ length(yrsprj)
  nyrs
  t <- unlist(res[[i]]$wt_obs[yrsfit-1990,])
  t <- data.frame(t,yr=yrsfit,             run=i,src="obs")
  names(t) <-c(3:15,names(t[14:16]))
  length(yrs)
  res[[i]]
  1991:2013
  dim(res[[i]]$wt_pre)
  r <- unlist(res[[i]]$wt_pre[(yrs-1969),])
  r <- data.frame(r,yr=yrs,run=i,src=c(rep("est",length(yrs)-3),rep("projected",length(yrsprj)) ) )
  # r <- bind_rows(r,data.frame(t(replicate(3,res[[i]]$mnwt)),yr=yrsprj,run=rep(i,3),src=rep("mean",3)) )
  names(r) <-c(3:15,names(r[14:16]))
  df <- rbind(df,gather(t,age,wt,1:13), gather(r,age,wt,1:13))
}
odf <- df %>% filter(src=="obs",age>3,age<9) # ,run<14) 
tdf <- df %>% filter(src=="projected",age>3,age<9) # ,run<14) 
edf <- df %>% filter(run==1,src=="est",age>3,age<9) # ,run<14) 
ggplot(odf, aes(x=yr,y=wt)) + geom_point() + geom_line(data=edf,aes(x=yr,y=wt)) + geom_line(data=tdf,aes(x=yr,y=wt,colour=as.factor(run)))  + 
         facet_grid(age ~ ., scale="free_y") + ylab("Weight (kg)") + xlab("Year") + scale_color_discrete(guide=FALSE) #+ geom_hline(data = df, yintercept=(mnwt))
ggplot(odf, aes(x=(age),y=wt)) + geom_line() 
+ geom_line(data=edf,aes(x=yr,y=wt)) + geom_line(data=tdf,aes(x=yr,y=wt,colour=as.factor(run)))  + facet_grid(age ~ .,scale="free_y") + ylab("Weight (kg)") + xlab("Year") + scale_color_discrete(guide=FALSE) #+ geom_hline(data = df, yintercept=(mnwt))
length(res)
iyr
names(res1)
res1
summary(dtmp)
# obj <- MakeADFun(data=data,parameters=parameters)
#map = list( log_alpha = factor(c(NA)),log_K=factor(c(NA)),log_sigma_yr=factor(c(NA)),log_sigma_coh=factor(c(NA)))
#map = list( log_alpha = factor(c(NA)),log_K=factor(c(NA)))
# opt <- do.call("optim",obj,opt$par)
#obj <- MakeADFun(data,parameters=parameters,random=c("coh_eff","yr_eff"))
# Do estimation 
	,rep(2015,13),rep(2016,13))
[1=="wt_prj",]
summary(obj)

?grep
res <- obj$report()
res$yrs = data$styr:data$endyr-1
res$wt_pre
res$mnlen
res$L1
res$likecomp

plot(res$yrs,res$wt_pre[,1],ylim=c(0,2.2),typ="l")
lines(res$yrs,res$wt_pre[,2])
lines(res$yrs,res$wt_pre[,3])
lines(res$yrs,res$wt_pre[,4])
lines(res$yrs,res$wt_pre[,6])
lines(res$yrs,res$wt_pre[,8])
lines(res$yrs,res$wt_pre[,10])
lines(res$yrs,res$wt_pre[,13])
points(1991:2014,res$wt_obs[,1],pch=19,col="yellow")
points(1991:2014,res$wt_obs[,2],pch=19,col="blue")
points(1991:2014,res$wt_obs[,3],pch=19,col="red")
points(1991:2014,res$wt_obs[,4],pch=19,col="green")
points(1991:2014,res$wt_obs[,6],pch=19,col="salmon")
points(1991:2014,res$wt_obs[,8],pch=19,col="straw")
points(1991:2014,res$wt_obs[,10],pch=19,col="grey")
points(1991:2014,res$wt_obs[,13],pch=19,col="black")
res$wt_pre
res$wt_hat
data$wt_obs