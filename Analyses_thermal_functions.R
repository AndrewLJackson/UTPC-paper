
#Thermal Functions Rezende

# Analyses_thermal_function.R

# Four functions to perform analyses and plot thermal performance curves

# 			compute.thermal.curve - computes performance based on parameters q10, CTE, thr and decay
# 			fit.thermal.curve - fits non-linear thermal curve to empirical data (ta vs performance)
# 			plot.thermal.curve - plot thermal curves with different parameters (q10,cte,thr and decay)
#			boot.thermal.curve - "bootstrap" to obtain error estimates of parameters from the error of measurements at each Ta 


# ---------------------------Function 'compute.thermal.curve'--------------------------------------- #

#  Computing estimated performance from q10, q10, CTE, thr and decay 

'compute.thermal.curve' <- 
  function(q10,CTE,thr,decay){	
    ta <- seq(10,50,by=5)
    pf.therm <- CTE*10^(log10(q10)/(10/ta))
    pf.state <- ifelse(ta<thr,1,1-decay*(thr-ta)^2)
    pf <- ifelse(pf.therm*pf.state > 0, pf.therm*pf.state, NA)
    data.frame(na.omit(cbind(ta,pf)))}


# ---------------------------Function 'fit.thermal.curve'--------------------------------------- #

# Nonlinear fitting to analyze empirical thermal performance curves. 

'fit.thermal.curve' <- 
  function(ta,pf,xlab="ta",ylab="pf",plot=TRUE){
    xx <- matrix(,1,13)
    raw.pf <- pf
    pf <- pf/max(pf)
    test <- data.frame(ta=ta,pf=pf)
    q10 <- seq(1.5,3.5,by=0.1) 
    CTE <- seq(0.00,0.1,by=0.02)
    thr <- seq(10,35,by=2)
    decay <- seq(0.000,0.01,by=0.002)
    out <- c(NA,NA,NA,NA,NA)
    for (i in 1:length(q10)){
      for (j in 1:length(CTE)){
        for (l in 1:length(thr)){
          for (k in 1:length(decay)){
            pf.therm <- CTE[j]*10^(log10(q10[i])/(10/ta))
            pf.state <- ifelse(ta<thr[l],1,1-decay[k]*(thr[l]-ta)^2)
            pf.tot <- pf.therm*pf.state
            out <- rbind(out,c(q10[i],CTE[j],thr[l],decay[k],sum((pf.tot-pf)^2)))
          }}}}
    out <- out[order(out[,5]),]
    out.min <- out[which(out[,5] == min(na.omit(out[,5]))),]
    opt <- function(ta,q10,CTE,thr,decay)
    {ifelse(ta<thr,(CTE*10^(log10(q10)/(10/ta))),(CTE*10^(log10(q10)/(10/ta)))*(1-decay*(thr-ta)^2))}
    m.opt <- nls(pf ~ opt(ta, q10, CTE, thr,decay), data = test, start = list(q10 = out.min[1],CTE = out.min[2], thr = out.min[3], decay = out.min[4]), trace = T,control = list(maxiter=5000,warnOnly = TRUE,minFactor = 0))
    r.square <- summary(lm(pf ~ predict(m.opt)))$r.squared	
    q10 <- summary(m.opt)$coefficients[1]
    CTE <- summary(m.opt)$coefficients[2]
    thr <- summary(m.opt)$coefficients[3]
    decay <- summary(m.opt)$coefficients[4]
    ta1 <- seq(0,50,by=0.1)
    pf.therm <- CTE*10^(log10(q10)/(10/ta1))
    pf.state <- ifelse(ta1<thr,1,1-decay*(thr-ta1)^2)
    pf.tot <- pf.therm*pf.state
    ta.opt <- ta1[which(pf.tot==max(pf.tot))]
    pf.tot <- ifelse(pf.tot<0,NA,pf.tot)
    k <- -10*(1-decay*(10/log(q10)+sqrt(1/decay + (10/log(q10))^2)))
    breadth.50 <- as.numeric(breadth.thermal.curve(q10,CTE,thr,decay,breadth=0.5)[3])
    breadth.80 <- as.numeric(breadth.thermal.curve(q10,CTE,thr,decay,breadth=0.8)[3])
    max.y <- max(raw.pf)*max(na.omit(pf.tot))
    
    if(plot==TRUE){
      plot(ta,raw.pf,type="b",xlab=xlab,ylab=ylab,xaxs="i",yaxs="i",bty="l",xlim=c(0,50),ylim=c(-0.02,max.y*1.25),las=1,bg="red",cex=1.3,pch=21)
      points(ta1,max.y*pf.tot/max(na.omit(pf.tot)),type="l",col="red")}
    
    if(m.opt$convInfo$stopMessage == 'converged'){
      xx[1,1] <- q10
      xx[1,2] <- summary(m.opt)$coefficients[2]*max(raw.pf)
      xx[1,3] <- thr
      xx[1,4] <- decay
      xx[1,5] <- thr + (1/decay)^0.5
      xx[1,6] <- max(na.omit(pf.tot))*max(raw.pf)
      xx[1,7] <- ta.opt
      xx[1,8] <- breadth.50
      xx[1,9] <- breadth.80
      xx[1,10] <- r.square
      xx[1,11] <- length(test$pf)
      xx[1,12] <- attr(logLik(m.opt),"df")
      xx[1,13] <- AIC(m.opt)}
    
    colnames(xx) <- c("q10","cte","thr","decay","ctmax","max.pf","ta.opt","breadth.50","breadth.80","r.square","n.samples","K.param","AIC")
    xx <- data.frame(xx)
    list(xx,paste(m.opt$convInfo$stopMessage))}



# -------------------------- Function 'plot.thermal.curve' ------------------------------------- #

# Function to plot different curves manipulating q10, CTE, thr and decay 

'plot.thermal.curve' <- 
  function(q10,CTE,thr,decay,points=FALSE,bound=FALSE,xlab="ta",ylab="pf",col.q10="red",col.th="blue",col="black",lty=1,...){
    q10 <- as.numeric(q10)
    CTE <- as.numeric(CTE)
    thr <- as.numeric(thr)
    decay <- as.numeric(decay)
    ta <- seq(0,80,by=0.1)
    pf.therm <- CTE*10^(log10(q10)/(10/ta))
    pf.state <- ifelse(ta<thr,1,1-decay*(thr-ta)^2)
    pf.tot <- pf.therm*pf.state
    if(bound==TRUE){pf.therm <- pf.therm/max(pf.tot);pf.tot <- ifelse(pf.tot<0,NA,pf.tot/max(pf.tot))}
    else{pf.tot <- ifelse(pf.tot<0,NA,pf.tot)}
    if(points==FALSE){plot(ta,pf.tot,type="l",las=1,xlab=xlab,ylab=ylab,ylim=c(-0.02,max(na.omit(pf.tot))),xaxs="i",yaxs="i",bty="l",col=col,lty=lty)}
    else{points(ta,pf.tot,type="l",col=col,lty=lty)}
    points(ta,pf.state*max(na.omit(pf.tot)),type="l",lty=2,col=col.th)
    points(ta,pf.therm,type="l",lty=2,col=col.q10)}



# -------------------------- Function 'plot.thermal.curve.1' ------------------------------------- #

# Function to plot different curves manipulating q10, CTE, thr and decay 
# Alternative function to return estimated performance (data = TRUE) 

'plot.thermal.curve.1' <- 
  function(q10,CTE,thr,decay,points=TRUE,bound=FALSE,xlab="ta",data=FALSE,ylab="pf",col.q10=NA,col.th=NA,col="black",lty=1,lwd=1,...){
    q10 <- as.numeric(q10)
    CTE <- as.numeric(CTE)
    thr <- as.numeric(thr)
    decay <- as.numeric(decay)
    ta <- seq(-70,70,by=0.1)
    pf.therm <- CTE*10^(log10(q10)/(10/ta))
    pf.state <- ifelse(ta<thr,1,1-decay*(thr-ta)^2)
    pf.tot <- pf.therm*pf.state
    if(data==FALSE){
      if(bound==TRUE){pf.therm <- pf.therm/max(pf.tot);pf.tot <- ifelse(pf.tot<0,NA,pf.tot/max(pf.tot))}
      else{pf.tot <- ifelse(pf.tot<0,NA,pf.tot)}
      if(points==FALSE){
        plot(ta,pf.tot,type="l",las=1,xlab=xlab,ylab=ylab,ylim=c(-0.02,max(na.omit(pf.tot))),xaxs="i",yaxs="i",bty="l",col=col,lty=lty,xaxt="n",yaxt="n",lwd=lwd)}
      else{points(ta,pf.tot,type="l",col=col,lty=lty,lwd=lwd)}
      points(ta,pf.state*max(na.omit(pf.tot)),type="l",lty=2,col=col.th,lwd=lwd)
      points(ta,pf.therm,type="l",lty=2,col=col.q10,lwd=lwd)}
    else{
      if(bound==TRUE){pf.tot <- ifelse(pf.tot<0,0,pf.tot/max(pf.tot))}
      list(ta,ifelse(pf.tot<0,0,pf.tot))}}



# --------------------------- Function 'breadth.thermal.curve' ------------------------------------ #

# Function to obtain the breadth of the thermal curve for a given fraction of maximum performance (e.g., 80% or 50%, etc)

'breadth.thermal.curve' <- function(q10,CTE,thr,decay,breadth){
  q10 <- as.numeric(q10)
  CTE <- as.numeric(CTE)
  thr <- as.numeric(thr)
  decay <- as.numeric(decay)
  ta.opt <- thr-10/(log(q10))+(sqrt(1/decay + (10/log(q10))^2))
  max.pf <- (CTE*10^(log10(q10)/(10/ta.opt)))*(1-decay*(thr-ta.opt)^2)
  ta <- seq(0,50,by=0.01)
  pf.therm <- (CTE/max.pf)*10^(log10(q10)/(10/ta))
  pf.state <- ifelse(ta<thr,1,1-decay*(thr-ta)^2)
  pf.tot <- pf.therm*pf.state
  ta <- ta[range(which(pf.tot > breadth))]
  ta.breadth <- ta[2] - ta[1]
  data.frame(ta.low = ta[1],ta.high = ta[2], ta.breadth = ta.breadth)}



# --------------------------- Function 'boot.thermal.curve' ------------------------------------ #

# Function to obtain parameter error estimates based on sd of empirical measurements
# Involves 2 steps:
# (1) simulating n times data points including error component ('rep' = number of replicates)
#	estimated as sd of empirical measurements (which assumes that the error is normally distributed with mean = 0)   
# (2) Fitting n times these curves and estimating the parameters with error

'boot.thermal.curve' <- 
  function(ta,pf,sd, xlab="ta",ylab="pf",rep=100){
    plot(ta,pf,xlim=c(0,50),ylim=c(0,max(pf*1.3)))
    fit <- matrix(,1,11)
    raw.pf <- pf
    pf <- pf/max(raw.pf)
    sd <- sd/max(raw.pf)
    output <- rep(NA,11) 
    {
      q10_ <- seq(1.5,3.5,by=0.1) 
      CTE_ <- seq(0.00,0.1,by=0.02)
      thr_ <- seq(10,35,by=2)
      decay_ <- seq(0.000,0.01,by=0.002)
      out <- c(NA,NA,NA,NA,NA)
      for (i in 1:length(q10_)){
        for (j in 1:length(CTE_)){
          for (l in 1:length(thr_)){
            for (k in 1:length(decay_)){
              pf.therm <- CTE_[j]*10^(log10(q10_[i])/(10/ta))
              pf.state <- ifelse(ta<thr_[l],1,1-decay_[k]*(thr_[l]-ta)^2)
              pf.tot <- pf.therm*pf.state
              out <- rbind(out,c(q10_[i],CTE_[j],thr_[l],decay_[k],sum((pf.tot-pf)^2)))
            }}}}
      out <- out[order(out[,5]),]
      out.min <- out[which(out[,5] == min(na.omit(out[,5]))),]
    }
    
    while(length(output) < (11*rep)+1)
    {
      pf1 <- rnorm(length(pf),pf,sd)
      test <- data.frame(ta=ta,pf1=pf1)
      
      opt <- function(ta,q10_,CTE_,thr_,decay_){ifelse(ta<thr_,(CTE_*10^(log10(q10_)/(10/ta))),(CTE_*10^(log10(q10_)/(10/ta)))*(1-decay_*(thr_-ta)^2))}
      m.opt <- print(try((as.data.frame(summary(nls(pf1 ~ opt(ta, q10_, CTE_, thr_,decay_), data = test, start = list(q10_ = out.min[1],CTE_ = out.min[2], thr_ = out.min[3],
                                                                                                                      decay_ = 		out.min[4]), trace = F,control = list(maxiter=1000,warnOnly = FALSE,minFactor = 0)))$coefficients[1:4])), TRUE))
      
      if(is.data.frame(m.opt)==TRUE)
      {
        m.opt1 <- nls(pf1 ~ opt(ta, q10_, CTE_, thr_,decay_), data = test, start = list(q10_ = out.min[1],CTE_ = out.min[2], thr_ = out.min[3], 
                                                                                        decay_ = 		out.min[4]), trace = T,control = list(maxiter=1000,warnOnly = FALSE,minFactor = 0))
        RSS <- sum(residuals(m.opt1)^2)
        TSS <- sum((test$pf-mean(test$pf))^2)
        r.square <- 1 - (RSS/TSS)
        q10 <- c(m.opt[[1]])[1]
        CTE <- c(m.opt[[1]])[2]
        thr <- c(m.opt[[1]])[3]
        decay <- c(m.opt[[1]])[4]
        ta1 <- seq(0,50,by=0.1)
        pf.therm <- CTE*10^(log10(q10)/(10/ta1))
        pf.state <- ifelse(ta1<thr,1,1-decay*(thr-ta1)^2)
        pf.tot <- pf.therm*pf.state
        ta.opt <- ta1[which(pf.tot==max(pf.tot))]
        pf.tot <- ifelse(pf.tot<0,NA,pf.tot)
        k <- -10*(1-decay*(10/log(q10)+sqrt(1/decay + (10/log(q10))^2)))
        breadth.50 <- as.numeric(breadth.thermal.curve(q10,CTE,thr,decay,breadth=0.5)[3])
        breadth.80 <- as.numeric(breadth.thermal.curve(q10,CTE,thr,decay,breadth=0.8)[3])
        
        fit[1,1] <- q10
        fit[1,2] <- CTE*max(raw.pf)
        fit[1,3] <- thr
        fit[1,4] <- decay
        fit[1,5] <- thr + (1/decay)^0.5
        fit[1,6] <- max(na.omit(pf.tot))*max(raw.pf)
        fit[1,7] <- ta.opt
        fit[1,8] <- breadth.50
        fit[1,9] <- breadth.80
        fit[1,10] <- r.square
        fit[1,11] <- length(test$pf)
        if(is.na(sum(as.numeric(fit)))==FALSE){
          output <- rbind(output,as.numeric(fit))}
        points(ta,pf1*max(raw.pf),col="black",bg="red",pch=21,cex=1.2)
      }
    }
    colnames(output) <- c("q10","cte","thr","decay","ctmax","max.pf","ta.opt","breadth.50","breadth.80","r.square","n.samples")
    output <- data.frame(output[-1,])
    for(i in 1:rep)
    {plot.thermal.curve(output$q10[i],output$cte[i],output$thr[i],output$decay[i],bound=FALSE,points=TRUE,col.q10="#20202000",col.th="#20202000",col="gray")}
    points(ta,raw.pf,type="b",col="black",bg="red",cex=1.5,pch=21)
    data.frame(output)}