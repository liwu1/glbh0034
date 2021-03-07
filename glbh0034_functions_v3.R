##############
## packages ##
##############

require(pROC)
require(mixtools)
require(scales)

###################################
## finite mixture model function ##
###################################

plot_fmm <- function(data) {
  
  plot_x <- seq( min(log(data$optical_density),na.rm=T),max(log(data$optical_density),na.rm=T),0.1)
  
  gauss1a <- dnorm(plot_x,mu1,sig1)
  gauss2a <- dnorm(plot_x,mu2,sig2)
  
  x_axis <- c(0.0001,0.001,0.01,0.05,0.1,0.5,1,1.5,2,4)
  log_axis <- log(x_axis)
  
  hist(log(data$optical_density),breaks=35,main="seropositivity threshold\n(finite mixture model)",cex.main=1,col="grey",border="white",
       xlab="optical density (OD)\n(log scale)",freq=F,ylim=c(0,0.5), xlim=c(-8,4), axes=F,cex.lab=0.9)
  axis(side=1,at=log_axis,labels=x_axis,cex.axis=0.8)
  axis(side=2,at=seq(0,0.5,0.1),cex.axis=0.8)
  abline(v=cutoff,col="red",lwd=2)
  abline(v=mu1,col="blue")
  lines(plot_x,gauss1a,col="blue")
  legend("topleft",legend=c(paste("sero-positive >",round(exp(cutoff),2), "OD"),"negative population"),
         col=c("red","blue"),lty=1, bty="n", cex=0.75, seg.len=1)
  
}

#####################################
## sero-conversion rate functions ##
#####################################

plot_seroprev <- function(plot.age.profile,district,colour) {
  
  max.age <- max(plot.age.profile$age.profiles[,2])+5
  age_vector <- c(1:max.age)
  
  seroprev_x <- plot.age.profile$age.profiles[,2]
  seroprev_y <- plot.age.profile$age.profiles[,3]
  seroprev_y_li <- plot.age.profile$age.profiles[,4]
  seroprev_y_ui <- plot.age.profile$age.profiles[,5]
  
  plot(age_vector,rep(1,max.age),pch=19,ylim=c(0,1.1),col="white",xlim=c(1,50),
       cex.axis=1.1,cex.main=1.5,cex.lab=1.1,axes=F,ylab="",xlab="")
  axis(side=1,at=c(0,10,20,30,40,50))
  axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1,1.1),labels=c(0,0.2,0.4,0.6,0.8,1,""))
  mtext(side=1,"Age, years",line=2.5,font=2,cex=1)
  mtext(side=2,"Sero-prevalence",line=2.5,font=2,cex=1)
  mtext(side=3,district,line=0,font=2,cex=1.5)
  
  points(seroprev_x,seroprev_y,pch=21,col=colour,bg=alpha(colour,0.4),cex=1.5)
  segments(x0=seroprev_x,x1=seroprev_x,y0=seroprev_y_li,seroprev_y_ui,col=colour)
  
}

plot_scr <- function(plot.age.profile,district,colour) {
  
  if (district=="Jinja") row <- 1
  if (district=="Kanungu") row <- 2
  if (district=="Tororo") row <- 3
  
  max.age <- max(plot.age.profile$age.profiles[,2])+5
  age_vector <- c(1:max.age)
  
  seroprev_x <- plot.age.profile$age.profiles[,2]
  seroprev_y <- plot.age.profile$age.profiles[,3]
  seroprev_y_li <- plot.age.profile$age.profiles[,4]
  seroprev_y_ui <- plot.age.profile$age.profiles[,5]
  
  p.lambda.rcm1 <- scr_fit$estimates[row,"lambda.est"]
  p.lambda.rcm1_li <- scr_fit$estimates[row,"lambda.lower"]
  p.lambda.rcm1_ui <- scr_fit$estimates[row,"lambda.upper"]
  p.rho.rcm1 <- scr_fit$estimates[row,"rho.est"]
  
  plot.pred <- theo.seroprev(rep(p.lambda.rcm1,max.age),rep(p.rho.rcm1,max.age))
  plot.pred_li <- theo.seroprev(rep(p.lambda.rcm1_li,max.age),rep(p.rho.rcm1,max.age))
  plot.pred_ui <- theo.seroprev(rep(p.lambda.rcm1_ui,max.age),rep(p.rho.rcm1,max.age))
  
  pol_y <- c(age_vector,rev(age_vector))
  pol_x <- c(plot.pred_li$seroprev,plot.pred_ui$seroprev[length(plot.pred_ui$seroprev):1])
  
  plot(age_vector,rep(1,max.age),pch=19,ylim=c(0,1.1),col="white",xlim=c(1,50),
       cex.axis=1.1,cex.main=1.5,cex.lab=1.1,axes=F,ylab="",xlab="")
  axis(side=1,at=c(0,10,20,30,40,50))
  axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1,1.1),labels=c(0,0.2,0.4,0.6,0.8,1,""))
  mtext(side=1,"Age, years",line=2.5,font=2,cex=1)
  mtext(side=2,"Sero-prevalence",line=2.5,font=2,cex=1)
  mtext(side=3,district,line=0,font=2,cex=1.5)
  
  points(seroprev_x,seroprev_y,pch=21,col=colour,bg=alpha(colour,0.4),cex=1.5)
  segments(x0=seroprev_x,x1=seroprev_x,y0=seroprev_y_li,seroprev_y_ui,col=colour)
  
  polygon(pol_y,pol_x, col=alpha(colour,0.2), border=NA )
  lines(age_vector,plot.pred$seroprev,col=colour,lwd=2)
  
  lambda_lab <- round(p.lambda.rcm1,4)
  text(10,0.16,bquote(lambda*": "*.(round(p.lambda.rcm1,3))*" ["*.(round(p.lambda.rcm1_li,3))*" - "*.(round(p.lambda.rcm1_ui,3))*"]"),adj=0,cex=1)
  text(10,0.1,bquote(rho*": "*.(round(p.rho.rcm1,3))),adj=0,cex=1)
  
}

plot_scr2 <- function(plot.age.profile,district,colour) {
  
  max.age <- max(plot.age.profile$age.profiles[,2])+5
  age_vector <- c(1:max.age)
  
  seroprev_x <- plot.age.profile$age.profiles[,2]
  seroprev_y <- plot.age.profile$age.profiles[,3]
  seroprev_y_li <- plot.age.profile$age.profiles[,4]
  seroprev_y_ui <- plot.age.profile$age.profiles[,5]
  
  p.lambda1 <- scr_fit2$estimates["lambda1"]
  p.lambda2 <- scr_fit2$estimates["lambda2"]
  p.rho <- scr_fit2$estimates["rho"]
  p.time <- scr_fit2$estimates["time.of.change"]
  
  plot.lambda <- unlist(c(rep(p.lambda1,max.age-p.time),rep(p.lambda2,p.time)))
  plot.rho <- unlist(rep(p.rho,max.age))
  
  plot.pred <- theo.seroprev(plot.lambda,plot.rho)
  
  plot(age_vector,rep(1,max.age),pch=19,ylim=c(0,1.1),col="white",xlim=c(1,50),
       cex.axis=1.1,cex.main=1.5,cex.lab=1.1,axes=F,ylab="",xlab="")
  axis(side=1,at=c(0,10,20,30,40,50))
  axis(side=2,at=c(0,0.2,0.4,0.6,0.8,1,1.1),labels=c(0,0.2,0.4,0.6,0.8,1,""))
  mtext(side=1,"Age, years",line=2.5,font=2,cex=1)
  mtext(side=2,"Sero-prevalence",line=2.5,font=2,cex=1)
  mtext(side=3,district,line=0,font=2,cex=1.5)
  
  points(seroprev_x,seroprev_y,pch=21,col=colour,bg=alpha(colour,0.4),cex=1.5)
  segments(x0=seroprev_x,x1=seroprev_x,y0=seroprev_y_li,seroprev_y_ui,col=colour)
  lines(age_vector,plot.pred$seroprev,col=colour,lwd=2)
  
  lambda1_lab <- as.numeric(round(p.lambda1,4))
  lambda2_lab <- round(p.lambda2,4)
  legend("topleft",legend=c(paste0("lambda1: ",lambda1_lab),paste0("lambda2: ",lambda2_lab)),bty="n")
  
}

#######################################
## reverse catalytic model functions ##
#######################################

# confidence interval #

find_ci <- function(results) { 
  
  x <- results$par
  n <- length(x)
  V <- solve(results$hessian)
  lci <- vector(length=n)
  uci <- vector(length=n)
  
  for(i in 1:n){
    
    se <- sqrt(-V[i,i])
    lci[i] <- x[i]-1.96*se
    uci[i] <- x[i]+1.96*se
    
  }
  
  M <- cbind(x, lci, uci)
  colnames(M) <- (c("estimate", "lower", "upper"))
  return(M)
  
}

# log likelihood #

log_lik <- function(x, calc_p, calc_deriv_p, pos, ...) {
  
  p <- calc_p(x, ...)
  ll <- pos*log(p) + (1-pos)*log(1-p)
  return(sum(ll))
  
}

# gradiant #

grad <- function(x, calc_p, calc_deriv_p, pos, ...) {
  
  p <- calc_p(x, ...)
  dp <- calc_deriv_p(x, ...)
  dl_dp <- (pos-p)/(p*(1-p))
  n <- length(x)
  d <- vector(length=n)
  
  for(i in 1:n) {
    di <- dl_dp*dp[,i]
    d[i]=sum(di)
  }
  
  return(d)
  
}

# probability of sero-positivity - function of lambda, rho, and individual age #

calc_p0 <- function(x, age) {
  
  lambda  <- exp(x[1]) 
  rho <- exp(x[2])
  p <- lambda/(lambda+rho)*(1-exp(-(lambda+rho)*age))
  return(p)
  
}

# derivative probability sero-positivity #

calc_deriv_p0 <- function(x, age) {
  
  lambda  <- exp(x[1]) 
  rho  	 <- exp(x[2])
  b <- lambda+rho
  a <- exp(-b*age)
  theta <- lambda/b
  d1 <- lambda*(rho/b^2*(1-a)+theta*age*a)
  d2 <- rho*(-lambda/b^2*(1-a)+theta*age*a)
  
  return(cbind(d1,d2))
  
}

predict_p <- function(ages, results, calc_p, calc_deriv_p, ...) {
  
  V <- -solve(results$hessian)
  x <- results$par
  n <- length(ages)
  pred <- vector(length=n)
  lpred <- vector(length=n)
  upred <- vector(length=n)
  
  for(i in 1:n) {
    
    age <- ages[i]
    
    if(age==0) {
      
      pred[i]=0
      lpred[i]=0
      upred[i]=0
      
    } else {
      
      p <- calc_p(x, age, ...)
      
      # find CI for logit(p) using the delta method #
      
      logit_p <- qlogis(p)
      dp <- calc_deriv_p(x, age, ...)
      deriv_logit <- 1/(p*(1-p))
      d <- vector(length=length(x))
      
      for(j in 1:length(x))
        d[j] <- deriv_logit*dp[j]
      V1 <- d %*% V  %*% d
      se <- sqrt(V1[1,1])
      pred[i] <- p
      lpred[i] <- plogis(logit_p - 1.96*se)
      upred[i] <- plogis(logit_p + 1.96*se)
      
    }
  }
  
  M <- cbind(ages, pred, lpred, upred)
  colnames(M) <- (c("age", "predicted", "lower", "upper"))
  return(M)
  
}

dbinom.mod <- function(x,n,p) lfactorial(n) - lfactorial(x) - lfactorial(n-x) + x*log(p) + (n-x)*log(1-p)

loglik.null <- function(n1,n0) {
  
  p.binom <- sum(n1)/(sum(n1)+sum(n0))
  cat(sum(n1),sum(n0),'\n',sep='\t')
  
  print(p.binom)
  
  lista <- which(p.binom==0)
  p.binom[lista] <- 0.0001
  
  lista <- which(p.binom==1)
  p.binom[lista] <- 0.9991
  
  data <- cbind(n1,n1+n0)
  data <- cbind(data,p.binom)
  
  loglik < -as.numeric(apply(data,1,function(x)dbinom.mod(x=x[1],n=x[2],p=x[3])))
  loglik.total <- sum(loglik)
  return(loglik.total)
  
}

# log likelihood SCR with fixed rho #

loglik.scm.fixed <- function(n1,n0,t,lambda,rho) {
  
  p <- lambda/(lambda+rho)
  g <- lambda+rho
  
  p.binom <- p*(1-exp(-g*t))
  lista <- which(p.binom==0)
  p.binom[lista] <- 0.0001
  
  lista <- which(p.binom==1)
  p.binom[lista] <- 0.9991
  
  data <- cbind(n1,n1+n0)
  data <- cbind(data,p.binom)
  loglik <- as.numeric(apply(data,1,function(x) dbinom.mod(x=x[1],n=x[2],p=x[3])))
  loglik.total <- sum(loglik)
  return(loglik.total)
  
}

# scr fit functions #

fit_model <- function(names, init, predict, calc_p, calc_deriv_p, pos, age, ...) {
  
  results <- optim(par=init, fn=log_lik,gr=grad, calc_p, calc_deriv_p, pos, age, ..., control=list(reltol=1E-15, fnscale=-1), method = "BFGS", hessian = TRUE)
  est <- exp(find_ci(results))
  rownames(est) <- names
  
  if(length(predict)>=1) {
    
    pred <- predict_p(predict, results, calc_p, calc_deriv_p, ...)
    return(list(estimates=est, LL=results$value, pred=pred))
    
  }
  else 
    return(list(estimates=est, LL=results$value))
  
}

# reversible catalytic model with a single force of infection (lambda), sero-reversion (rho) #

rev_cat <- function(age, pos, predict=vector()){
  return(fit_model(names=c("lambda","rho"), init=c(-2, -3), predict=predict, calc_p=calc_p0, calc_deriv_p=calc_deriv_p0, pos=pos, age=age))
}

null.model.analysis <- function(scm.object,analysis='overall') {
  
  if(analysis=='overall') {
    
    tabela <- table(scm.object$age,scm.object$seropos)
    loglik <- loglik.null(n1=tabela[,2],n0=tabela[,1])
    output <- list(loglik.total=loglik,df=1)
    return(output)
    
  }
  
  if(analysis=='split') {
    
    loglik <- 0
    df <- 0
    groups <- unique(scm.object$group)
    
    for(i in groups) {
      
      lista <- which(scm.object$group==i)
      tabela <- table(scm.object$age[lista],scm.object$seropos[lista])
      loglik <- loglik+loglik.null(tabela[,2],tabela[,1])
      
    }
    
    output <- list(loglik.total=loglik,df=length(groups))
    
    return(output)
    
  }
  
  if(analysis!='split'&analysis!='overall')cat("'Analysis' option unknown! Try 'overall' or 'split'\n")
  
}

# creates matrix of individual level data - age, seropositivity, optional stratification (group) #

create.data.object <- function(age,seropos,group=NULL) {
  
  n.age <- length(age)
  n.seropos <- length(seropos)
  
  if(is.null(group)==T) group <- rep(1,n.age)
  if(n.age==n.seropos) out <- list(age=round(age),seropos=seropos,group=as.character(group),model='null') else out <- NULL
  
  return(out)
  
}

# creates matrix with sero-prevalence by age category #

create.age.profile <- function(scm.object,analysis='overall',lag=0.05,cl=0.95) {
  
  age <- scm.object$age
  sero.pos <- scm.object$seropos
  
  if(analysis=='overall') {
    
    output <- age.profile.group(age,sero.pos,lag,cl)
    group <- rep('overall',dim(output)[1])
    output <- data.frame(group=group,output)
    output <- list(age.profiles=output,analysis='overall')
    
    return(output)
    
  }	
  
  if(analysis=='split') {
    
    groups <- unique(scm.object$group)
    output <- c()
    group <- c()
    
    for(i in groups) {
      
      lista <- which(scm.object$group==i)	
      results <- age.profile.group(age[lista],sero.pos[lista],lag,cl)
      group <- c(group,rep(i,dim(results)[1]))
      output <- rbind(output,results)
      
    }
    
    output <- data.frame(group=group,output)
    output <- list(age.profiles=output,analysis='split')
    
    return(output)
    
  }	
  
  if(analysis!='overall'&analysis!='split') {
    
    exist.group <- which(unique(scm.object$group)==analysis)
    
    if(length(exist.group)>0) {
      
      lista <- which(scm.object$group==analysis)	
      output <- age.profile.group(age[lista],sero.pos[lista],lag,cl)
      group <- rep(analysis,dim(output)[1])
      output <- data.frame(group=group,output)
      output <- list(age.profiles=output,analysis=analysis)
      
      return(output)
      
    } else {
      
      cat("ERROR: Group name does not exist in the data set!\t")
      
    }
  }	
}

age.profile.group <- function(age,sero.pos,lag,cl){
  
  list.prop <- seq(lag,1,lag)
  list.quantiles <- unique(round(quantile(age,list.prop)))
  num.points <- length(list.quantiles)
  
  output <- matrix(NA,ncol=4,nrow=num.points)	
  
  list.aux <- which(age<=list.quantiles[1])
  output[1,1] <- median(age[list.aux])
  output[1,2] <- mean(sero.pos[list.aux])	
  output[1,3:4] <- binom.test(sum(sero.pos[list.aux],na.rm=T),length(list.aux),conf.level=cl)$conf.int
  
  for(i in 2:num.points){
    
    list.aux <- which(age>list.quantiles[i-1]&age<=list.quantiles[i])
    output[i,1] <- median(age[list.aux])
    output[i,2] <- mean(sero.pos[list.aux])	
    output[i,3:4] <- binom.test(sum(sero.pos[list.aux],na.rm=T),length(list.aux),conf.level=cl)$conf.int
    
  }
  
  output <- data.frame(output)
  colnames(output) <- c('age','sero.prev','lower','upper')
  return(output)
  
}

simple.rcm.analysis <- function(scm.object,analysis='overall',int.lambda=c(0,1),int.rho=c(0.001,0.250),lag=0.01,age.predict=1:60){
  
  age <- scm.object$age
  sero.pos <- scm.object$seropos
  
  if(analysis=='overall') {
    
    results <- rev_cat(age=age,pos=sero.pos,predict=age.predict)
    output <- data.frame(results$pred)
    tabela <- table(age,sero.pos)
    age.values <- sort(unique(age))
    
    lambda <- results$estimates[1,1]
    rho <- results$estimates[2,1]
    loglik <- loglik.scm.fixed(n1=tabela[,2],n0=tabela[,1],t=age.values,lambda=lambda,rho=rho)
    
    fitted.values <- sapply(age.values,function(x,lambda,rho)lambda/(lambda+rho)*(1-exp(-(lambda+rho)*x)),lambda=lambda,rho=rho)
    fitted.values2 <- sapply(age,function(x,lambda,rho)lambda/(lambda+rho)*(1-exp(-(lambda+rho)*x)),lambda=lambda,rho=rho)
    
    fit.roc <- roc(sero.pos,fitted.values2)		
    
    output <- list(loglik.total=loglik,estimates=results$estimates,df=2,model='M1',analysis='overall',age=age.values,
                 fitted.values=fitted.values,fitted.values2=fitted.values2,roc=auc(fit.roc))
    
    return(output)
    
  }
  
  if(analysis=='split-shared-lambda') {
    
    results <- c()
    groups <- unique(scm.object$group)
    
    for(lambda1 in seq(int.lambda[1],int.lambda[2],lag)) {
      
      loglik <- 0
      est.rho <- c()
      
      for(i in groups) {
        
        lista <- which(scm.object$group==i)
        seropos <- scm.object$seropos[lista]
        age <- scm.object$age[lista]
        
        fit <- my.mle.function.lambda(age,seropos,lambda1)
        loglik <- loglik+fit[1]
        est.rho <- c(est.rho,fit[2])	
        
      }
      
      results <- rbind(results,c(loglik,lambda1,est.rho))		
      
    }
    
    aux <- which.max(results[,1])
    results <- results[aux,]
    estimates <- data.frame(Group=groups,lambda=rep(results[2],length(groups)),rho=results[3:(length(groups)+2)])
    results <- list(loglik.total=results[1],estimates=estimates,df=length(groups)+1,model='M1',analysis='split-shared-lambda')
    return(results)
    
  }
  
  if(analysis=='split-shared-rho') {
    
    results <- c()
    groups <- unique(scm.object$group)
    
    for(rho1 in seq(int.rho[1],int.rho[2],lag)) {
      
      loglik <- 0
      est.lambda <- c()
      
      for(i in groups) {
        
        lista <- which(scm.object$group==i)
        seropos <- scm.object$seropos[lista]
        age <- scm.object$age[lista]
        
        fit <- my.mle.function.rho(age,seropos,rho1)
        loglik <- loglik+fit[1]
        est.lambda <- c(est.lambda,fit[2])	
        
      }
      
      results <- rbind(results,c(loglik,est.lambda,rho1))		
      
    }
    
    aux <- which.max(results[,1])
    results <- results[aux,]
    
    estimates <- data.frame(Group=groups,lambda=results[2:(length(groups)+1)],rho=rep(results[length(results)],length(groups)))
    results <- list(loglik.total=results[1],estimates=estimates,df=length(groups)+1,model='M1',analysis='split-shared-rho')
    return(results)
    
  }
  
  if(analysis=='split-unshared-rho') {
    
    groups <- unique(scm.object$group)
    output <- c()
    
    for(i in groups) {
      
      lista <- which(scm.object$group==i)
      results <- rev_cat(age=age[lista],pos=sero.pos[lista],predict=age.predict)
      tabela <- table(age[lista],sero.pos[lista])
      loglik <- loglik.scm.fixed(n1=tabela[,2],n0=tabela[,1],t=as.numeric(rownames(tabela)),lambda=results$estimates[1,1],rho=results$estimates[2,1])
      output <- rbind(output,c(loglik,matrix(t(results$estimates),nrow=1,byrow=T)))
      
    }
    
    loglik.total <- sum(output[,1])
    output <- data.frame(Group=groups,output[,-1])
    colnames(output) <- c('Group','lambda.est','lambda.lower','lambda.upper','rho.est','rho.lower','rho.upper')
    output <- list(loglik.total=loglik.total,estimates=output,df=2*length(groups),model='M1',analysis='split-unshared-rho')
    return(output)
    
  }
  
  if(analysis!='overall'&analysis!='split-unshared-rho'&analysis!='split-shared-rho') {
    
    results <- c()
    group.exist <- which(unique(scm.object$group)==analysis)
    
    if(length(group.exist)>0) {
      
      lista <- which(scm.object$group==analysis)
      age <- age[lista]
      sero.pos <- sero.pos[lista]
      
      results <- rev_cat(age=age,pos=sero.pos,predict=age.predict)
      output <- data.frame(results$pred)
      tabela <- table(age,sero.pos)
      age.values <- sort(unique(age))
      lambda <- results$estimates[1,1]
      rho <- results$estimates[2,1]
      loglik <- loglik.scm.fixed(n1=tabela[,2],n0=tabela[,1],t=age.values,lambda=lambda,rho=rho)
      
      fitted.values <- sapply(age.values,function(x,lambda,rho)lambda/(lambda+rho)*(1-exp(-(lambda+rho)*x)),lambda=lambda,rho=rho)
      fitted.values2 <- sapply(age,function(x,lambda,rho)lambda/(lambda+rho)*(1-exp(-(lambda+rho)*x)),lambda=lambda,rho=rho)
      fit.roc <- roc(sero.pos,fitted.values2)	
    
      output <- list(loglik.total=loglik,estimates=results$estimates,df=2,model='M1',analysis='overall',age=age.values,fitted.values=fitted.values,fitted.values2=fitted.values2,roc=auc(fit.roc))
      return(output)
      
    } else {
      
      cat('ERROR: Group does not exist in the data set!')	
      
    }
  }
}

two.rcm.analysis<-function(scm.object,analysis='overall',rho.int=c(0.001,0.250),time.int=c(1,10),time.step=1,int.rho=c(0.001,0.250),lag.rho=0.01,time.common=NULL,age.predict=1:60,trace=F){
  
  age <- scm.object$age
  sero.pos <- scm.object$seropos
  
  results <- c()
  
  if(analysis=='overall') {
    
    output<-two.rcm.analysis.overall(age,sero.pos,time.int,time.step,trace)
    
    age.values<-sort(unique(age))	
    lambda<-output$estimates[1:2,1]
    rho<-output$estimates[3,1]
    
    fitted.values <- as.numeric(sapply(age.values,seromodel2.time.of.change,change=output$time.of.change,x=c(lambda,rho)))
    fitted.values2 <- as.numeric(sapply(age,seromodel2.time.of.change,change=output$time.of.change,x=c(lambda,rho)))
    
    fit.roc2 <- roc(sero.pos,fitted.values2)
    output <- list(loglik.total=output$loglik.total,estimates=output$estimates,time.of.change=output$time.of.change,df=4,model='M2',analysis='overall',age=age.values,fitted.values=fitted.values,fitted.values2=fitted.values2,roc=auc(fit.roc2))
    
    return(output)
    
  }
  
  
  if (analysis=='split-shared-rho') {
    
    loglik.max<-(-1)*(10^6)
    
    for (rho1 in seq(int.rho[1],int.rho[2],lag.rho)) {
      
      loglik<-0
      estimates<-c()
      
      j<-1
      
      for (i in unique(scm.object$group)) {	
        
        lista <- which(scm.object$group==i)
        fit <- my.mle.scm2(scm.object$age[lista],scm.object$seropos[lista],lambda=c(0.01,0.005),rho=rho1,time=time.common[j])
        loglik <- loglik+fit$value[1]
        
        estimates<-rbind(estimates,c(i,fit$value[1],fit$par,rho1,time.common[j]))
        j<-j+1	
        
      }
      
      if(loglik.max<loglik){
        
        loglik.max<-loglik
        results<-estimates
        
      }
    }
    
    results<-data.frame(Group=as.character(results[,1]),loglik=as.numeric(results[,2]),lambda1=as.numeric(results[,3]),lambda2=as.numeric(results[,4]),rho=as.numeric(results[,5]),time.of.change=as.numeric(results[,6]))
    return(list(loglik.total=sum(results[,2]),estimates=results[,-2],df=2*dim(results)[1]+1,model='M2',analysis='split-shared-rho'))
    
  }
  
  if (analysis=='split-unshared-rho') {
    
    groups <- unique(scm.object$group)
    output<-c()
    
    for (i in groups) {
      
      lista<-which(scm.object$group==i)
      results<-two.rcm.analysis.overall(age[lista],sero.pos[lista],time.int,time.step,trace)		
      output<-rbind(output,c(results$loglik.total,matrix(t(results$estimates),nrow=1,byrow=T),results$time.of.change))
      
    }
    
    loglik.total<-sum(output[,1])
    output<-data.frame(Group=groups,output[,-1])
    colnames(output)<-c('Group','lambda1.est','lambda1.lower','lambda1.upper','lambda2.est','lambda2.lower','lambda2.upper','rho.est','rho.lower','rho.upper','time.of.change')
    output<-list(loglik.total=loglik.total,estimates=output,df=3*length(groups),model='M2',analysis='split-unshared-rho')
    return(output)
    
  }
  
  if(analysis!='overall'&analysis!='split-shared-rho'&analysis!='split-unshared-rho') {
    
    age <- scm.object$age
    sero.pos <- scm.object$seropos
    results <- c()
    aux <- which(unique(scm.object$group)==analysis)
    
    if (length(aux)>0) {
      
      lista <- which(scm.object$group==analysis)
      age <- age[lista]
      sero.pos <- sero.pos[lista]
      output <- two.rcm.analysis.overall(age,sero.pos,time.int,time.step,trace)
      age.values <- sort(unique(age))	
      lambda <- output$estimates[1:2,1]
      rho <- output$estimates[3,1]
      
      fitted.values <- as.numeric(sapply(age.values,seromodel2.time.of.change,change=output$time.of.change,x=c(lambda,rho)))
      fitted.values2 <- as.numeric(sapply(age,seromodel2.time.of.change,change=output$time.of.change,x=c(lambda,rho)))
      fit.roc <- roc(sero.pos,fitted.values2)	
      output <- list(loglik.total=output$loglik.total,estimates=output$estimates,time.of.change=output$time.of.change,df=4,model='M2',analysis='overall',age=age.values,fitted.values=fitted.values,fitted.values2=fitted.values2,roc=auc(fit.roc))
      return(output)
      
    } else {
      cat('ERROR: Group does not exist in the data set\n') 
    }
  }
}

two.rcm.analysis.overall <- function(age,sero.pos,time.int,time.step,trace) {
  times <- seq(time.int[1],time.int[2],time.step)
  loglik <- (-1)*(10^6)
  tabela <- table(age,sero.pos)
  
  for(i in 1:length(times)){
    results=rev_cat2(age=age , pos=sero.pos, change=times[i], predict=1:60)
    loglik.new <- loglik.scm2.fixed(n1=tabela[,2],n0=tabela[,1],t=as.numeric(rownames(tabela)),lambda=results$estimates[1:2,1],rho=results$estimates[3,1],time.of.change=times[i])
    
    if(loglik<loglik.new) {
      loglik <- loglik.new
      output <- list(loglik.total=loglik,estimates=results$estimates,time.of.change=times[i],df=3)	
    }
    if(trace==T)cat('time of change=',lista[i],', log.likelihood=',loglik,'\n',sep='')
  }
  return(output)
}


# maximum likelihood function for rho #

my.mle.function.rho <- function(age,seropos,rho) {
  
  tabela <- table(age,seropos)
  age <- as.numeric(rownames(tabela))
  mle.fit <- optimize(loglik.scm.fixed,interval=c(10^(-15),10),n1=tabela[,2],n0=tabela[,1],t=age,rho=rho,maximum=T)
  estimates.mle <- c(mle.fit$objective,mle.fit$maximum,rho)
  
  return(c(estimates.mle))
  
}

matrix.prob <- function(t,lambda,rho) {
  
  p1 <- lambda/(lambda+rho)*(1-exp(-(lambda+rho)*t))
  p2 <- rho/(lambda+rho)*(1-exp(-(lambda+rho)*t))
  out <- matrix(c(1-p1,p1,p2,1-p2),ncol=2,nrow=2,byrow=T)
  return(out)	
  
}

theo.seroprev <- function(lambda,rho) {
  
  age.max <- length(lambda)
  prob <- c()
  predict.prob <- c()
  
  for(age in 1:age.max) {
    
    prob <- matrix(c(1,0),ncol=2,nrow=1)
    l1 <- lambda[age.max-age+1]
    r1 <- rho[age.max-age+1]
    
    t.aux <- 1
    
    if(age==1) {
      
      prob <- prob%*%matrix.prob(1,l1,r1)	
      
    } else {
      
      for(x in 1:(age-1)){
        
        if(l1==lambda[age.max-age+x+1]&r1==rho[age.max-age+x+1]) {
      
          t.aux <- t.aux+1				
          
        } else {
          
          if(x!=age) {
            
            prob <- prob%*%matrix.prob(t.aux,l1,r1)
            t.aux <- 1
            l1 <- lambda[age.max-age+1+x]
            r1 <- rho[age.max-age+1+x]
            
          } 
        }	 
      }
      
      prob <- prob%*%matrix.prob(t.aux,l1,r1)				
      
    }
    
    predict.prob<-rbind(predict.prob,prob)
    
  }
  
  out <- cbind(1:age.max,lambda[age.max:1])
  out <- cbind(out,rho[age.max:1])
  out <- cbind(out,predict.prob)
  out <- data.frame(out[,-4])
  colnames(out) <- c('age','lambda','rho','seroprev')
  return(out)
  
}

my.mle.scm2 <- function(age,seropos,lambda,rho,time) {
  tabela <- table(age,seropos)
  age <- as.numeric(rownames(tabela))	
  fit <- optim(par=lambda,fn=loglik.scm2.fixed,n1=tabela[,2],n0=tabela[,1],t=age,time.of.change=time,rho=rho,control=list(fnscale=-1,pgtol=1E-15))
  return(fit)
}

loglik.scm2.fixed <- function(n1,n0,t,lambda,rho,time.of.change) {
  
  p.binom<-as.numeric(sapply(t,seromodel2.time.of.change,change=time.of.change,x=c(lambda,rho)))
  lista<-which(p.binom==0)
  p.binom[lista]<-0.0001
  lista<-which(p.binom==1)
  p.binom[lista]<-0.9991
  
  data<-cbind(n1,n1+n0)
  data<-cbind(data,p.binom)
  loglik<-as.numeric(apply(data,1,function(x)dbinom.mod(x=x[1],n=x[2],p=x[3])))
  loglik.total<-sum(na.omit(loglik))
  
  return(loglik.total)
}

seromodel2.time.of.change <- function(age, change,x) {
  
  lambda1  <- x[1] 
  lambda2  <- x[2] 
  rho  	 <- x[3]
  small <- 1E-6
  b0 <- age-change
  theta1 <- lambda1/(lambda1+rho)
  theta2 <- lambda2/(lambda2+rho)
  p_b <- (age > change)*theta1*(1-exp(-(lambda1+rho)*b0)) 
  
  return((age <= change+small)*(theta2*(1-exp(-(lambda2+rho)*age))) + (age > change+small)*((p_b-theta2)*exp(-(lambda2+rho)*change)+theta2))
}
