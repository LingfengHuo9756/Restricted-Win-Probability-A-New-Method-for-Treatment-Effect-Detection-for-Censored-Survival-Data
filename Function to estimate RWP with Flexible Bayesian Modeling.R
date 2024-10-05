
## Function to estimate different metrics
## Arguments
# hmc.t: fitted model for treatment group(using fit.models from 'survHE' package)
# hmc.c: fitted model for control group(using fit.models from 'survHE' package)
# distr.t: distribution of treatment group
# distr.c: distribution of control group
# tau: end point for rwp and hazard ratio
# tau1: end point for restricted survival time
# delta: threshold for wp or rwp ex. define a win when Y_t>Y_c+delta
# alpha: tie value for RWP (usually it is 1/2)
# iter: number of posterior draw of the distribution parameters for both treatment and control group
# only_gof: when set as 0, output all the metrics; when set as 1, output only AIC, BIC, DIC and DIC2 of the fitted model
# knots: number of internal knots if the treatment or the control group model is rps

## output
# hmc.t: fitted model for treatment group
# hmc.c: fitted model for control group
# wp.mean: estimate wp by plugging in the posterior mean of the model parameters into the wp formula
# wp.vec: estimate wp by plugging in each posterior draw of the model parameters into the wp formula (hence it will be a vector with 'iter' number of values)
# rwp.mean: estimate wp by plugging in the posterior mean of the model parameters into the rwp formula
# rwp.vec: estimate rwp by plugging in each posterior draw of the model parameters into the wp formula (hence it will be a vector with 'iter' number of values)
# wp.nb.mean: net benefit derived from wp
# wp.prob.mean: probability that wp is larger than 0.5
# wp.slo.mean: standardized log oddds from wp
# rwp.nb.mean: net benefit derived from rwp
# rwp.prob.mean: probability that rwp is larger than 0.5 
# rwp.slo.mean: standardized log oddds from rwp 
# wp.delta.mean: estimate wp with delta as the threshold by plugging in the posterior mean of the model parameters into the wp with delta formula
# wp.delta.vec: estimate wp with threshold delta by plugging in each posterior draw of the model parameters into the wp with delta formula (hence it will be a vector with 'iter' number of values)
# wp.delta.nb.mean: net benefit derived from wp with delta
# wp.delta.prob.mean: probability that wp with delta is larger than 0.5 
# wp.delta.slo.mean: standardized log oddds from wp with delta
# rwp.delta.mean: estimate rwp with delta as the threshold by plugging in the posterior mean of the model parameters into the rwp with delta formula
# rwp.delta.vec: estimate rwp with delta by plugging in each posterior draw of the model parameters into the rwp with delta formula (hence it will be a vector with 'iter' number of values)
# rwp.delta.nb.mean: net benefit derived from rwp with delta
# rwp.delta.prob.mean: probability that rwp with delta is larger than 0.5 
# rwp.delta.slo.mean: standardized log oddds from rwp with delta
# hr.mean: estimate the average hazard ratio within tau by plugging in the posterior mean of the model parameters
# hr.vec: estimate the average hazard ratio by plugging in each posterior draw of the model parameters (hence it will be a vector with 'iter' number of values)
# Diff.RMST.mean: estimate the difference of the restricted mean survival time between treatment and control group by plugging in the posterior mean of the model parameters
# Diff.RMST.vec: estimate the difference of the restricted mean survival time between treatment and control group by plugging in each posterior draw of the model parameters (hence it will be a vector with 'iter' number of values)
# AIC, BIC, DIC, DIC2: goodness of fit of the treatment and control group model to the data
# func_T: survival function of the treatment group
# func_C: survival function of the control group

# hr.mean=hr.mean, hr.vec=hr.vec, Diff.RMST.mean=Diff.RMST.mean, Diff.RMST.vec=Diff.RMST.vec,AIC=AIC,BIC=BIC,DIC=DIC,DIC2=DIC2,func_T=func_T,func_C=func_C
surv.metrics = function(hmc.t,hmc.c,distr.t,distr.c,tau,tau1,delta,alpha,iter,only_gof=0,knots=NULL){
  # initialization
  
  if(only_gof==0)
  {
    
    
    # construct the distribution to parameter mapping list
    mapping=list()
    mapping$loglogistic=c("alpha","rate")
    mapping$rps=paste0("gamma[", 1:(knots+2), "]")
    mapping$exponential=c("rate")
    mapping$weibull=c("alpha","scale")
    mapping$weibullPH=c("alpha","scale")
    mapping$lognormal=c("meanlog","alpha")
    mapping$gamma=c("alpha","rate")
    mapping$gompertz=c("alpha","rate")
    mapping$gengamma=c("beta[1]","sigma","Q")
    mapping$genf=c("beta[1]","sigma","Q","P")
    
    # extract the posterior draw of the model parameters
    # for treatment group
    pdat = as.matrix(hmc.t$models[[1]])
    para.t=matrix(nrow=iter,ncol=length(mapping[[distr.t]]))
    for(i in 1:length(mapping[[distr.t]]))
    {
      para.t[,i]=pdat[,mapping[[distr.t]][i]]
    }
    if(distr.t=='weibullPH')
    {
      para.t[,2]=(1/para.t[,2])^(1/para.t[,1])
    }
    
    # extract the knots
    if(distr.t=='rps')
    {
      knot.t = hmc.t$misc$data.stan[[1]]$knots
    }else{
      knot.t=NULL
    }
    
    # for control group
    pdat = as.matrix(hmc.c$models[[1]])
    para.c=matrix(nrow=iter,ncol=length(mapping[[distr.c]]))
    for(i in 1:length(mapping[[distr.c]]))
    {
      para.c[,i]=pdat[,mapping[[distr.c]][i]]
    }
    
    if(distr.c=='weibullPH')
    {
      para.c[,2]=(1/para.c[,2])^(1/para.c[,1])
    }
    
    # extract the knots
    if(distr.c=='rps')
    {
      knot.c = hmc.c$misc$data.stan[[1]]$knots
    }else{
      knot.c=NULL
    }
    
    # WP
    wp.mean=integrate(function(x) exp(logS(x,distr.t,colMeans(para.t),knot.t)+logf(x,distr.c,colMeans(para.c),knot.c)),lower=0,upper=Inf)$value
    wp.vec=sapply(1:iter,function(m) integrate(function(x) exp(logS(x,distr.t,para.t[m,],knot.t)+logf(x,distr.c,para.c[m,],knot.c)),lower=0,upper=Inf)$value)
    # Net benefit
    wp.nb.mean=mean(2*wp.vec-1)
    # Prob that WP is larger than 0.5
    wp.prob.mean=mean(wp.vec>0.5)
    # Standardized log odds
    wp.slo.mean=mean(log(wp.vec/(1-wp.vec)))/sd(log(wp.vec/(1-wp.vec)))
    
    # RWP    
    rwp.mean=integrate(function(x) exp(logS(x,distr.t,colMeans(para.t),knot.t)+logf(x,distr.c,colMeans(para.c),knot.c)),lower=0,upper=tau)$value+alpha*exp(logS(tau,distr.t,colMeans(para.t),knot.t))*exp(logS(tau,distr.c,colMeans(para.c),knot.c))
    rwp.vec=sapply(1:iter,function(m) integrate(function(x) exp(logS(x,distr.t,para.t[m,],knot.t)+logf(x,distr.c,para.c[m,],knot.c)),lower=0,upper=tau)$value+alpha*exp(logS(tau,distr.t,para.t[m,],knot.t))*exp(logS(tau,distr.c,para.c[m,],knot.c)))
    # Net benefit
    rwp.nb.mean=mean(2*rwp.vec-1)
    # Prob that WP is larger than 0.5
    rwp.prob.mean=mean(rwp.vec>0.5)
    # Standardized log odds
    rwp.slo.mean=mean(log(rwp.vec/(1-rwp.vec)))/sd(log(rwp.vec/(1-rwp.vec)))
    
    # WP with delta
    wp.delta.mean=integrate(function(x) exp(logS(x+delta,distr.t,colMeans(para.t),knot.t)+logf(x,distr.c,colMeans(para.c),knot.c)),lower=0,upper=Inf)$value+0.5*(integrate(function(x) exp(logF(x+delta,distr.t,colMeans(para.t),knot.t)+logf(x,distr.c,colMeans(para.c),knot.c)),lower=0,upper=Inf)$value-integrate(function(x) exp(logF(x-delta,distr.t,colMeans(para.t),knot.t)+logf(x,distr.c,colMeans(para.c),knot.c)),lower=delta,upper=Inf)$value)
    wp.delta.vec=sapply(1:iter,function(m) integrate(function(x) exp(logS(x+delta,distr.t,para.t[m,],knot.t)+logf(x,distr.c,para.c[m,],knot.c)),lower=0,upper=Inf)$value+0.5*(integrate(function(x) exp(logF(x+delta,distr.t,para.t[m,],knot.t)+logf(x,distr.c,para.c[m,],knot.c)),lower=0,upper=Inf)$value-integrate(function(x) exp(logF(x-delta,distr.t,para.t[m,],knot.t)+logf(x,distr.c,para.c[m,],knot.c)),lower=delta,upper=Inf)$value))
    # Net benefit
    wp.delta.nb.mean=mean(2*wp.delta.vec-1)
    # Prob that WP.DELTA is larger than 0.5
    wp.delta.prob.mean=mean(wp.delta.vec>0.5)
    # Standardized log odds
    wp.delta.slo.mean=mean(log(wp.delta.vec/(1-wp.delta.vec)))/sd(log(wp.delta.vec/(1-wp.delta.vec)))
    
    
    # RWP with delta
    rwp.delta.mean=integrate(function(x) exp(logS(x+delta,distr.t,colMeans(para.t),knot.t)+logf(x,distr.c,colMeans(para.c),knot.c)),lower=0,upper=tau-delta)$value+0.5*(integrate(function(x) exp(logF(x+delta,distr.t,colMeans(para.t),knot.t)+logf(x,distr.c,colMeans(para.c),knot.c)),lower=0,upper=tau)$value-integrate(function(x) exp(logF(x-delta,distr.t,colMeans(para.t),knot.t)+logf(x,distr.c,colMeans(para.c),knot.c)),lower=delta,upper=delta+tau)$value+exp(logF(tau,distr.t,colMeans(para.t),knot.t))*(exp(logF(tau+delta,distr.c,colMeans(para.c),knot.c))- exp(logF(tau,distr.c,colMeans(para.c),knot.c))))+0.5*exp(logS(tau,distr.t,colMeans(para.t),knot.t))*exp(logS(tau,distr.c,colMeans(para.c),knot.c))+0.5*(integrate(function(x) exp(logS(x+delta,distr.t,colMeans(para.t),knot.t)+logf(x,distr.c,colMeans(para.c),knot.c)),lower=tau-delta,upper=tau)$value+integrate(function(x) exp(logF(x-delta,distr.t,colMeans(para.t),knot.t)+logf(x,distr.c,colMeans(para.c),knot.c)),lower=tau,upper=delta+tau)$value-exp(logF(tau-delta,distr.t,colMeans(para.t),knot.t))* exp(logS(tau,distr.c,colMeans(para.c),knot.c))+exp(logF(tau,distr.t,colMeans(para.t),knot.t))* exp(logS(tau+delta,distr.c,colMeans(para.c),knot.c)))
    rwp.delta.vec=sapply(1:iter,function(m) integrate(function(x) exp(logS(x+delta,distr.t,para.t[m,],knot.t)+logf(x,distr.c,para.c[m,],knot.c)),lower=0,upper=tau-delta)$value+0.5*(integrate(function(x) exp(logF(x+delta,distr.t,para.t[m,],knot.t)+logf(x,distr.c,para.c[m,],knot.c)),lower=0,upper=tau)$value-integrate(function(x) exp(logF(x-delta,distr.t,para.t[m,],knot.t)+logf(x,distr.c,para.c[m,],knot.c)),lower=delta,upper=delta+tau)$value+exp(logF(tau,distr.t,para.t[m,],knot.t))*(exp(logF(tau+delta,distr.c,para.c[m,],knot.c))- exp(logF(tau,distr.c,para.c[m,],knot.c))))+0.5*exp(logS(tau,distr.t,para.t[m,],knot.t))*exp(logS(tau,distr.c,para.c[m,],knot.c))+0.5*(integrate(function(x) exp(logS(x+delta,distr.t,para.t[m,],knot.t)+logf(x,distr.c,para.c[m,],knot.c)),lower=tau-delta,upper=tau)$value+integrate(function(x) exp(logF(x-delta,distr.t,para.t[m,],knot.t)+logf(x,distr.c,para.c[m,],knot.c)),lower=tau,upper=delta+tau)$value-exp(logF(tau-delta,distr.t,para.t[m,],knot.t))* exp(logS(tau,distr.c,para.c[m,],knot.c))+exp(logF(tau,distr.t,para.t[m,],knot.t))* exp(logS(tau+delta,distr.c,para.c[m,],knot.c))))
    # Net benefit
    rwp.delta.nb.mean=mean(2*rwp.delta.vec-1)
    # Prob that RWP.DELTA is larger than 0.5
    rwp.delta.prob.mean=mean(rwp.delta.vec>0.5)
    # Standardized log odds
    rwp.delta.slo.mean=mean(log(rwp.delta.vec/(1-rwp.delta.vec)))/sd(log(rwp.delta.vec/(1-rwp.delta.vec)))
    
    # Difference of RMST
    Diff.RMST.mean=integrate(function(x) exp(logS(x,distr.t,colMeans(para.t),knot.t)),lower=0,upper=tau1)$value-integrate(function(x) exp(logS(x,distr.c,colMeans(para.c),knot.c)),lower=0,upper=tau1)$value
    Diff.RMST.vec=sapply(1:iter,function(m) integrate(function(x) exp(logS(x,distr.t,para.t[m,],knot.t)),lower=0,upper=tau1)$value-integrate(function(x) exp(logS(x,distr.c,para.c[m,],knot.c)),lower=0,upper=tau1)$value)    
    
    # Hazard Ratio
    hr.mean=tryCatch(
      {
        integrate(function(t) exp(log(h(t,distr.t,colMeans(para.t),knot.t))-log(h(t,distr.c,colMeans(para.c),knot.c))),lower=0,upper=tau)$value/tau
      },
      error = function(e) {
        NA
      }
    )
    
    hr.vec=sapply(1:iter, function(m) {
      result <- tryCatch(
        {
          integrate(function(t) exp(log(h(t,distr.t,para.t[m,],knot.t))-log(h(t,distr.c,para.c[m,],knot.c))),lower=0,upper=tau)$value/tau
        },
        error = function(e) {
          NA
        }
      )
      return(result)
    })  
    
    # Store the survival function for plotting
    func_T=function(x) {
      return(exp(logS(x,distr.t,colMeans(para.t),knot.t)))
    }
    
    func_C=function(x) {
      return(exp(logS(x,distr.c,colMeans(para.c),knot.c)))
    }

  }
  
  
  # goodness of fit
  AIC=hmc.t$model.fitting$aic+hmc.c$model.fitting$aic
  BIC=hmc.t$model.fitting$bic+hmc.c$model.fitting$bic
  DIC=hmc.t$model.fitting$dic+hmc.c$model.fitting$dic
  DIC2=hmc.t$model.fitting$dic2+hmc.c$model.fitting$dic2
  
  
  if(only_gof==0)
  {
    return(list(hmc.t=hmc.t,hmc.c=hmc.c,wp.mean=wp.mean,wp.vec=wp.vec,rwp.mean=rwp.mean,rwp.vec=rwp.vec,wp.nb.mean=wp.nb.mean, wp.prob.mean=wp.prob.mean, wp.slo.mean=wp.slo.mean, rwp.nb.mean=rwp.nb.mean, rwp.prob.mean=rwp.prob.mean, rwp.slo.mean=rwp.slo.mean, wp.delta.mean=wp.delta.mean, wp.delta.vec=wp.delta.vec, wp.delta.nb.mean=wp.delta.nb.mean, wp.delta.prob.mean= wp.delta.prob.mean, wp.delta.slo.mean=wp.delta.slo.mean, rwp.delta.mean=rwp.delta.mean, rwp.delta.vec=rwp.delta.vec, rwp.delta.nb.mean=rwp.delta.nb.mean, rwp.delta.prob.mean= rwp.delta.prob.mean, rwp.delta.slo.mean= rwp.delta.slo.mean, hr.mean=hr.mean, hr.vec=hr.vec, Diff.RMST.mean=Diff.RMST.mean, Diff.RMST.vec=Diff.RMST.vec,AIC=AIC,BIC=BIC,DIC=DIC,DIC2=DIC2,func_T=func_T,func_C=func_C))
    
  }else{
    return(list(AIC=AIC,BIC=BIC,DIC=DIC,DIC2=DIC2))
  }
}


### Only calculate rwp
surv.rwp = function(hmc.t,hmc.c,distr.t,distr.c,tau,alpha,iter){
  # initialization
  
    # construct the distribution to parameter mapping list
    mapping=list()
    mapping$loglogistic=c("alpha","rate")
    mapping$rps=c("gamma[1]","gamma[2]","gamma[3]","gamma[4]","gamma[5]")
    mapping$exponential=c("rate")
    mapping$weibull=c("alpha","scale")
    mapping$weibullPH=c("alpha","scale")
    mapping$lognormal=c("meanlog","alpha")
    mapping$gamma=c("alpha","rate")
    mapping$gompertz=c("alpha","rate")
    mapping$gengamma=c("beta[1]","sigma","Q")
    mapping$genf=c("beta[1]","sigma","Q","P")
    
    # extract the posterior draw of the model parameters
    # for treatment group
    pdat = as.matrix(hmc.t$models[[1]])
    para.t=matrix(nrow=iter,ncol=length(mapping[[distr.t]]))
    for(i in 1:length(mapping[[distr.t]]))
    {
      para.t[,i]=pdat[,mapping[[distr.t]][i]]
    }
    if(distr.t=='weibullPH')
    {
      para.t[,2]=(1/para.t[,2])^(1/para.t[,1])
    }
    
    # extract the knots
    if(dist.t=='rps')
    {
      knot.t = hmc.t$misc$data.stan[[1]]$knots
    }else{
      knot.t=NULL
    }
    
    # for control group
    pdat = as.matrix(hmc.c$models[[1]])
    para.c=matrix(nrow=iter,ncol=length(mapping[[distr.c]]))
    for(i in 1:length(mapping[[distr.c]]))
    {
      para.c[,i]=pdat[,mapping[[distr.c]][i]]
    }
    
    if(distr.c=='weibullPH')
    {
      para.c[,2]=(1/para.c[,2])^(1/para.c[,1])
    }
    
    # extract the knots
    if(distr.c=='rps')
    {
      knot.c = hmc.c$misc$data.stan[[1]]$knots
    }else{
      knot.c=NULL
    }
    
    
    # RWP    
    rwp.mean=integrate(function(x) exp(logS(x,distr.t,colMeans(para.t))+logf(x,distr.c,colMeans(para.c))),lower=0,upper=tau)$value+alpha*exp(logS(tau,distr.t,colMeans(para.t)))*exp(logS(tau,distr.c,colMeans(para.c)))
    rwp.vec=sapply(1:iter,function(m) integrate(function(x) exp(logS(x,distr.t,para.t[m,],knot.t)+logf(x,distr.c,para.c[m,],knot.c)),lower=0,upper=tau)$value+alpha*exp(logS(tau,distr.t,para.t[m,],knot.t))*exp(logS(tau,distr.c,para.c[m,],knot.c)))
    # Net benefit
    rwp.nb.mean=mean(2*rwp.vec-1)
    # Prob that WP is larger than 0.5
    rwp.prob.mean=mean(rwp.vec>0.5)
    # Standardized log odds
    rwp.slo.mean=mean(log(rwp.vec/(1-rwp.vec)))/sd(log(rwp.vec/(1-rwp.vec)))

    
    
  
  
  return(list(rwp.mean=rwp.mean,rwp.vec=rwp.vec,rwp.nb.mean=rwp.nb.mean,rwp.prob.mean=rwp.prob.mean,rwp.slo.mean=rwp.slo.mean))
}

logS = function(y,distr,para,knot=NULL){ 
  if(distr=='loglogistic')
  {
    return(pllogis(y,para[1],para[2],lower.tail=F,log.p=TRUE))
  }
  
  if(distr=='rps')
  {
    return(lSrp(y,para,knot))
  }
  
  if(distr=='exponential')
  {
    return(pexp(y,para[1],lower.tail=F,log.p=TRUE))
  }
  
  if(distr=='weibull')
  {
    return(pweibull(y,para[1],para[2],lower.tail=F,log.p=TRUE))
  }
  
  if(distr=='weibullPH')
  {
    return(pweibull(y,para[1],para[2],lower.tail=F,log.p=TRUE))
  }
  
  if(distr=='lognormal')
  {
    return(plnorm(y,para[1],para[2],lower.tail=F,log.p=TRUE))
  }
  
  if(distr=='gamma')
  {
    return(pgamma(y,para[1],para[2],lower.tail=F,log.p=TRUE))
  }
  
  if(distr=='gompertz')
  {
    return(pgompertz(y,para[1],para[2],lower.tail=F,log.p=TRUE))
  }
  
  if(distr=='gengamma')
  {
    return(pgengamma(y,para[1],para[2],para[3],lower.tail=F,log.p=TRUE))
  }
  
  if(distr=='genf')
  {
    return(pgenf(y,para[1],para[2],para[3],para[4],lower.tail=F,log.p=TRUE))
  }
}

logF = function(y,distr,para,knot=NULL){ 
  if(distr=='loglogistic')
  {
    return(pllogis(y,para[1],para[2],lower.tail=T,log.p=TRUE))
  }
  if(distr=='rps')
  {
    return(lCrp(y,para,knot))
  }
  if(distr=='exponential')
  {
    return(pexp(y,para[1],lower.tail=T,log.p=TRUE))
  }
  
  if(distr=='weibull')
  {
    return(pweibull(y,para[1],para[2],lower.tail=T,log.p=TRUE))
  }
  
  if(distr=='weibullPH')
  {
    return(pweibull(y,para[1],para[2],lower.tail=T,log.p=TRUE))
  }
  
  if(distr=='lognormal')
  {
    return(plnorm(y,para[1],para[2],lower.tail=T,log.p=TRUE))
  }
  
  if(distr=='gamma')
  {
    return(pgamma(y,para[1],para[2],lower.tail=T,log.p=TRUE))
  }
  
  if(distr=='gompertz')
  {
    return(pgompertz(y,para[1],para[2],lower.tail=T,log.p=TRUE))
  }
  
  if(distr=='gengamma')
  {
    return(pgengamma(y,para[1],para[2],para[3],lower.tail=T,log.p=TRUE))
  }
  
  if(distr=='genf')
  {
    return(pgenf(y,para[1],para[2],para[3],para[4],lower.tail=T,log.p=TRUE))
  }
}

logf = function(y,distr,para,knot=NULL){
  if(distr=='loglogistic')
  {
    return(dllogis(y,para[1],para[2],log=TRUE))
  }
  
  if(distr=='rps')
  {
    return(lhrp(y,para,knot)+lSrp(y,para,knot))
  }
  if(distr=='exponential')
  {
    return(dexp(y,para[1],log=TRUE))
  }
  if(distr=='weibull')
  {
    return(dweibull(y,para[1],para[2],log=TRUE))
  }
  
  if(distr=='weibullPH')
  {
    return(dweibull(y,para[1],para[2],log=TRUE))
  }
  
  if(distr=='lognormal')
  {
    return(dlnorm(y,para[1],para[2],log=TRUE))
  }
  
  if(distr=='gamma')
  {
    return(dgamma(y,para[1],para[2],log=TRUE))
  }
  
  if(distr=='gompertz')
  {
    return(dgompertz(y,para[1],para[2],log=TRUE))
  }
  
  if(distr=='gengamma')
  {
    return(dgengamma(y,para[1],para[2],para[3],log=TRUE))
  }
  
  if(distr=='genf')
  {
    return(dgenf(y,para[1],para[2],para[3],para[4],log=TRUE))
  }
}

h = function(y,distr,para,knot=NULL){
  if(distr=='loglogistic')
  {
    return(hllogis(y,para[1],para[2]))
  }
  if(distr=='rps')
  {
    return(exp(lhrp(y,para,knot)))
  }
  if(distr=='exponential')
  {
    return(hexp(y,para[1]))
  }
  if(distr=='weibull')
  {
    return(hweibull(y,para[1],para[2]))
  }
  
  if(distr=='weibullPH')
  {
    return(hweibull(y,para[1],para[2]))
  }
  
  if(distr=='lognormal')
  {
    return(hlnorm(y,para[1],para[2]))
  }
  
  if(distr=='gamma')
  {
    return(hgamma(y,para[1],para[2]))
  }
  
  if(distr=='gompertz')
  {
    return(hgompertz(y,para[1],para[2]))
  }
  
  if(distr=='gengamma')
  {
    return(hgengamma(y,para[1],para[2],para[3]))
  }
  
  if(distr=='genf')
  {
    return(hgenf(y,para[1],para[2],para[3],para[4]))
  }
}





### define the distribution functions for RP model
## Log Survival function
lSrp = function(x,gamma,knots) 
{
  gamma=cbind(gamma)
  sapply(x,function(t) -exp(flexsurv::basis(knots,log(t))%*%gamma))
}

## Log CDF
lCrp = function(x,gamma,knots) 
{
  gamma=cbind(gamma)
  sapply(x,function(t) log(1-exp(-exp(flexsurv::basis(knots,log(t))%*%gamma))))
}
## Log Hazard function
lhrp = function(x,gamma,knots) 
{
  gamma=cbind(gamma)
  sapply(x,function(t) -log(t) + log(flexsurv::dbasis(knots,log(t))%*%gamma) + flexsurv::basis(knots,log(t))%*%gamma)
}  

