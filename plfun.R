# Simulates prior samples of lambda from a normal distributon
my.rprior<-function(ppar) {
  #simulate prior density 
  sim.state=list()
  sim.state$lambda=rnorm(ppar$N,mean=ppar$lmu,sd=ppar$lsig) 
  return(sim.state)
}


# Evaluates the log likelhood of the current lambda values from the prior 
my.dprior<-function(pstate,ppar) {
  #evaluate log prior density given state/parameter values
    pars=pstate$pars
    lp=sum(dnorm(pstate$lambda,mean=ppar$lmu,sd=ppar$lsig,log=TRUE))
    return(lp)
}


my.dobs<-function(state,dat) {
  #evaluate log-likelihood
  nc=length(dat)
  lla=vector('list',nc) # nc = number of competitions
  for (i in 1:nc) {
    lam=state$lambda[dat[[i]]$o] #lambda vector for the items participating in the ith race, we subset by indexes in race
                                  # The lambdas are in the order that they placed in the race
    
    lp=lpf(lam)    # the log probabilities to evaluate the likelihood for the given order
    lla[[i]]$lp=lp 
    lla[[i]]$ll=sum(lp) # ll = log-likelihood for race o
  }
  ll=sum(sapply(lla, function (x) x$ll)) # sum up the loglikelihoods for all the races to get total log-likelihood
  return(list(ll=ll,lla=lla)) # So we get log likelihood and we get a list of the specific log-likelihoods per race
}

# Log Probability Factor (We use this in the function above)
lpf<-function(lam) { #lambda values are in the order the people placed in the race
  np=length(lam) # Number of players in the race
  lp=numeric(np) 
  for (j in 1:np) {
    lp[j]=lam[j]-log(sum(exp(lam[j:np]))) #  log probability of the jth person placing jth 
  }
  return(lp)
}

#Normalised probabilities
pn<-function(lam) {
  elam=exp(lam)
  return(elam/sum(elam))
}


# Simulated dataset given the current lambda state
my.robs<-function(state,dat) { #need dat as we will simulate with same players in each competition
  #simulate data given pars
  nc=length(dat)
  sdat=vector('list',nc) #simulated data intialisation
  for (i in 1:nc) {
    lam=state$lambda[dat[[i]]$o]# current state lambda in the order of how the participating players placed
    np=length(lam) # number of players
    o=ho=oe=numeric(np) #vector initialised with same length as number of players in competition i
    for (j in 1:np) { # Sample data according to our current state - lambda
      oj=sample(1:np,1,prob=pn(lam)) 
      lam[oj]=-Inf #dont pick this one again
      o[j]=dat[[i]]$o[oj]; ho[j]=dat[[i]]$ho[oj]; oe[j]=dat[[i]]$e[oj] # o = order of players from sampled data
                                                                        # oe = order of senioritys corresponding to o
                                                                         # ignore ho
    }
    sdat[[i]]$o=o; sdat[[i]]$e=oe; #sdat[[i]]$ho=ho; 
  }
  return(sdat)
}

#Proposal distribution- how we are going to move from lambda_t to lambda_t+1
# Geoff here just created a clever symmetric (q(x'|x) = q(x|x')) proposal. So it cancels out in alpha for MCMC
my.q<-function(state,qpars) {
  nstate=list()
  N=length(state$lambda)
  if (runif(1)<0.5) {
    nstate$lambda=rnorm(N,mean=state$lambda,sd=qpars$sda) #normal centers around current lambda_t(i) for each player i
    logQQ=0;
  } else { #otherwise we edit only one lambda value
    i=sample(1:N,1)
    nstate$lambda=state$lambda
    nstate$lambda[i]=rnorm(1,mean=state$lambda[i],sd=qpars$sd1)
    logQQ=0;
  }
  return(list(nstate=nstate,qq=logQQ))
}

my.mcmc<-function(mcmc.pars,model,dat) {
  
  T=mcmc.pars$samples*mcmc.pars$sample.interval # No of total iterations. - Samples = How many to collect, Interval = thinning how many you skip till you get another
  state=mcmc.pars$init.state # Initialize state
  
  pars.in.state=match(model$parnames,names(unlist(state))) # the indices we care about for our current in the state list
  n.pars=length(pars.in.state) # How many params we are tracking in the current state

  X=matrix(NA,mcmc.pars$samples+1,n.pars,dimnames=list(NULL,model$parnames)) # Matrix of samples
  X[1,]=unlist(state)[pars.in.state] #make sure same params always go into same columns of X (avoiding label switching)

  post=npost=list()
  post$ll=model$dobs(state,dat)$ll # same log-likelihood to start off with
  post$lp=model$dprior(state,model$prior.param) # same log-probability of prior to start off with


  for (i in 1:T) {

    prop = mcmc.pars$proposal(state,mcmc.pars$prop.pars) # generate proposal (using q function)
    nstate=prop$nstate # The state generated
    npost$ll=model$dobs(nstate,dat)$ll #log-likelihood given proposal state
    npost$lp=model$dprior(nstate,model$prior.param)# log probabilities of proposal from prior on lmbda

    logMHR=npost$ll+npost$lp-post$ll-post$lp+prop$qq # calculating the log alpha in metropolis hastings

    if (log(runif(1))<logMHR) {
	state=nstate
      post=npost # Otherwise we keep the same state and posterior probabilities
    }
   
    if (!(i%%mcmc.pars$sample.interval)) X[i/mcmc.pars$sample.interval+1,]=unlist(state)[pars.in.state] #basically the thinning (0 is false in R)
  }
  return(as.mcmc(X)) #as.mcmc is the coda function #mcmc(data=X,start=0,end=T,thin=mcmc$sample.interval)
}

# Autocorrelation function
my.acf<-function(X,lag,cols=5) {
  nv=dim(X)[2] # number of variables
  par(mfrow=c(ceiling(nv/cols),cols),oma=c(1,1,1,1));                   #acf plots
  for (i in 1:nv) {
    par(mai=0.2*c(1,1,1,1)); 
    plot(acf(X[,i],lag.max=lag,plot=F),type='l',ann=F,xaxp=c(0,lag,2),yaxp=c(0,1,1)); 
    text(lag/2,0.8,colnames(X)[i])
  }
}

##################################################### (Then below is the same things as above but with covariate params)

# Simulate lambda and beta priors from a normal distribution
my.rprior2<-function(ppar) {
  #simulate prior density 
  sim.state=list()
  sim.state$lambda=rnorm(ppar$N,mean=ppar$lmu,sd=ppar$lsig) 
  sim.state$beta=rnorm(ppar$R,mean=ppar$bmu,sd=ppar$bsig) 
  return(sim.state)
}

# Evaluates log-likelihood of the current lambda/beta values in the MCMC
my.dprior2<-function(pstate,ppar) {
  #evaluate log prior density given state/parameter values
    pars=pstate$pars
    lpl=sum(dnorm(pstate$lambda,mean=ppar$lmu,sd=ppar$lsig,log=TRUE))
    lpb=sum(dnorm(pstate$beta,mean=ppar$bmu,sd=ppar$bsig,log=TRUE))
    return(lpl+lpb)
}

my.dobs2<-function(state,dat) {
  #evaluate log-likelihood
  nc=length(dat)
  lla=vector('list',nc)
  for (i in 1:nc) {
    lam=state$lambda[dat[[i]]$o]
    bet=state$beta[dat[[i]]$e]
    lp=lpf2(lam,bet)    
    lla[[i]]$lp=lp
    lla[[i]]$ll=sum(lp)
  }
  ll=sum(sapply(lla, function (x) x$ll))
  return(list(ll=ll,lla=lla))
}

lpf2<-function(lam,bet) {
  np=length(lam) #should equal length bet
  lp=numeric(np)
  for (j in 1:np) {
    lp[j]=lam[j]+bet[j]-log(sum(exp(lam[j:np]+bet[j:np])))
  }
  return(lp)
}

pn2<-function(lam,bet) {
  elam=exp(lam+bet)
  return(elam/sum(elam))
}


# This just generates random data, dat is taken here so we have the same participating players
my.robs2<-function(state,dat) { #need dat as we will simulate with same players in each competition
  #simulate data given pars
  nc=length(dat)
  sdat=vector('list',nc)
  for (i in 1:nc) {
    lam=state$lambda[dat[[i]]$o]
    bet=state$beta[dat[[i]]$e]
    np=length(lam)
    o=ho=oe=numeric(np)
    for (j in 1:np) {
      oj=sample(1:np,1,prob=pn2(lam,bet)) 
      lam[oj]=-Inf #dont pick this one again
      o[j]=dat[[i]]$o[oj]; ho[j]=dat[[i]]$ho[oj]; oe[j]=dat[[i]]$e[oj]
    }
    sdat[[i]]$o=o; sdat[[i]]$e=oe; #sdat[[i]]$ho=ho;  
  }
  return(sdat)
}

my.q2<-function(state,qpars) {
  nstate=list()
  if (runif(1)<0.5) {#update lambda
    nstate$beta=state$beta
    N=length(state$lambda)
    if (runif(1)<0.5) {
      nstate$lambda=rnorm(N,mean=state$lambda,sd=qpars$sda)
      logQQ=0;
    } else {
      i=sample(1:N,1)
      nstate$lambda=state$lambda
      nstate$lambda[i]=rnorm(1,mean=state$lambda[i],sd=qpars$sd1)
      logQQ=0;
    } 
  } else {
    nstate$lambda=state$lambda
    R=length(state$beta)
    if (runif(1)<0.5) {
      nstate$beta=rnorm(R,mean=state$beta,sd=qpars$sdab)
      logQQ=0;
    } else {
      i=sample(1:R,1)
      nstate$beta=state$beta
      nstate$beta[i]=rnorm(1,mean=state$beta[i],sd=qpars$sd1b)
      logQQ=0;
    }
  }
  return(list(nstate=nstate,qq=logQQ))
}
