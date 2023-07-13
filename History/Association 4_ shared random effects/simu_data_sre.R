#########################################
### Simulation: shared random effects ###
#########################################

library(mvtnorm) # for multivariate normal generation (random effects)

N = 200 # number of individuals

############################
### Longitudinal dataset ###
############################

beta_0 = 2 # Intercept
beta_1 = 0.2 # slope (time)
beta_2 = -0.3 # baseline effect (hypertension)
beta_3 = -0.3 # baseline effect x time
sigma_e = 0.3 # error term (standard error)
sigma_b = 0.1 # random intercept (standard error)
sigma_bt = 0.1 # random slope (standard error)
rho = -0.2 # correlation continuous intercept/slope
cov <- sigma_b*sigma_bt*rho # covariance
D = matrix(c(sigma_b^2,cov, # variance-covariance matrix
             cov,sigma_bt^2),ncol=2,nrow=2)

gap = 0.5 # gap between longitudinal repeated measurements
followup = 12 # follow-up duration
mestime = seq(0, followup, gap) # measurement times
time = rep(mestime, N) # time column
mes_i = followup/gap+1 # number of individual measurements
mes_t = mes_i*N # number of longi measurements
id <- rep(1:N, each = mes_i) # patient id

# random effects generation
MVnorm <- mvtnorm::rmvnorm(N, rep(0, 2), D)
b_i = MVnorm[,1] # random intercept
b_i.L <- rep(b_i, each=mes_i) # random intercept (repeated for longi dataset)
bt_i = MVnorm[,2] # random slope
bt_i.L <- rep(bt_i, each=mes_i)

baseEffect <- sample(c(-1,0,1), N, replace=T) # baseline information
be.L = rep(baseEffect, each=mes_i) # baseline information repeated for longi

trt <- sample(c(-1,0,1), N, replace=T) # treatment covariate

## linear predictor
linPred <- beta_0+b_i.L+(beta_1+bt_i.L)*time+beta_2*be.L+beta_3*time*be.L
# observed biomarker values
Y <- linPred + rnorm(length(linPred), mean = 0, sd = sigma_e) 
LongDat <- data.frame(id=id, time, be.L, Y) # longitudinal dataset


########################
### Survival dataset ###
########################

linPred.tc <- b_i # time constant part of linear predictor

integral.f = function(t, alpha, bt, gamma_2){ 
  # integral within Survival function
  exp((alpha-1)*log(t) + gamma_2*bt*exp(log(t)))
}

surv.func = function(t, u, linPred.tc, bt, trt){
  # time-depend cox model, weibull baseline hazard
  # u, sampled from uniform(0, 1)
  # current value association between the longitudinal and survival
  associ = 0.3 # association coefficient
  baseScale = 0.25 # baseline hazard scale
  baseShape = 0.5 # baseline hazard shape
  gamma_1 = -0.2 # treatment effects on survival
  
  if (t == 0){
    return(1 - u)
  }else{
    exp(-baseScale*baseShape*exp(gamma_1*trt+associ*linPred.tc)*
          integrate(function(t){integral.f(t, bt=bt, alpha=baseScale,
                                           gamma_2=associ)},
                    0, t)$value) - u
  }
}

survsim = function(N=N){
  # simulate event time for 200 patients
  output = numeric(N) # to save event time
  u <- runif(N) # uniform
  for (i in 1:N){
    # find t let survival(t) - u == 0
    output[i] = 
      tryCatch(uniroot(function(t) surv.func(t, u = u[i], 
                                             linPred.tc=linPred.tc[i], 
                                             bt=bt_i[i], trt=trt[i]), 
                       interval=c(1e-20, 12))$root, # avoid 0
               error=function(e) return(13)) # if t is out of 12, then set 13
  }
  return(output)
}

eventTimes <- survsim(N) # event time (inversion method)
event <- as.numeric(eventTimes < followup) # event indicator
## censoring individuals at end of follow-up (not at random)
endpt <- pmin(eventTimes, 12) # event time > 12 -- censoring
SurvDat <- data.frame(id=1:N, endpt, event, trt, baseEffect) # survival dataset


##################################
### renew Longitudinal dataset ###
##################################

## removing longi measurements after death
ind <- rep(NA, N*length(mes_i))
for (i in 1:N){
  for(j in 1:length(mes_i)){
    if(LongDat[(i-1)*length(mes_i)+j, "time"]<=SurvDat[i,"endpt"]){
      ind[(i-1)*length(mes_i)+j]=1
    } } }
LongDat <- LongDat[!is.na(ind),]

