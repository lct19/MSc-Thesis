###############################################
### Simulation: two longitudinal biomarkers ###
###############################################

library(mvtnorm) # for multivariate normal generation (random effects)

N = 200 # number of individuals

############################
### Longitudinal dataset ###
############################

#### Biomarker 1
beta.L1.0 = 2 # Intercept
beta.L1.1 = 0.2 # slope (time)
beta.L1.2 = -0.3 # baseline effect (hypertension ...)
beta.L1.3 = -0.3 # baseline effect x time
sigma.L1.e = 0.3 # error term (standard error)
sigma.L1.b = 0.4 # random intercept (standard error)
sigma.L1.bt = 0.2 # random slope (standard error)
rho.L1 = -0.2 # correlation continuous intercept/slope
cov.L1 <- sigma.L1.b*sigma.L1.bt*rho.L1 # covariance
D.L1 = matrix(c(sigma.L1.b^2,cov.L1, # variance-covariance matrix
             cov.L1,sigma.L1.bt^2),ncol=2,nrow=2)

#### Biomarker 2
beta.L2.0 = -2 # Intercept
beta.L2.1 = -0.3 # slope (time)
beta.L2.2 = 0.2 # baseline effect (hypertension ...)
beta.L2.3 = 0.3 # baseline effect x time
sigma.L2.e = 0.2 # error term (standard error)
sigma.L2.b = -0.2 # random intercept (standard error)
sigma.L2.bt = 0.3 # random slope (standard error)
rho.L2 = 0.3 # correlation continuous intercept/slope
cov.L2 <- sigma.L2.b*sigma.L2.bt*rho.L2 # covariance
D.L2 = matrix(c(sigma.L2.b^2,cov.L2, # variance-covariance matrix
                cov.L2,sigma.L2.bt^2),ncol=2,nrow=2)

associ.1 = 0.25 # association coefficient
associ.2 = 0.3 # association coefficient
gamma_1 = -0.2 # treatment effects on survival

gap = 0.5 # gap between longitudinal repeated measurements
followup = 12 # follow-up duration
mestime = seq(0, followup, gap) # measurement times
time = rep(mestime, N) # time column
mes_i = followup/gap+1 # number of individual measurements
mes_t = mes_i*N # number of longi measurements
id <- rep(1:N, each = mes_i) # patient id

# random effects generation
MVnorm.L1 <- mvtnorm::rmvnorm(N, rep(0, 2), D.L1)
MVnorm.L2 <- mvtnorm::rmvnorm(N, rep(0, 2), D.L2)

b_i.1 = MVnorm.L1[,1] # random intercept
b_i.L1 <- rep(b_i.1, each=mes_i) # random intercept (repeated for longi dataset)
bt_i.1 = MVnorm.L1[,2] # random slope
bt_i.L1 <- rep(bt_i.1, each=mes_i)

b_i.2 = MVnorm.L2[,1] # random intercept
b_i.L2 <- rep(b_i.2, each=mes_i) # random intercept (repeated for longi dataset)
bt_i.2 = MVnorm.L2[,2] # random slope
bt_i.L2 <- rep(bt_i.2, each=mes_i)

Effect1 <- sample(c(-1,0,1), N, replace=T) # baseline information
be.L1 = rep(Effect1, each=mes_i) # baseline information repeated for longi
Effect2 <- sample(c(-1,0,1), N, replace=T) # baseline information
be.L2 = rep(Effect2, each=mes_i) # baseline information repeated for longi

trt <- sample(c(-1,0,1), N, replace=T) # treatment covariate

## linear predictor
linPred.L1 <- beta.L1.0 + b_i.L1 + (beta.L1.1+bt_i.L1)*time +
  beta.L1.2*be.L1 + beta.L1.3*time*be.L1

linPred.L2 <- beta.L2.0 + b_i.L2 + (beta.L2.1+bt_i.L2)*time +
  beta.L2.2*be.L2 + beta.L2.3*time*be.L2

# observed biomarker values
Y.L1 <- linPred.L1 + rnorm(length(linPred.L1), mean = 0, sd = sigma.L1.e)
Y.L2 <- linPred.L2 + rnorm(length(linPred.L2), mean = 0, sd = sigma.L2.e)
LongDat <- data.frame(id=id, time, be.L1, be.L2, Y.L1, Y.L2) # longitudinal dataset


########################
### Survival dataset ###
########################

# time constant part of linear predictor 1 & 2
linPred.tc.L1 <- beta.L1.0 + b_i.L1 + beta.L1.2*be.L1
linPred.tc.L2 <- beta.L2.0 + b_i.L2 + beta.L2.2*be.L2

integral.f = function(t, alpha, bt1, betaL1.1, betaL1.3, Effect1, associ1,
                      bt2, betaL2.1, betaL2.3, Effect2, associ2){ 
  # integral within Survival function
  exp((alpha-1)*log(t) +
        associ1*((betaL1.1+bt1)*exp(log(t)) + betaL1.3*exp(log(t))*Effect1) + 
        associ2*((betaL2.1+bt2)*exp(log(t)) + betaL2.3*exp(log(t))*Effect2))
}

surv.func = function(t, u, trt, 
                     linPred.tc.L1, betaL1.1, betaL1.3, bt1, Effect1,
                     linPred.tc.L2, betaL2.1, betaL2.3, bt2, Effect2
                     ){
  # time-depend cox model, weibull baseline hazard
  # u, sampled from uniform(0, 1)
  # current value association between the longitudinal and survival
  associ.1 = 0.25 # association coefficient
  associ.2 = 0.3 # association coefficient
  baseScale = 0.25 # baseline hazard scale
  baseShape = 0.5 # baseline hazard shape
  gamma_1 = -0.2 # treatment effects on survival
  
  if (t == 0){
    return(1 - u)
  }else{
    exp(-baseScale*baseShape*
          exp(gamma_1*trt+associ.1*linPred.tc.L1+associ.2*linPred.tc.L2)*
          integrate(function(t){integral.f(t, alpha=baseScale, bt1=bt1,
                                           betaL1.1=betaL1.1, betaL1.3=betaL1.3,
                                           Effect1=Effect1, associ1=associ.1,
                                           bt2=bt2, betaL2.1=betaL2.1,
                                           betaL2.3=betaL2.3, Effect2=Effect2,
                                           associ2=associ.2)},
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
      tryCatch(uniroot(function(t) surv.func(t, u = u[i], trt=trt[i],
                                             linPred.tc.L1=linPred.tc.L1[i], 
                                             betaL1.1=beta.L1.1, 
                                             betaL1.3=beta.L1.3,
                                             bt1=bt_i.L1[i], 
                                             Effect1=Effect1[i],
                                             linPred.tc.L2=linPred.tc.L2[i], 
                                             betaL2.1=beta.L2.1, 
                                             betaL2.3=beta.L2.3,
                                             bt2=bt_i.L2[i], 
                                             Effect2=Effect2[i]), 
                       interval=c(1e-20, 12))$root, # avoid 0
               error=function(e) return(13)) # if t is out of 12, then set 13
  }
  return(output)
}


eventTimes <- survsim(N) # event time (inversion method)
event <- as.numeric(eventTimes < followup) # event indicator
## censoring individuals at end of follow-up (not at random)
endpt <- pmin(eventTimes, 12) # event time > 12 -- censoring
SurvDat <- data.frame(id=1:N, endpt, event, trt, Effect1, Effect2) # survival dataset


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

