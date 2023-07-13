source('simu_data_sre.R')

##############################
### JMs: parametric hazard ###
##############################

library(INLA)

##########################################
### construct new data set for JM INLA ###
##########################################

n.l <- nrow(LongDat) # Number of observations in Longi dataset
n.s <- nrow(SurvDat) # Number of observations in expanded surv dataset
NAs.l <- rep(NA, n.l) # create NAs vector
NAs.s <- rep(NA, n.s)
Ones.l <- rep(1, n.l) # create ones vector
Ones.s <- rep(1, n.s)
Zeros.s <- rep(0, n.s) # create zeros vector

y.surv <- inla.surv(time=c(NAs.l, SurvDat$endpt),
                    event=c(NAs.l, SurvDat$event)) # surv response

Y.joint <- list(
  # Longi, Surv
  c(LongDat$Y, NAs.s),
  y.surv
)

# unique id for each patients (same id, same coefficient)
unique.id <- unique(SurvDat$id)
n.patients <- length(unique.id) # number of patients

# Longitudinal
idx.inte <- match(LongDat$id, unique.id) # 1 : N (iid2d in INLA)
idx.slope <- idx.inte + n.patients # N+1 : 2N 
# Survival (same rules as Longi to copy coefficient for each patient)
idx.inte.c <- match(SurvDat$id, unique.id)
idx.slope.c <- idx.inte.c + n.patients


# Longitudinal part
covariate.l <- data.frame(
  # Longitudinal part
  intercept.l = c(Ones.l, NAs.s),
  time.l = c(LongDat$time, NAs.s),
  baseEffect.l = c(LongDat$be.L, NAs.s),
  id.inte = c(idx.inte, NAs.s), # random intercept (coefficient)
  id.slope = c(idx.slope, NAs.s) # random slope (coefficient)
)

# Survival part and copied terms
covariate.s <- data.frame(
  # covariates only Survival part
  intercept.s = c(NAs.l, Ones.s), # intercept is necessary for poisson
  trt = c(NAs.l, SurvDat$trt),
  time.s = c(NAs.l, SurvDat$endpt), # event time
  
  # shared random effects
  inte.c = c(NAs.l, idx.inte.c),
  slope.c = c(NAs.l, idx.slope.c)
)

joint.data  <- c(covariate.s, covariate.l)
joint.data$Y <- Y.joint # save all information in joint.data


#############
### prior ###
#############

fixed.prior <- list(expand.factor.strategy="inla",
                    mean = 0, prec = 0.01, # prior for fixed effects
                    mean.intercept = 0, prec.intercept = 0.01)
associ.prior <- list(beta = list(fixed = FALSE, param = c(0, 0.01),
                                 initial=0.1)) # prior for association
basehaz.prior <- list(prec = list(hyperid = 54001, 
                                  name = "log precision", short.name = "prec",
                                  initial = 3, fixed = FALSE, 
                                  prior = "pc.prec", param = c(0.5, 0.01), 
                                  to.theta = function(x) log(x), 
                                  from.theta = function(x) exp(x)))

k = 2 # prior for random effects (iidkd, k=2)
rd.prior <- list(theta1 = list(param = c(10, rep(1, k), rep(0, (k*k-k)/2))))


#############
### Model ###
#############

jm.fit <- inla(Y ~ -1 +
                 # Longitudinal submodel
                 intercept.l + time.l * baseEffect.l +
                 f(id.inte, model='iid2d', n=2*n.patients,
                   constr = F, hyper=rd.prior) +
                 f(id.slope, time.l, copy='id.inte') +
                 
                 #Survival submodel
                 trt +
                 f(inte.c, model='iid2d', copy = 'id.inte',
                   hyper=associ.prior) +
                 f(slope.c, time.l, copy='inte.c'),
               
               data = joint.data, family = c("gaussian", "weibull.surv"),
               E = joint.data$expect,
               control.compute = list(dic = TRUE,
                                      waic = TRUE,
                                      cpo = TRUE,
                                      config = TRUE),
               control.fixed = fixed.prior,
               control.inla=list(int.strategy="eb"))


