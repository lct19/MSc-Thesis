source('simu_data_cv.R')

##########################
### JMs: current value ###
##########################

library(INLA)

####################################
### convert surv set for poisson ###
####################################

y.surv <- inla.surv(time=SurvDat$endpt, event=SurvDat$event) # surv response

basehaz.prior <- list(model='rw1', diagonal=1e-2, # prior of baseline hazard
                      constr=FALSE, # random walk to fit basehaz
                      n.intervals=24, # num of intervals for each patient
                      scale.model=TRUE,
                      hyper=list(prec=list(prior='pc.prec',
                                           param=c(0.5,0.01),
                                           initial=3)))

Surv.info = inla.coxph(y.surv ~ -1 + trt + time, # expand surv dataset
                       data=list(y.surv=y.surv, 
                                 trt=SurvDat$trt,
                                 time=SurvDat$endpt, # prepared for biomarker
                                 baseEffect=SurvDat$baseEffect,
                                 id=SurvDat$id),
                       control.hazard=basehaz.prior)

SurvNew <- Surv.info$data # save expended surv data

# use median t for each each interval
SurvNew$time <- SurvNew$baseline.hazard.time+0.5*SurvNew$baseline.hazard.length
SurvNew$pt.id <- 1:nrow(SurvNew) # unique id for each patient at each time

##########################################
### construct new data set for JM INLA ###
##########################################

n.l <- nrow(LongDat) # Number of observations in Longi dataset
n.s <- nrow(SurvNew) # Number of observations in expanded surv dataset
NAs.l <- rep(NA, n.l) # create NAs vector
NAs.s <- rep(NA, n.s)
Ones.l <- rep(1, n.l) # create ones vector
Ones.s <- rep(1, n.s)
Zeros.s <- rep(0, n.s) # create zeros vector

Y.joint <- list(
  # Longi, pseudo, Surv
  c(LongDat$Y, NAs.s, NAs.s),
  c(NAs.l, Zeros.s, NAs.s),
  c(NAs.l, NAs.s, SurvNew$y..coxph)
)

covariate.s <- data.frame(
  # covariates only Survival part
  intercept.s = c(NAs.l, NAs.s, Ones.s), # intercept is necessary for poisson
  trt = c(NAs.l, NAs.s, SurvNew$trt),
  
  # start point for each intervals
  basehaz = c(NAs.l, NAs.s, SurvNew$baseline.hazard), # fit basehaz
  cv = c(NAs.l, NAs.s, SurvNew$pt.id) #unique id for each patient at each time
)

# unique id for each patients (same id, same coefficient)
unique.id <- unique(SurvDat$id)
n.patients <- length(unique.id) # number of patients

# Longitudinal
idx.inte <- match(LongDat$id, unique.id) # 1 : N (iid2d in INLA)
idx.slope <- idx.inte + n.patients # N+1 : 2N 
# Survival (same rules as Longi to copy coefficient for each patient)
idx.inte.c <- match(SurvNew$id, unique.id)
idx.slope.c <- idx.inte.c + n.patients


# Longitudinal part and copied terms
covariate.l <- data.frame(
  # Longitudinal part
  intercept.l = c(Ones.l, Ones.s, NAs.s),
  time = c(LongDat$time, SurvNew$time, NAs.s),
  baseEffect.l = c(LongDat$be.L, SurvNew$baseEffect, NAs.s),
  id.inte = c(idx.inte, idx.inte.c, NAs.s), # random intercept (coefficient)
  id.slope = c(idx.slope, idx.slope.c, NAs.s), # random slope (coefficient)
  
  rd.inte.l = c(Ones.l, Ones.s, NAs.s), # covariate
  rd.time.l = c(LongDat$time, SurvNew$time, NAs.s), # covariate
  
  # copied terms
  id.cv = c(NAs.l, SurvNew$pt.id, NAs.s), #save c.v. for each patient at each t
  weight.c = c(NAs.l, rep(-1, n.s), NAs.s)
)

expectation <- list(expect = c(NAs.l, NAs.s, SurvNew$E..coxph)) # expectation

joint.data  <- c(covariate.s, covariate.l, expectation, Surv.info$data.list)
joint.data$Y <- Y.joint # save all information in joint.data


#############
### prior ###
#############

fixed.prior <- list(expand.factor.strategy="inla",
                    mean = 0, prec = 0.0001, # prior for fixed effects
                    mean.intercept = 0, prec.intercept = 0.0001)
associ.prior <- list(beta = list(fixed = FALSE, param = c(0, 0.0001),
                                 initial=0.1)) # prior for association
basehaz.prior <- list(prec = list(hyperid = 54001, 
                                  name = "log precision", short.name = "prec",
                                  initial = 3, fixed = FALSE, 
                                  prior = "pc.prec", param = c(0.5, 0.01), 
                                  to.theta = function(x) log(x), 
                                  from.theta = function(x) exp(x)))

k = 2 # prior for random effects (iidkd, k=2)
rd.prior <- list(theta1 = list(param = c(4, rep(1, k), rep(0, (k*k-k)/2))))


#############
### Model ###
#############

jm.fit <- inla(Y ~ -1 +
                 # Longitudinal submodel
                 intercept.l + time * baseEffect.l +
                 f(id.inte, rd.inte.l, model='iid2d', n=2*n.patients,
                   constr = F, hyper=rd.prior) +
                 f(id.slope, rd.time.l, copy='id.inte') +
                 
                 # persudo submodel
                 f(id.cv, weight.c, model = "iid", # infinite percison
                   hyper = list(prec = list(initial = -6, fixed = TRUE)),
                   constr = F) +
                 
                 #Survival submodel
                 intercept.s + trt +
                 f(basehaz, 
                   model = "rw1", values = baseline.hazard.values, 
                   hyper = basehaz.prior, constr = FALSE, 
                   diagonal = 0.01, scale.model = TRUE) +
                 
                 # copied terms
                 f(cv, copy = "id.cv", 
                   hyper = associ.prior, initial = 0.1),
               
               data = joint.data, family = c("gaussian", "gaussian", "poisson"),
               E = joint.data$expect, inla.mode = "experimental",
               control.compute = list(dic = FALSE,
                                      waic = FALSE,
                                      cpo = FALSE,
                                      config = FALSE,
                                      control.gcpo = list(enable = FALSE,
                                                          num.level.sets = -1,
                                                          correct.hyperpar = TRUE)),
               control.fixed = fixed.prior,
               control.family = list(list(),
                                     list(hyper = list(prec = list(initial = 6, fixed=TRUE))),
                                     list()),
               safe=TRUE,
               control.inla = list(int.strategy = "ccd"))


coxfit <- coxph(Surv(tstart, tstop, status) ~ trt + biomarker,
                data = SurvTdc, cluster = id, timefix=FALSE)
        