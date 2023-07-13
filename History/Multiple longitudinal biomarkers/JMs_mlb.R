source('simu_data_mlb.R')

#############################################
### JMs: multiple longitudinal biomarkers ###
#############################################

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

Surv.info = inla.coxph(y.surv ~ -1 + trt + time + time*0, # expand surv dataset
                       data=list(y.surv=y.surv, 
                                 trt=SurvDat$trt,
                                 time=SurvDat$endpt, # prepared for biomarker
                                 Effect1=SurvDat$Effect1,
                                 Effect2=SurvDat$Effect2,
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
  # Longi1, pseudo1, Longi2, pseudo2 Surv
  c(LongDat$Y.L1, NAs.s, NAs.l, NAs.s, NAs.s),
  c(NAs.l, Zeros.s, NAs.l, NAs.s, NAs.s),
  c(NAs.l, NAs.s, LongDat$Y.L2, NAs.s, NAs.s),
  c(NAs.l, NAs.s, NAs.l, Zeros.s, NAs.s),
  c(NAs.l, NAs.s, NAs.l, NAs.s, SurvNew$y..coxph)
)

covariate.s <- data.frame(
  # covariates only Survival part
  intercept.s = c(NAs.l, NAs.s, NAs.l, NAs.s, Ones.s), # intercept is necessary for poisson
  trt = c(NAs.l, NAs.s, NAs.l, NAs.s, SurvNew$trt),
  
  # start point for each intervals
  basehaz = c(NAs.l, NAs.s, NAs.l, NAs.s, SurvNew$baseline.hazard), # fit basehaz
  cv.1 = c(NAs.l, NAs.s, NAs.l, NAs.s, SurvNew$pt.id), #unique id for each patient at each time
  cv.2 = c(NAs.l, NAs.s, NAs.l, NAs.s, SurvNew$pt.id) #unique id for each patient at each time
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
  # Longitudinal biomarker 1
  intercept.l1 = c(Ones.l, Ones.s, NAs.l, NAs.s, NAs.s),
  time1 = c(LongDat$time, SurvNew$time, NAs.l, NAs.s, NAs.s),
  Effect1 = c(LongDat$be.L1, SurvNew$Effect1, NAs.l, NAs.s, NAs.s),
  id.inte.1 = c(idx.inte, idx.inte.c, NAs.l, NAs.s, NAs.s), #rd.inte coefficient
  id.slope.1 = c(idx.slope, idx.slope.c, NAs.l, NAs.s, NAs.s), # coefficient
  rd.inte.l1 = c(Ones.l, Ones.s, NAs.l, NAs.s, NAs.s), # random intercept
  rd.time.l1 = c(LongDat$time, SurvNew$time, NAs.l, NAs.s, NAs.s), #random slope
  
  # Longitudinal biomarker 1
  intercept.l2 = c(NAs.l, NAs.s, Ones.l, Ones.s, NAs.s),
  time2 = c(NAs.l, NAs.s, LongDat$time, SurvNew$time, NAs.s),
  Effect2 = c(NAs.l, NAs.s, LongDat$be.L2, SurvNew$Effect2, NAs.s),
  id.inte.2 = c(NAs.l, NAs.s, idx.inte, idx.inte.c, NAs.s), #rd.inte coefficient
  id.slope.2 = c(NAs.l, NAs.s, idx.slope, idx.slope.c, NAs.s), # coefficient
  rd.inte.l2 = c(NAs.l, NAs.s, Ones.l, Ones.s, NAs.s), # random intercept
  rd.time.l2 = c(NAs.l, NAs.s, LongDat$time, SurvNew$time, NAs.s), #random slope
  
  # copied terms
  id.cv1 = c(NAs.l, SurvNew$pt.id, NAs.l, NAs.s, NAs.s), 
  id.cv2 = c(NAs.l, NAs.s, NAs.l, SurvNew$pt.id, NAs.s),
  
  weight.c1 = c(NAs.l, rep(-1, n.s), NAs.l, NAs.s, NAs.s),
  weight.c2 = c(NAs.l, NAs.s, NAs.l, rep(-1, n.s), NAs.s)
)

# expectation
expectation <- list(expect = c(NAs.l, NAs.s, NAs.l, NAs.s, SurvNew$E..coxph)) 

joint.data  <- c(covariate.s, covariate.l, expectation, Surv.info$data.list)
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
                 # Longitudinal submodel 1
                 intercept.l1 + time1 * Effect1 +
                 f(id.inte.1, rd.inte.l1, model='iid2d', n=2*n.patients,
                   constr = F, hyper=rd.prior) +
                 f(id.slope.1, rd.time.l1, copy='id.inte.1') +
                 
                 # persudo submodel 1
                 f(id.cv1, weight.c1, model = "iid", # infinite percison
                   hyper = list(prec = list(initial = -6, fixed = TRUE)),
                   constr = F) +
                 
                 # Longitudinal submodel 2
                 intercept.l2 + time2 * Effect2 +
                 f(id.inte.2, rd.inte.l2, model='iid2d', n=2*n.patients,
                   constr = F, hyper=rd.prior) +
                 f(id.slope.2, rd.time.l2, copy='id.inte.2') +
                 
                 # persudo submodel 2
                 f(id.cv2, weight.c2, model = "iid", # infinite percison
                   hyper = list(prec = list(initial = -6, fixed = TRUE)),
                   constr = F) +
                 
                 #Survival submodel
                 intercept.s + trt +
                 f(basehaz, 
                   model = "rw1", values = baseline.hazard.values, 
                   hyper = basehaz.prior, constr = FALSE, 
                   diagonal = 0.01, scale.model = TRUE) +
                 
                 # copied terms
                 f(cv.1, copy = "id.cv1", 
                   hyper = associ.prior, initial = 0.1) +
                 f(cv.2, copy = "id.cv2", 
                   hyper = associ.prior, initial = 0.1),
               
               data = joint.data, 
               family = c("gaussian", "gaussian", 
                          "gaussian", "gaussian",
                          "poisson"),
               E = joint.data$expect,
               control.compute = list(dic = TRUE,
                                      waic = TRUE,
                                      cpo = TRUE,
                                      config = TRUE),
               control.fixed = fixed.prior,
               control.inla=list(int.strategy="eb"))


summary(jm.fit)
