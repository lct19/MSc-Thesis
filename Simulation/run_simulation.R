##########################
### Run the simulation ###
##########################
set.seed(2023) # set seed for simulaion
library(tidyverse) # ggplot2


#################################
### current value association ###
#################################
library(survival) # compared with Cox
Nsim <- 500 # repeated simulation times
direction <- 'JMs_cv.R' # direction of the simulation file


means = ci.l = ci.r = numeric(Nsim) # save simulation results
means.cox = se.cox = numeric(Nsim) # save simulation results

for (index in 1:Nsim){
  source(direction, local = TRUE) # one time simulation
  
  #association
  position <- nrow(jm.fit$summary.hyperpar)
  means[index] = jm.fit$summary.hyperpar[position, 1]
  ci.l[index] = jm.fit$summary.hyperpar[position, 3]
  ci.r[index] = jm.fit$summary.hyperpar[position, 5]
  means.cox[index] = summary(coxfit)$coefficients[2, 1]
  se.cox[index] = summary(coxfit)$coefficients[2, 3]
  print(index)
}

output_cv = data.frame(i = 1:Nsim, means, ci.l, ci.r,
                    means.cox, se.cox) # save simulation results
save(output_cv, file='output_cv.Rdata')

#################################
### current slope association ###
#################################
set.seed(2023)
library(JMbayes2) # compared with JMbayes2
Nsim <- 100 # repeated simulation times
direction <- 'JMs_cs.R' # direction of the simulation file

means = ci.l = ci.r = numeric(Nsim) # save simulation results
mean.jmb = cil.jmb = cir.jmb = numeric(Nsim) # save simulation results

for (index in 1:Nsim){
  source(direction, local = TRUE) # one time simulation
  
  # association, INLA
  position <- nrow(jm.fit$summary.hyperpar)
  means[index] = jm.fit$summary.hyperpar[position, 1]
  ci.l[index] = jm.fit$summary.hyperpar[position, 3]
  ci.r[index] = jm.fit$summary.hyperpar[position, 5]
  
  # association, JMbayes2
  mean.jmb[index] = summary(JMP)$Survival[2, 1]
  cil.jmb[index] = summary(JMP)$Survival[2, 3]
  cir.jmb[index] = summary(JMP)$Survival[2, 4]
  print(index)
}

output_cs = data.frame(i = 1:Nsim, means, ci.l, ci.r,
                        mean.jmb, cil.jmb, cir.jmb) # save simulation results
save(output_cs, file='output_cs.Rdata')

######################################
### cumulative effects association ###
######################################
set.seed(2023) # set seed for simulaion
Nsim <- 200 # repeated simulation times
direction <- 'JMs_cumulative.R' # direction of the simulation file

means = rep(NA, Nsim)
ci.l = rep(NA, Nsim)
ci.r = rep(NA, Nsim)

for (naind in 1:Nsim){
  source(direction, local = TRUE) # one time simulation
  # association, INLA
  position <- nrow(jm.fit$summary.hyperpar)
  means[naind] = jm.fit$summary.hyperpar[position, 1]
  ci.l[naind] = jm.fit$summary.hyperpar[position, 3]
  ci.r[naind] = jm.fit$summary.hyperpar[position, 5]
  print(naind)
}
output_cumu = data.frame(i = 1:Nsim, means, ci.l, ci.r) # save simulation results
save(output_cumu, file='output_cumu.Rdata')