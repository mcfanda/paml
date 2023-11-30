source("R/constants.R")
source("R/functions.R")
source("R/operations.R")
source("R/Design.R")






d<-s_data()
d<-d + s_integer("x")
d<-d + s_factor("f1",levels=3)
d<-d + s_dependent("y")
d<-d + s_params(c(1,.5))
d<-d + s_model(y~.b(1) * 1+.b(2)*x+ .b(3,4)*f1)
d$print()
dd<-d$get_sample(N = 100)

m<-d$root_model
stats:::deviance.lm

m$residuals<-(m$residuals/sd(m$residuals))
deviance(m)
anova(m)
stats:::simulate.lm
yy<-simulate(m,1)[[1]]
hist(dd$x)
summary(lm(yy~x+f1,dd))
lm(y~x,dd)
m2<-glm(f1~x,family = binomial(),data=dd)
summary(m2)
bb<-binomial()
bb$simulate
deviance(m2)
fitted(m2)



b0 = .5
b1 = 2.5
b2 = 5

set.seed(16)
x = runif(100, min = 0, max = 1)
head(x) # First six values of x

lambda = exp(b0 + b1*x )
head(lambda)

y = rpois(100, lambda = lambda)
m<-glm(y~x,family = poisson())
performance::r2(m)

library(simglm)
library(tidyr)
sim_arguments <- list(
  formula = y ~ 1 +x ,
  reg_weights = c(4, 3),
  fixed = list(x = list(var_type = 'continuous', mean = 180, sd = 30)),
  sample_size = 100,
  outcome_type="binary"
)

data<-simglm(sim_arguments)

model<-  glm(y~x,family = binomial(),data=data)


library(lmerTest)
mod<-lmer(y~x+(1|cluster),data=d)
summary(mod)
table(d$x_1,d$z_1)
model.matrix(y ~ 1 + x+z,mdata)
hist(d$y)
m<-glm(y~x,family = poisson(),data=d)
deviance(m)
performance::r2(m)
