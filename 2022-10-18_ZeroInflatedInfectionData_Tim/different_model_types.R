#!bin/bash

library(pscl)
library(quantreg)
library(MASS)
library(tidyverse)

## ZIPoisson is better when data not over-dispersed (aka when mean is > than the variance). Otherwise try ZINB

dat <- read.csv("FigDataBd5.csv") %>%
  mutate(logDaysSince=log(DaysSince)
         ,DaysFactor = factor(DaysSince)
         ,logJliv=log(JlivNum+1)
         ,Bd=round(exp(LogBd)-1))

dat$Bd
### Look at distribution

ggplot(dat) +
  geom_point(aes(x=log(JlivNum+1), y=LogBd, col=DaysFactor))
# Day 345 doesn't actually have ANY Jliv on anyone; so let's remove those.

dat_filt <- dat %>% filter(DaysSince!=345)

ggplot(dat_filt) +
  geom_point(aes(x=log(JlivNum+1), y=LogBd, col=DaysFactor))

dat_filt %>%
  ggplot() +
  geom_histogram(aes(x=LogBd, group=DaysFactor, fill=DaysFactor), position = "dodge")


# Mean is smaller than variance; go with ZINB

dat_filt %>% #filter(DaysSince<100, !is.na(Treatment)) %>%
ggplot() +
  geom_point(aes(x=(DaysSince+1), y=LogBd, col=Treatment)) +
  scale_x_log10()


##### Normal regression #####

linear_ml <- s
summary(linear_ml)
# Plotting residuals
# plot(linear_ml$fitted.values, linear_ml$residuals)
plot(linear_ml)

newdata <- data.frame(logJliv=rep(dat_filt$logJliv, 3)
                      , DaysFactor = rep(unique(dat_filt$DaysFactor), each=length(dat_filt$logJliv)))
linear_ml_predict <- predict(linear_ml, newdata = newdata,interval = "confidence")
LinearModelFit <- cbind(newdata, linear_ml_predict)

ggplot(dat_filt) +
  geom_point(aes(x=logJliv, y=LogBd, col=DaysFactor)) +
  geom_line(data=LinearModelFit, aes(x=logJliv, y=fit, col=DaysFactor)) +
  geom_ribbon(data=LinearModelFit, aes(x=logJliv, ymin=lwr, ymax=upr, fill=DaysFactor), alpha=0.25)

##### Poisson regression #####
nb_ml <- glm.nb(Bd ~ logJliv*DaysFactor, data=dat_filt)
summary(nb_ml)
# Plotting residuals
plot(nb_ml$fitted.values, nb_ml$residuals)
plot(nb_ml)

newdata <- data.frame(logJliv=rep(dat_filt$logJliv, 3)
                      , DaysFactor = rep(unique(dat_filt$DaysFactor), each=length(dat_filt$logJliv)))
nb_ml_predict <- predict(nb_ml, newdata = newdata,interval = "confidence")
NBModelFit <- cbind(newdata, fit=nb_ml_predict)

ggplot(dat_filt) +
  geom_point(aes(x=logJliv, y=LogBd, col=DaysFactor)) +
  geom_line(data=NBModelFit, aes(x=logJliv, y=fit, col=DaysFactor)) 

##### Zero inflated negative binomial #####
# logit part: time
# negbin part: Jliv
zinb_ml_daysF2 <- zeroinfl(Bd ~ logJliv*DaysFactor | DaysFactor,
                          data = dat_filt, dist = "negbin")
zinb_ml_daysF <- zeroinfl(Bd ~ logJliv+DaysFactor | DaysFactor,
               data = dat_filt, dist = "negbin")
zinb_ml <- zeroinfl(Bd ~ logJliv | DaysFactor,
                    data = dat_filt, dist = "negbin")

summary(zinb_ml_daysF2)
summary(zinb_ml_daysF)
summary(zinb_ml)

test <- summary(zeroinfl(Bd ~ -1 +logJliv*DaysFactor | DaysSince,
                           data = dat_filt, dist = "negbin"))

##### Zero inflated poisson #####
# logit part: time
# negbin part: Jliv
zip_ml_daysF <- zeroinfl(Bd ~ logJliv+DaysFactor | DaysFactor,
                          data = dat_filt, dist = "poisson", link="log")
zip_ml <- zeroinfl(Bd ~ logJliv | DaysFactor,
                    data = dat_filt, dist = "poisson", link="log")
summary(zeroinfl(Bd ~ logJliv*DaysFactor | DaysFactor,
                 data = dat_filt, dist = "poisson", link="log"))

summary(zip_ml_daysF)
summary(zip_ml)

### Quantile Regression #####
quantreg_ml <- rq(LogBd ~ logJliv, data=dat_filt, tau = 0.90)
quantreg_ml_days <- rq(LogBd ~ logJliv+DaysFactor, data=dat_filt, tau = 0.90)
quantreg_ml_days2 <- rq(LogBd ~ logJliv*DaysFactor, data=dat_filt, tau = 0.50)

summary(quantreg_ml, se="boot")
summary(quantreg_ml_days, se="boot")
summary(quantreg_ml_days2, se="boot")


plot(log(Bd+1) ~ logJliv, data = dat, pch = 16)
abline(zinb_ml_daysF$coefficients$count[1:2], col="red")
abline(zinb_ml$coefficients$count, col="purple")
abline(quantreg_ml$coefficients, col="green")
abline(quantreg_ml_days$coefficients[1:2], col="darkgreen")
