# Load packages -----------------------------------------------------------
library("nlme")
# version 3.1-120 is required
library("dplyr")
library("ggplot2")
library("tidyr")
library("corrplot")
library("grid")
library("gridExtra")

# Import dataset ----------------------------------------------------------
cows <- read.table("data/cattle_mes dades.txt", header = TRUE, 
                   sep = "\t", dec = ",", na.strings = "")

names(cows) <- c("id", "dose", "pcv", "time", "nbirth")
cows$dose <- factor(cows$dose, levels = c("L", "M", "H"))
str(cows)
head(cows)
summary(cows)
cows$pcv.b <- as.numeric(cows$pcv > 20)
cows$pcv.f <- factor(cows$pcv > 20, labels = c("Unhealthy", "Healthy"))
cows$time.f <- paste("Time", cows$time)
cows.com <- na.omit(cows)
cows.com$idDose <- as.factor(paste(cows.com$id, cows.com$dose, sep = "_"))

table(cows$pcv.f , cows$time.f, useNA = "ifany")
table(cows.com$pcv.f , cows.com$time.f)

qplot(factor(time), pcv.b, data = cows.com, group = id, geom = c("line", "point"), facets = id~ dose, 
      colour = factor(id)) + scale_color_discrete(guide = 'none')

# Covariance and Correlation Structure ------------------------------------

cows.w <- spread(cows, time, pcv)
cor(cows.w[, c("1", "2", "3")], use = "pairwise.complete.obs")

# GEE ---------------------------------------------------------------------

library("geepack")

# In help it is pointed in id the following
# "Data are assumed to be sorted so that observations 
# on a cluster are contiguous rows for all entities in the formula."
cows.com <- cows.com %>% arrange(idDose)

# model1 <- geeglm(pcv.b ~ time, id = idDose, data = cows.com,
#                family = binomial, corstr = "unstructured", scale.fix = TRUE)
# 
# model12 <- geeglm(pcv.b ~ dose*time + nbirth, id = idDose, data = cows.com,
#                   family = binomial, corstr = "unstructured", scale.fix = TRUE)
# Not converging with the unstructured

model00 <- geeglm(pcv.b ~ dose, id = idDose, data = cows.com,
                 family = binomial, corstr = "exch", scale.fix = TRUE)

summary(model00)

model01 <- update(model00, formula = ~. + time)
model02 <- update(model00, formula = ~. + nbirth)

summary(model01)
summary(model02)

# Criteria for comparision -> QIC
library("MuMIn") # Package for GEE model selection (QIC)
model.sel(model00, model01, model02, rank = QIC)

anova(model00, model01)

# Iteration 1
model11 <- update(model01, formula = ~. + dose:time)
model12 <- update(model01, formula = ~. + nbirth)

summary(model11)
summary(model12)

anova(model01, model11)
anova(model01, model12)

# Model 01 is the selected one

# Choose best working correlation structure

model01ind <- update(model01, corstr = "independence")
model01ar1 <- update(model01, corstr = "ar1")
# model01uns <- update(model01, corstr = "unstructured") 
# Error en `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
#   contrasts can be applied only to factors with 2 or more levels

model.sel(model01, model01ar1, model01ind, rank = QIC)

# All models quite the same
summary(model01)
summary(model01ind)

# The independence model is the one with better QIC. 
# TODO(Lluis): Think about interpretation.
# Independence model smaller variance in non significative coefficient doseM
# Exchangability smaller variance in significative coefficnets intercept, doseH and time

# GLMM --------------------------------------------------------------------

# TODO(Mathieu)

library("lme4")

modelmm1 <- glmer(pcv.b~dose+time+(0+dose+time|id),data=cows.com,family=binomial)
summary(modelmm1)

modelmm2 <- glmer(pcv.b~dose+time+(0+time|idDose),data=cows.com,family=binomial)
summary(modelmm2)

modelmm3 <- glmer(pcv.b~dose*time+(0+time|id/dose),data=cows.com,family=binomial)
summary(modelmm3)
# Warning + non significant

# A lot of change in the st. dev depending on the group used.

model.mm.1 <- glmer(pcv.b~dose*time+nbirth+(time|id/dose),data=cows.com,family=binomial)
summary(model.mm.1)
# cor -1

model.mm.2 <- glmer(pcv.b~dose*time+nbirth+(0+time|id/dose),data=cows.com,family=binomial)
summary(model.mm.2)
# huge difference between the st dev of fixed effects and the one for random effects.
# Is it a good sign when it comes to the importance of random effects?
# every covariate is significant.
# Meaning of the estimates?


# Is the model stable?
modelmm1 <- glmer(pcv.b~dose*time+nbirth+(0+time|id/dose),data=cows.com,family=binomial)
summary(modelmm1)

modelmm12 <- glmer(pcv.b~dose+(0+time|id/dose),data=cows.com,family=binomial)
summary(modelmm12)

modelmm13 <- glmer(pcv.b~time+(0+time|id/dose),data=cows.com,family=binomial)
summary(modelmm13)

modelmm14 <- glmer(pcv.b~nbirth+(0+time|id/dose),data=cows.com,family=binomial)
summary(modelmm14)
# Correlation between covariates.

# Missingnes --------------------------------------------------------------

# Sensitivity analysis
# Pattern

# Try to impute




