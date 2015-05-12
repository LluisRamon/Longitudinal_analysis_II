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
cows$pcv.f <- factor(cows$pcv > 20, labels = c("Healthy", "Unhealthy"))
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

?geese
?geeglm

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
model01 <- geeglm(pcv.b ~ dose + time, id = idDose, data = cows.com,
                 family = binomial, corstr = "exch", scale.fix = TRUE)
model02 <- geeglm(pcv.b ~ dose + nbirth, id = idDose, data = cows.com,
                 family = binomial, corstr = "exch", scale.fix = TRUE)
model03 <- geeglm(pcv.b ~ dose + time:dose, id = idDose, data = cows.com,
                 family = binomial, corstr = "exch", scale.fix = TRUE)

summary(model00)
summary(model01)
summary(model02)
summary(model03)

# Criteria for comparision -> QIC
# Scales changes a lot from one pseudo-likelihood from one to other

library("MuMIn") # Package for GEE model selection (QIC)
model.sel(model00, model01, model02, model03, rank = QIC)

anova(model00, model01)

# Iteration 1
model11 <- geeglm(pcv.b ~ dose*time, id = idDose, data = cows.com,
                  family = binomial, corstr = "exch", scale.fix = TRUE)

model12 <- geeglm(pcv.b ~ time + dose + nbirth, id = idDose, data = cows.com,
               family = binomial, corstr = "exch", scale.fix = TRUE)

model13 <- geeglm(pcv.b ~ time + dose + dose:nbirth, id = idDose, data = cows.com,
                  family = binomial, corstr = "exch", scale.fix = TRUE)

# nonstructure (weaker)

summary(model11)
summary(model12)
summary(model13)


anova(model1, model2)
# Pseudo likelihood function

library("gee")

model <- geese(pcv.b ~ dose*time + nbirth, id = idDose, data = cows.com,
               family = binomial, corstr = "exch", scale.fix = TRUE)


# GLMM --------------------------------------------------------------------

library("lme4")

# Random effect on the intercep and the slope.

model.re.11 <- glmer(pcv.b~dose+(time|id/dose),data=cows.com,family=binomial,control=glmerControl(optimizer="bobyqa"))
summary(model.re.11)
# Correlation of -1.

start11 <- unlist(getME(model.re.11,name="ST"))
start12 <- unlist(getME(model.re.11,name="theta"))
start13 <- unlist(getME(model.re.11,name="beta"))
start1 <- list(c(start11,start12,start13)) 

model.re.12 <- glmer(pcv.b~time+(time|id/dose),data=cows.com,family=binomial,start=start1,control=glmerControl(optimizer="bobyqa"))
summary(model.re.11)
# The model failed to converge and corr of -1.

# We will not use a random intercept.


# Random effect only on the slope.

model.res.11 <- glmer(pcv.b~dose+(0+time|id/dose),data=cows.com,family=binomial)
summary(model.res.11)
# AIC = 57.9

start11 <- unlist(getME(model.res.11,name="ST"))
start12 <- unlist(getME(model.res.11,name="theta"))
start13 <- unlist(getME(model.res.11,name="beta"))
start1 <- list(c(start11,start12,start13))

model.res.12 <- glmer(pcv.b~time+(0+time|id/dose),data=cows.com,family=binomial)
summary(model.res.12)
# AIC = 63.7

model.res.13 <- glmer(pcv.b~nbirth+(0+time|id/dose),data=cows.com,family=binomial)
summary(model.res.13)
# AIC = 65.9

# We use the model with dose.
model.res.21 <- glmer(pcv.b~dose+time+(0+time|id/dose),data=cows.com,family=binomial,start=start1,control=glmerControl(optimizer="bobyqa"))
summary(model.res.21)
# AIC = 36.59

anova(model.res.11,model.res.21)

model.res.22 <- glmer(pcv.b~dose+nbirth+(0+time|id/dose),data=cows.com,family=binomial,start=start1,control=glmerControl(optimizer="bobyqa"))
summary(model.res.22)
# AIC = 59.9

start21 <- unlist(getME(model.res.21,name="ST"))
start22 <- unlist(getME(model.res.21,name="theta"))
start23 <- unlist(getME(model.res.21,name="beta"))
start2 <- list(c(start21,start22,start23)) 

# We use the model with dose+time.
model.res.31 <- glmer(pcv.b~dose+time+nbirth+(0+time|id/dose),data=cows.com,family=binomial,start=start2,control=glmerControl(optimizer="bobyqa"))
summary(model.res.31)
# AIC = 38.3

anova(model.res.21,model.res.31)
# model with dose+time is better.

model.res.41 <- glmer(pcv.b~dose*time+(0+time|id/dose),data=cows.com,family=binomial,start=start2,control=glmerControl(optimizer="bobyqa"))
summary(model.res.41)
# AIC = 39

anova(model.res.21,model.res.41)
# We keep the model with dose+time.

finalModel <- model.re.21
summary(finalModel)

orfixed <- c("intercept"=exp(coef(finalModel)$id[[1]][1]),"doseM"=exp(coef(finalModel)$id[[2]][1]),
             "doseH"=exp(coef(finalModel)$id[[3]][1]))

orrandom <- exp(c(coef(finalModel)$`dose:id`[[4]],coef(finalModel)$id[[4]]))

# Is it ok to have Odds Ratio either very large or very small (*10^16;*10^-37)?

library("caret")

pred.bin <- function(model){
  pred <- predict(model)
  p <- exp(pred)/(1 + exp(pred))
  round(p)
}

prediction <- pred.bin(finalModel)

confusionMatrix(data=prediction,reference=cows.com$pcv.b)

# Is it ok to have a sensitivity and a specificity of 1?


# Missingnes --------------------------------------------------------------

# Sensitivity analysis
# Pattern

# Try to impute




