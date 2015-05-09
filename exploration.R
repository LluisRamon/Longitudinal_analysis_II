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

# Some more models
model13 <- update(model01, formula = ~. + nbirth*time)
model14 <- update(model01, formula = ~. + nbirth*dose)

summary(model13)
summary(model14)

model.sel(model01, model13, model14, rank = QIC)
anova(model01, model13)
anova(model01, model14)

# Model 01 is the selected one

# Choose best working correlation structure

model01ind <- update(model01, corstr = "independence")
model01ar1 <- update(model01, corstr = "ar1")
# model01uns <- update(model01, corstr = "unstructured") 
# Error en `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
#   contrasts can be applied only to factors with 2 or more levels

model.sel(model01, model01ar1, model01ind, rank = QIC)
# All models quite the similar

summary(model01)
summary(model01ar1)

# The ar1 model is the one with better QIC. 
# ar1 model smaller variance in dose coefficients
# Exchangability smaller variance in intercept and time

# Note: If cows.com was not sorted results in working correlation matrix were different!!!

# Classification table

pred.bin <- function(model){
  
  pred <- predict(model)
  p <- exp(pred)/(1 + exp(pred))
  round(p)
  
}

table(pred.bin(model01ar1), cows.com$pcv.b)

library("caret")
sensitivity(as.factor(pred.bin(model01ar1)), as.factor(cows.com$pcv.b))
specificity(as.factor(pred.bin(model01ar1)), as.factor(cows.com$pcv.b))

# GLMM --------------------------------------------------------------------

library("lme4")

model.mm.1 <- glmer(pcv.b~dose*time+nbirth+(time|id/dose),data=cows.com,family=binomial)
summary(model.mm.1)
# Warning: model failed to converge.

# The only model without any warnings regarding the convergence is the following one. 
# It is necessary to use idDose as a groupin factor for the random effects because id/dose
# does not work.

model.mm.2 <- glmer(pcv.b~dose+time+(0+time|idDose),data=cows.com,family=binomial)
summary(model.mm.2)

# We cannot take into account the second random effects (first by id then by dose).




# Missingnes --------------------------------------------------------------

aggregate(is.na(cows$pcv.f), list(Time=cows$time, Dose=cows$dose), sum)[,c(2,1,3)]

a <- barplot(table(is.na(cows$pcv.f), paste(cows$dose, cows$time, sep="_"))[2,c(4:9, 1:3)]/nrow(cows)*100, beside=T, space=c(0,0,0,1,0,0,1,0,0), xaxt="n", ylab="%")
abline(h=0)
axis(1, at=a[c(2,5,8),], paste("Dose", levels(cows$dose)), tick=F, line=1.5)
axis(1, at=a, rep(paste("T", 1:3, sep=""),3), tick=F, line=-0.5)
## Percentatge de missings en cada dosi i temps respecte del total de dades de la mostra.

fisher.test(table(is.na(cows$pcv.f),cows$dose))
fisher.test(table(is.na(cows$pcv.f),cows$time))

# oagg <- order((agg <- aggregate(is.na(cows$pcv.f), list(cows$time, cows$dose, cows$id), sum))[,"x"])
# agg <- agg[oagg,]
# agg <- agg[order(agg[[3]]),]

# Sensitivity analysis
# Pattern

# Try to impute



