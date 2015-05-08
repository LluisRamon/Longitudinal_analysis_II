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

TODO(Mathieu)

# Missingnes --------------------------------------------------------------

# Sensitivity analysis
# Pattern

# Try to impute




