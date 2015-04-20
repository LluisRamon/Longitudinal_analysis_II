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
cows$pcv.b <- cows$pcv > 22
cows.com <- na.omit(cows)
cows.com$idDose <- paste(cows.com$id, cows.com$dose, sep = "_")

table(cows.com$pcv.b , cows.com$time)

qplot(factor(time), pcv.b, data = cows.com, group = id, geom = c("line", "point"), facets = id~ dose, 
      colour = factor(id)) + scale_color_discrete(guide = 'none')

# Covariance and Correlation Structure ------------------------------------

cows.w <- spread(cows, time, pcv)
cor(cows.w[, c("1", "2", "3")], use = "pairwise.complete.obs")

# GEE ---------------------------------------------------------------------

library("geepack")

?geese

model <- geese(pcv.b ~ dose*time + nbirth, id = idDose, data = cows.com,
                family = binomial, corstr = "exch", scale.fix = TRUE)

model <- geeglm(pcv.b ~ dose*time + nbirth, id = idDose, data = cows.com,
               family = binomial, corstr = "exch", scale.fix = TRUE)

model2 <- geeglm(pcv.b ~ dose*time, id = idDose, data = cows.com,
               family = binomial, corstr = "exch", scale.fix = TRUE)

# nonstructure (weaker)

summary(model)

anova(model, model2)
# Pseudo likelihood function

library("gee")

model <- geese(pcv.b ~ dose*time + nbirth, id = idDose, data = cows.com,
               family = binomial, corstr = "exch", scale.fix = TRUE)


# Missingnes --------------------------------------------------------------

# Sensitivity analysis
# Pattern

# Try to impute




