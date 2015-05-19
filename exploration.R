# Load packages -----------------------------------------------------------
library("nlme")
# version 3.1-120 is required
library("dplyr")
library("ggplot2")
library("tidyr")
library("corrplot")
library("grid")
library("gridExtra")
library("lme4")

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

qplot(factor(time), pcv.b, data = cows.com, group = id, xlab = "time", geom = c("line", "point"), facets = id~ dose, 
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

# Random effect on the intercep and the slope.

model.re.11 <- glmer(pcv.b~time+(time|id/dose),data=cows.com,family=binomial,control=glmerControl(optimizer="bobyqa"))
summary(model.re.11)
# Correlation of -1.

start11 <- unlist(getME(model.re.11,name="ST"))
start12 <- unlist(getME(model.re.11,name="theta"))
start13 <- unlist(getME(model.re.11,name="beta"))
start1 <- list(c(start11,start12,start13)) 

model.re.12 <- glmer(pcv.b~time+dose+(time|id/dose),data=cows.com,family=binomial,start=start1,control=glmerControl(optimizer="bobyqa"))
summary(model.re.12)
# The model failed to converge and corr of -1.

# We will not use a random intercept.


# Random effect only on the slope. We start with time.

model.res.11 <- glmer(pcv.b~dose+time+(0+time|id/dose),data=cows.com,family=binomial,control=glmerControl(optimizer="bobyqa"))
summary(model.res.11)
# AIC = 36.59

start11 <- unlist(getME(model.res.11,name="ST"))
start12 <- unlist(getME(model.res.11,name="theta"))
start13 <- unlist(getME(model.res.11,name="beta"))
start1 <- list(c(start11,start12,start13)) 

# We use the model with dose+time as a start.

model.res.21 <- glmer(pcv.b~dose+time+nbirth+(0+time|id/dose),data=cows.com,family=binomial,start=start1,control=glmerControl(optimizer="bobyqa"))
summary(model.res.21)

anova(model.res.11,model.res.21)
# model with dose+time is better.

model.res.31 <- glmer(pcv.b~dose*time+(0+time|id/dose),data=cows.com,family=binomial,start=start1,control=glmerControl(optimizer="bobyqa"))
summary(model.res.31)

anova(model.res.11,model.res.31)
# We keep the model with dose+time.

# No point of adding random effect since there are no other fixed effects and the sample size
# is too small so we would have too many parameters to estimate.

finalModel <- model.res.11
summary(finalModel)

se <- sqrt(diag(vcov(finalModel)))

tab <- exp(cbind(inf = fixef(finalModel) - 1.96 * se, est = fixef(finalModel), 
             sup = fixef(finalModel) + 1.96 *se))

# Missings ----------------------------------------------------------------

fisher.test(is.na(cows$pcv), cows$dose)
fisher.test(is.na(cows$pcv), cows$time)

cows$doseExt <- factor(cows$dose, labels = c("Low", "Medium", "High"))
cows$missing <- is.na(cows$pcv.f)

cows2 <- cows %>% group_by(doseExt, time) %>% dplyr::summarise(n = sum(missing))
ggplot(cows2, aes(factor(time), n, fill = doseExt)) + 
  geom_bar(position = "dodge", stat="identity") + facet_grid(.~ doseExt) + 
  guides(fill = FALSE) + xlab("time") + ylab("Number of missings")

aggregate(list(Missing=is.na(cows$pcv.f), Available=!is.na(cows$pcv.f)), list(Time=cows$time, Dose=cows$dose), sum)[,c(2,1,3:4)]


cows$doseExt2 <- cows$doseExt
cows$doseExt2[cows$missing] <- NA

ggplot(cows, aes(x = factor(time),  y = factor(id), fill = doseExt2)) + facet_grid(.~ doseExt) + 
  geom_tile(colour = "black") + guides(fill = FALSE) + xlab("time") + ylab("Cow ID") +
  scale_fill_discrete(na.value = "white")


a <- -6:6
a[7] <- 6
sums <- sapply(a, function(k) {
  m <- glm(is.na(pcv.f)~time+dose, offset=k*pcv.b, data=cows)
  sm <- exp(summary(m)$coef[-(1:2), 1])
  cbind(OR=as.matrix(sm), exp(confint(m)[-(1:2),]))
}, simplify=F)
a <- -6:6
m <- glm(is.na(pcv.f)~time+dose, data=cows)
sm <- exp(summary(m)$coef[-(1:2), 1])
sums[[7]] <- cbind(OR=as.matrix(sm), exp(confint(m)[-(1:2),]))
taula <- sapply(1:length(sums), function(x) {
  paste(round(sums[[x]][,1],3), " (", round(sums[[x]][,2],3), "-", round(sums[[x]][,3],3), ")", sep="")
})
rownames(taula) <- rownames(sums[[7]])

a <- seq(-2,2,by=0.1)
sums <- sapply(a, function(k) summary(geeglm(is.na(pcv.b) ~ dose + time, offset=k*pcv.b, id = idDose, data = cows.com,
                                             family = binomial, corstr="ar1" , scale.fix = TRUE)
)$coef[-1,-3], simplify=F)
sums[[21]] <- summary(geeglm(is.na(pcv.b) ~ dose + time, id = idDose, data = cows.com,
                             family = binomial, corstr="ar1", scale.fix = TRUE))$coef[-1,-3]
taula <- sapply(1:length(sums), function(x) {
  paste(round(sums[[x]][-3,1],3), " (", round(sums[[x]][-3,2],3), ")", sep="")
})
a[20] <- 0.1
rownames(taula) <- rownames(sums[[21]])[-1]