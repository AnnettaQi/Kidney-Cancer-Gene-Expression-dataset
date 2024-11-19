#Download a dataset from bioconductor
source("http://bioconductor.org/biocLite.R")
biocLite("kidpack") #kidney microarray data
data(eset)
#Subtypes of kidney cancer has three levels
subtype <- as.factor(eset$type) #Change it to be a level variable
summary(subtype) #52 ccRCC, 9 chRCC, 13 pRCC, total 74
#I want to do a mulitple linear regression to check if age, gender and subtypes
will affect patient's survival time.
#Step1:EDA for survival time
survival <- eset$survival.time
hist(survival,xlab="Survival Time",ylab="Frequency of Survival
Time",ylim=c(0,25),col="green")
boxplot(survival,main="Boxplot of survival time for 74 patients with kidney
cancer",ylab="Survival Time",col="grey")
qqnorm(survival,main="qqplot for survival time")
qqline(survival,col="red")
#Step2:EDA for age
age <-phenoData(eset)$age
summary(age) #Ages ranges from 26 to 85 with mean 59.77
boxplot(age,main="Boxplot of age for 74 patients with kidney
cancer",ylab="Age",col="grey") #Summetric, no skewness
#Step3:EDA for metastasis
m <- phenoData(eset)$m
summary(m) #27 NA's
table(m) #24 without metastasis, 23 with metastasis
#Step4:EDA for gender
gender <-ifelse(phenoData(eset)$sex == 'm',1,0)
table(gender)
#Plot of survival time versus age under three subtypes
plot(age,survival,type='n',xlab="Age",ylab="Survival
Time",xlim=c(0,90),main="Plot of survival time for different ages under three
cancer subtypes",cex.lab=1,cex.axix=1)
points(age[subtype=="ccRCC"],survival[subtype=="ccRCC"],col="green",pch=19,
       xlim=c(0,90))
points(age[subtype=="pRCC"],survival[subtype=="pRCC"],col="red",pch=19)
points(age[subtype=="chRCC"],survival[subtype=="chRCC"],col="blue",pch=19)
legend("topleft",c("ccRCC","pRCC","chRCC"),fill=c("green","red",'blue'))
#Statistical modelling: multiple linear regression
ccRCC <- ifelse(subtype =="ccRCC",1,0)
pRCC <- ifelse(subtype =="pRCC",1,0)
chRCC <- ifelse(subtype =="chRCC",1,0)
gender <-ifelse(phenoData(eset)$sex == 'm',1,0)
lmod <- lm(survival~age+gender+ccRCC+pRCC+m)
summary(lmod)
plot(lmod)
lmod2 <- lm(survival~m)
summary(lmod2)
plot(m[which((!is.na(m)) & (!is.na(survival)))],survival[which((!is.na(m)) &
                                                                 (!is.na(survival)))],ylab="Survival Time",xlab="Metastasis",main="Plot of
Survival Time Versis Metastasis")
points(m[which((!is.na(m)) & (!is.na(survival)))],lmod2$fitted.values,col="red")
#Predicted values are very close to 13 when m=1, 29 when m=0
cor(m[which((!is.na(m)) & (!is.na(survival)))],age[which((!is.na(m)) &
                                                           (!is.na(survival)))]) #-0.058 no correlation
#Fitting proportional hazard model for survival time as a function of age, gender,
subtypes and metastasis
library(survival)
death <- phenoData(eset)$died
phmod <-
  coxph(Surv(time=survival,event=death)~age+gender+ccRCC+pRCC+m)
summary(phmod)