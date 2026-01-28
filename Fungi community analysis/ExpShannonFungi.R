# Shannon exp funghi
# Load library
library(psych)

# Load the data
FunghiExpH<-read.table("EXPShannon_funghi.csv",h=T, sep=",")
# Visualize
FunghiExpH
# Set as factor the landcover
FunghiExpH$Landcover<-as.factor(FunghiExpH$Landcover)
# Verify
summary(FunghiExpH)
# Calculate average
tapply(FunghiExpH$expH_fungi,FunghiExpH$Landcover, mean)
# Calculate STD
tapply(FunghiExpH$expH_fungi, FunghiExpH$Landcover, sd)
# Or smply
describeBy(FunghiExpH, group=FunghiExpH$Landcover)

# KW test and Dunn's test
dunn.test(FunghiExpH$expH_fungi, FunghiExpH$Landcover, method="bonferroni")


