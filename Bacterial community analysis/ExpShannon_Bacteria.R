# Shannon exp Bacteria
# Load library
library(psych)

# Load the data
BacteriaExpH<-read.table("ExpShannon_Bacteria.csv",h=T, sep=",")
# Visualize
BacteriaExpH
# Set as factor the landcover
BacteriaExpH$Landcover<-as.factor(BacteriaExpH$Landcover)
# Verify
summary(BacteriaExpH)
# Calculate average
tapply(BacteriaExpH$expShannon,BacteriaExpH$Landcover, mean)
# Calculate STD
tapply(BacteriaExpH$expShannon, BacteriaExpH$Landcover, sd)
# Or smply
describeBy(BacteriaExpH$expShannon, group=BacteriaExpH$Landcover)

# KW test and Dunn's test
dunn.test(BacteriaExpH$expShannon, BacteriaExpH$Landcover, method="bonferroni")


