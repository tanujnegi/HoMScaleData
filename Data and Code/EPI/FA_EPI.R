#### EFA - EPI ####

#set working directory
setwd("F:/Research Work/Survey/Data Analysis/OBJ1/Completed/EPI")

#Invoking Psych package
library(psych)
library(GPArotation)

#Get all the survey data from the CSV file
homSurveyData <- read.csv("EPI.csv")

write.csv(describe(homSurveyData), "DescriptEPI.csv")

#Make a 'working variable' so that originalsurvey data isn't disturbed
homEFA1 <- homSurveyData #[,c(2,3,5,10,11)]  #Q10 MSA lower than 0.5; others do not load eventually

#Remove missing items
homEFA1 <- na.omit(homEFA1) 

#Check if any null values left
apply(homEFA1,2,function(x)sum(is.na(x)))

#Get correlation matrix
homCor <- cor(homEFA1)

write.csv(homCor, "EPICor.csv")

#Get KMO
homKMO <- KMO(homCor)

#Bartlett test
homBartlett <- data.frame(cortest.bartlett(homCor,nrow(homEFA1)))  #p value comes out ~0

#Get eigenvalues
evHoM <- eigen(homCor)

#Get only eigenvalues and not the vector
corMatHoM <- evHoM$values

#Collect output of EFA
pa.out <- fa(homCor, nfactors=sum(corMatHoM > 1), n.obs=nrow(homEFA1), rotate="oblimin", scores="regression", residuals=FALSE, SMC=TRUE, covar=FALSE, missing=TRUE, impute="median", min.err=0.01, max.iter=50, symmetric=TRUE, fm="ml", alpha=.05, p=.05, oblique.scores=TRUE)

#Plot the eigen-value graph
plot(corMatHoM,main = "Factors Extracted",ylab = "Eigen value",xlab="No.of factors")
lines(corMatHoM,col="blue")
abline(h=1,col="red")

#Draw factor diagram
fa.diagram(pa.out,col="blue",digits = 3)

print(pa.out$loadings,sort=T, cutoff=0.3, digits=3)

#Write all the loadings in a CSV file
write.csv(pa.out$loadings, "EFALoadings.csv")

#### CFA Ord 1 - EPI ####

#invoke relevant libraries
library(lavaan) #model estimation
library(semPlot) #sem diagram
library(semTools) #reliability

#Import the dataset
homSurveyData <- read.csv("EPI.csv")

#Model specification: ML1 -> taken from EFA
cfamodel <- '
ML1 =~ Q5 + Q8 + Q12 + Q35 + Q46
'

#Connecting data and model
fit <- cfa(cfamodel, data = homSurveyData)

#Plot the model
semPlot::semPaths(fit,title=FALSE,'std', style="lisrel", rotation = 2, color = "green",
                  sizeMan = 5,sizeLats = 8,residuals = T,edge.label.cex = 1)

#See the summary of fit measures; pay special attention to significance of paths
summary(fit, fit.measures=TRUE, standardized = TRUE)

#select useful fitness measures
cfaoutput<-fitmeasures(fit,fit.measures=c("npar","chisq","df","pvalue","cfi","rmsea","tli", "srmr", "nnfi"))

#Make it more presentable
cfaoutput1 <- data.frame(cfaoutput)
round(cfaoutput1, 3)

#calculate CMIN/df
cminByDf <- (cfaoutput1[2,])/(cfaoutput1[3,])

#Adding CMIN/Df to table
cfaoutput2 <- rbind(cfaoutput1,cminByDf)

#rename the row
row.names(cfaoutput2)[10] <- "CMinByDf"

# specify reference ranges
threshold <- c("NIL", "Lower the better", "> 1", "> 0.05", "> 0.90", "< 0.08", "> 0.90", "< 0.05", "> 0.90", "< 5")
data.frame(cfaoutput2, threshold) #binding

#Get modification indicies
modindices(fit, sort = TRUE, minimum.value = 4)

#Get Reliability scores
relval <- semTools::reliability(fit)

#Transform this table into something more readable
relval <- t(relval)
relval <- relval[-8,c(1,2,5)]
colnames(relval) = c("Alpha", "CR", "AVE")

#Moving to establish discriminant validity

#Method 1: FL criterion
#avsqrt <- sqrt(relval[,3])  #Take squareroot of AVE
#lcl <- inspect(fit,"cor.lv")
#diag(lcl) <- c(avsqrt[1],avsqrt[2])  #the diagonal values should be highest amongst column
#avsqrt is greater than lcl >> good discriminant validity