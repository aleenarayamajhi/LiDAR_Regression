###INSTALLING MULTIPLE PACKAGES
pkgs <- c('car', 'MASS', 'leaps', 'stats', 'locfits', 'Hmisc', 'qpcR')
#install.packages(pkgs, dependencies=TRUE)

###### Multiple linear regression diagnostics

#set up the data (read the data)
#install.packages("readxl")
library(readxl)
demo <- read_excel("C:/Spring_2024/LiDAR Remote Sensing/Assignments/Assignment_4_Regression/LiDAR_Regression/Data/Data_LAI-excel.xlsx")
#demo <- read.table("C:/Spring_2024/LiDAR Remote Sensing/Assignments/Assignment_4_Regression/LiDAR_Regression/Data/Data_vol.txt", header = TRUE)
print(demo)
names(demo) #it shows a list of all variables in your table
head(demo) #you can use this to show the first 5 observations
tail(demo) #you can use this to show the last 5 observations

#####Add BA and Biomass to the table with lidar metrics
attach(demo)
demo$BA <- demo$DBH_mean*demo$DBH_mean*0.005454
demo <- data.frame(demo)
names(demo)

######Clean the table and keep only variables that can be used in the model 
#newdemo <- demo[, c(21:53, 76:79, 180)]
newdemo <- demo[, c(4, 37, 27, 9, 8, 12, 22)]
names(newdemo)

###first we plot the observations with respect to our Y (dependent variable)
library(carData)
library(car)
scatterplotMatrix(~ BA + all_90+ veg_stdv + veg_mean + veg_20  + all_20, data=newdemo)

#Fit the model with all variables
model <- lm(BA ~ all_90+ veg_stdv + veg_mean + veg_20  + all_20, data=newdemo)
summary(model)

#graph of the residuals vs the predicted variables
model.res = rstandard(model)
apre = fitted.values(model)
model.pre = (apre-mean(apre))/sd(apre)
plot(model.pre, model.res,ylab="Res Std", xlab="Pre Std") 
abline(0, 0)   

# distribution of studentized residuals
library(MASS)
sresid <- studres(model) 
hist(sresid, freq=FALSE, 
   main="Distribution of Studentized Residuals")
xmodel<-seq(min(sresid),max(sresid),length=40) 
ymodel<-dnorm(xmodel) 
lines(xmodel, ymodel) 

#graph of the residuals
qqnorm(model.res,ylab="Res Sdt", xlab="Normal Scores") 
abline(0,1)

#Pearsons correlation coefficients
cor(newdemo)

library("PerformanceAnalytics")

# Exclude DBH_mean and BA columns from corrtable
corrtable <- newdemo[, c(2, 3, 4, 5, 6, 7)]

# Generate correlation chart
chart.Correlation(corrtable, histogram = TRUE)

# CORRELATION PANEL
panel.cor <- function(x, y) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- round(cor(x, y), digits = 2)
  txt <- paste0("R = ", r)
  cex.cor <- 0.8 / strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

####CUSTOMIZE UPPER PANEL
upper.panel<-function(x, y){
  points(x,y, pch = 19, cex = 0.3)
}


###CREATE THE PLOTS
pairs(newdemo[, c(2, 3, 4, 5, 6, 7)], 
      lower.panel = panel.cor,
      upper.panel = upper.panel)


# Other useful functions 
coefficients(model) # model coefficients
confint(model, level=0.95) # CIs for model parameters 
fitted(model) # predicted values
residuals(model) # residuals
anova(model) # anova table  

#removed veg_stdv and all_20 based on correlation values
# Stepwise Regression
library(MASS)
full.model <- lm(BA ~ all_90 + veg_mean + veg_20, data=newdemo)
step1 <- stepAIC(full.model, direction="both")
summary(step1)
step2 <- stepAIC(full.model, direction="forward")
summary(step2)
step3 <- stepAIC(full.model, direction="backward")
summary(step3)
# display results
step1$anova 
step2$anova 
step3$anova  

# All Subsets Regression
library(leaps)
leaps<-regsubsets(BA ~ all_90 + veg_mean + veg_20, data=newdemo,nbest=5, method="exhaustive")
summary(leaps)
# view results 
best_summary=summary(leaps)
names(best_summary)
par(mfrow = c(2, 2))
plot(best_summary$rsq, xlab = "Number of Variables", ylab = "R2")
plot(best_summary$rss, xlab = "Number of Variables", ylab = "RSS")
plot(best_summary$adjr2, xlab = "Number of Variables", ylab = "Adjusted RSq")
plot(best_summary$bic, xlab = "Number of Variables", ylab = "BIC")
plot(best_summary$cp, xlab = "Number of Variables", ylab = "Cp")

# We will now plot a red dot to indicate the model with the largest adjusted R^2 statistic.
# The which.max() function can be used to identify the location of the maximum point of a vector
adj_r2_max = which.max(best_summary$adjr2) # 11

# The points() command works like the plot() command, except that it puts points 
# on a plot that has already been created instead of creating a new plot
points(adj_r2_max, best_summary$adjr2[adj_r2_max], col ="red", cex = 2, pch = 20)

# We'll do the same for C_p and BIC, this time looking for the models with the SMALLEST statistic
plot(best_summary$cp, xlab = "Number of Variables", ylab = "Cp")
cp_min = which.min(best_summary$cp) # 10
points(cp_min, best_summary$cp[cp_min], col = "red", cex = 2, pch = 20)

plot(best_summary$bic, xlab = "Number of Variables", ylab = "BIC")
bic_min = which.min(best_summary$bic) # 6
points(bic_min, best_summary$bic[bic_min], col = "red", cex = 2, pch = 20)

summary.out <- summary(leaps)
as.data.frame(summary.out$outmat)

library(car)
layout(matrix(1:2, ncol = 2))
## Adjusted R2
res.legend <-
  subsets(leaps, statistic="adjr2", legend = FALSE, min.size = 2, main = "Adjusted R^2")
## Mallow Cp
res.legend <-
  subsets(leaps, statistic="cp", legend = FALSE, min.size = 2, main = "Mallow Cp")
abline(a = 1, b = 1, lty = 2)

# plot a table of models showing variables in each model.
# models are ordered by the selection statistic.
plot(leaps,scale="r2")
plot(leaps,scale="Cp")
plot(leaps,scale="adjr2")
# plot statistic by subset size 

# Evaluate extreme observations on the chosen model
model <- lm(BA ~ veg_mean + all_90 + veg_mean + veg_20, data=newdemo)
summary(model)
car::outlierTest(model) # Bonferonni p-value for most extreme observations
cooksd <- cooks.distance(model) # Cook's Distance
# Set up margin values
par(mar=c(5, 5, 2, 2))  # Adjust margins as needed



plot(cooksd, pch="*", cex=2, main="Influential Obs by Cooks distance")  # plot cook's distance
abline(h = 4*mean(cooksd, na.rm=T), col="red")  # add cutoff line
text(x=1:length(cooksd)+1, y=cooksd, labels=ifelse(cooksd>4*mean(cooksd, na.rm=T),names(cooksd),""), col="red")  # add labels
###Outlier with respect to X - Leverage
plot(hatvalues(model), pch=23, bg='blue', cex=0.5, ylab='Hat values') ###Outlier with respect to Y

# Infuential points
# added variable plots 
avPlots(model)
# Cook's D plot
# identify D values > 4/(n-k-1) 
cutoff <- 4/((nrow(newdemo)-length(model$coefficients)-2)) 
plot(model, which=4, cook.levels=cutoff)
# Influence Plot Leverage
influencePlot(model, scale=10, xlab="Hat-Values", ylab="Studentized Residuals")

###DFBETAS
#find number of observations
n <- nrow(newdemo)
#calculate DFBETAS threshold value
thresh <- 2/sqrt(n)
thresh
#calculate DFBETAS for each observation in the model
dfbetas <- as.data.frame(dfbetas(model))
dfbetas
#specify 2 rows and 1 column in plotting region
par(mfrow=c(2,1))
#plot DFBETAS for a variable with threshold lines
plot(dfbetas$veg_mean, type='h')
abline(h = thresh, lty = 2)
abline(h = -thresh, lty = 2)
#plot DFBETAS for hp with threshold lines 
plot(dfbetas$veg_20, type='h')
abline(h = thresh, lty = 2)
abline(h = -thresh, lty = 2)

###DFFITS
#calculate DFFITS for each observation in the model
dffits <- as.data.frame(dffits(model))
dffits
#find number of predictors in model
p <- length(model$coefficients)-1
#find number of observations
n <- nrow(newdemo)
#calculate DFFITS threshold value
thresh <- 2*sqrt(p/n)
thresh
#sort observations by DFFITS, descending
dffits[order(-dffits['dffits(model)']), ]
#plot DFFITS values for each observation
plot(dffits(model), type = 'h')
#add horizontal lines at absolute values for threshold
abline(h = thresh, lty = 2)
abline(h = -thresh, lty = 2)

#Criteria indicators
library(stats)
AIC(model)
BIC(model)
library(locfit)
cp(model)

#Predictions
library(stats)
predict(model,se.fit=TRUE,interval="confidence",level=0.95)

#PRESS statistic (compare values among several models)
library(qpcR)
model <- lm(BA ~ veg_mean + veg_20 + all_90, data=newdemo)
PRESS(model)
barplot(model$residuals)

#Colinearity Variance Inflation Factor - VIF 
library(car)
vif(model)
sqrt(vif(model))>2  #indicates which are greater than 2

#Colinearity Condition Index - CI
library(klaR)
condition_index <- cond.index(BA ~ veg_mean+ veg_20 + all_90, data=newdemo)
condition_index

# Load required library
library(klaR)

# Compute condition index
condition_index <- cond.index(BA ~ veg_mean + veg_20 + all_90, data = newdemo)

# Compute correlation matrix
correlation_matrix <- cor(newdemo[, c("veg_mean", "veg_20", "all_90")])

# Create correlation heatmap
library(ggplot2)
library(reshape2)
melted_corr_matrix <- melt(correlation_matrix)
ggplot(data = melted_corr_matrix, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "yellow", high = "pink") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Heatmap of Predictor Variables",
       x = "Predictor Variables",
       y = "Predictor Variables",
       fill = "Correlation Coefficient")

