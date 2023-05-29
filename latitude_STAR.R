library(tidyverse)




#Load data

df <- read.csv("Data/lat-feed-data.csv")%>%
  mutate(across(c(site, chamber, bath, pond, id), as.factor))


#Calculate background mortality with binomial glmer


library(lme4)

##Create a subset of the data with just control treatments. Remove the bath over 35C because that appears to be outlier
controlDf <- df%>%
  filter(pond == "Control")%>%
  filter(temp < 35)

#Run glmer. Remove the water bath above 35C because it had very high background mortality of prey     
backMortGlmer <- controlDf%>%
  glmer(data = ., cbind(surv, dead) ~ temp  + (1|site) + (1|bath), family = "binomial")


##Obtain predicted values of background mortality
library(bootpredictlme4)

controlFit <- controlDf%>%
  predict(backMortGlmer, newdata = ., re.form = NA, se.fit = T, nsim = 1000, type = "response")

controlFit2 <- controlDf%>%
  mutate(fit = unname(data.frame(controlFit$fit)$controlFit.fit),
         lwr = unname(t(data.frame(controlFit$ci.fit))[,1]),
         upr = unname(t(data.frame(controlFit$ci.fit))[,2]),
         backMort = 100 - fit*100)

##Plot with data
controlFit2%>%
  ggplot(aes(x = temp, y = dead)) +
  geom_line(aes(y = 100 - fit*100 ), linewidth = 2) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymax = 100 - upr*100, ymin = 100 - lwr*100), alpha = 0.25, linetype = 0) +
  theme_classic() +
  geom_point(data = filter(df, pond == "Control", temp >35), shape = 1, size = 3) + #Include the outlier point that was excluded from analysis
  labs(y = expression(paste(italic("Daphnia pulex"), " background mortality ± CI")),
       x = "Temperature (C)") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size=16))

#Uncomment to save   
#ggsave("Figures/backMort.jpg")

##Create a final data set 
feedRateDf <- df%>%
  filter(temp < 35)%>%
  left_join(dplyr::select(controlFit2, site, chamber, bath, backMort), by = c("site", "chamber", "bath"))%>%
  mutate(numEaten = dead - backMort)%>%
  mutate(numEaten = replace(numEaten, numEaten < 0, 0))%>%
  filter(pond != "Control")





#STAR method of analysis
library(rSTAR)

#Make model matrix
X <- model.matrix(numEaten ~ temp*site + I(temp^2)*site + scale(predmass), data = feedRateDf)

#Dimensions
n = nrow(X); p = ncol(X)

#Define the estimator function
estimator = function(numEaten) lm(numEaten ~ X - 1)

#Select a transformation
transformation = 'log'

#EM algorithm for STAR

fit = star_EM(y = numEaten,
                   estimator = estimator,
                   transformation = transformation)

y = round(feedRateDf$numEaten)
#x = feedRateDf$temp


#Function to calculate p values for each predictor
get_p <- function(predictor_column){
  j = predictor_column #corresponds to the column number of the effect to be tested
  fit_0 = star_EM(y = y,
                        estimator = function(y) lm(y ~ X[,-j] - 1),
                        transformation = transformation)
  
  p_val = pchisq(-2*(fit_0$logLik - fit$logLik),
                 df = length(predictor_column), lower.tail = FALSE)
  
  return(p_val)
}

#Get p values


get_p(c(2,5,7:10)) #temp
get_p(c(5,9:10)) #temp^2
get_p(c(3,4,7:10)) #site
get_p(7:10) #site:temp
get_p(9:10) #site:temp^2
get_p(6) #predmass







# Compute the predicted values. In order to get separate curves for particular sites, I'll add an out-of-sample matrix that only has one site and  repeat for each one
##Matrix of MI values

MI.X <- feedRateDf%>%
  filter(site == "MI")%>%
  model.matrix(numEaten ~ temp*site + I(temp^2) + scale(predmass), data = .)

y.pred.MI = star_pred_dist(y, X, MI.X,  transformation = 'log')


# Using these draws, compute prediction intervals for STAR:
MI.PI = t(apply(y.pred.MI, 2, quantile, c(0.05, 1 - 0.05)))


# Plot the results: PIs and CIs
plot(MIdf$temp, MIdf$numEaten, ylim = range(MIdf$numEaten, y.pred.MI), main = 'STAR: 95% prediction Intervals')
lines(fit$fitted.values)
lines(MIdf$temp, MI.PI[,1], col='darkgray', type='s', lwd=4);
lines(MIdf$temp, MI.PI[,2], col='darkgray', type='s', lwd=4)





