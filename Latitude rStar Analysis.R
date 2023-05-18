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


##Obtain siteicted values of background mortality
library(bootsiteictlme4)

controlFit <- controlDf%>%
  siteict(backMortGlmer, newdata = ., re.form = NA, se.fit = T, nsim = 1000, type = "response")

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

numEaten <- round(feedRateDf$numEaten)

#Make model matrix

feedRateDf$pond <- droplevels(feedRateDf$pond)
X <- model.matrix(numEaten ~ temp*site + I(temp^2) + scale(predmass), data = feedRateDf)

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


#Function to calculate p values for each siteictor
get_p <- function(siteictor_column){
  j = siteictor_column #corresponds to the column number of the effect to be tested
  fit_0 = star_EM(y = y,
                        estimator = function(y) lm(y ~ X[,-j] - 1),
                        transformation = transformation)
  
  p_val = pchisq(-2*(fit_0$logLik - fit$logLik),
                 df = 1, lower.tail = FALSE)
  
  return(p_val)
}

#Get p values

get_p(2) #Temp
get_p(3) #MO
get_p(4) #TX
get_p(5) #Temp^2
get_p(6) #predmass
get_p(7) #temp:MO
get_p(8) #temp:TX



# p-value for Temp (likelihood ratio test):

j = 2 #temperature column
fit_0 = star_EM(y = y,
                      estimator = function(y) lm(y ~ X[,-j] - 1),
                      transformation = transformation)

p_val_temp = pchisq(-2*(fit_0$logLik - fit$logLik),
                    df = 1, lower.tail = FALSE)

p_val_temp
# p-value for Temp2 (likelihood ratio test):

j = 3 #temperature squared column
fit_0 = star_EM(y = y,
                      estimator = function(y) lm(y ~ X[,j] - 1),
                      transformation = transformation)

p_val_temp2 = pchisq(-2*(fit_0$logLik - fit$logLik),
                     df = 1, lower.tail = FALSE)

p_val_temp2

# p-value for site (likelihood ratio test):

j = 1 #temperature column
fit_0 = star_EM(y = y,
                      estimator = function(y) lm(y ~ X[,-j] - 1),
                      transformation = transformation)

p_val_temp = pchisq(-2*(fit_0$logLik - fit$logLik),
                    df = 1, lower.tail = FALSE)





# Compute the siteicted values. In order to get separate curves for particular site x comp combinations, I'll add an out-of-sample matrix that only has one combination and iteratively repeat for each combo
##Matrix of no site x no comp siteictors

PulexNN.X <- feedRateDf%>%
  filter(site == "N", Species == "Pulex")%>%
  model.matrix(numEaten ~ poly(temp, 2)*site*Species, data = .)

PulexNN <- feedRateDf%>%
  filter(site == "N", Species == "Pulex")


y = feedRateDf$numEaten
X = feedRateDf$temp

y.site = star_site_dist(y, X,  transformation = 'log')


# Using these draws, compute siteiction intervals for STAR:
Pulex.PI = t(apply(y.site, 2, quantile, c(0.05, 1 - 0.05)))


# Plot the results: PIs and CIs
plot(PulexNN$temp, PulexNN$numEaten, ylim = range(PulexNN$numEaten, Pulex.PI), main = 'STAR: 95% siteiction Intervals')
lines(fit$fitted.values)
lines(PulexNN$temp, Pulex.PI[,1], col='darkgray', type='s', lwd=4);
lines(PulexNN$temp, Pulex.PI[,2], col='darkgray', type='s', lwd=4)





