
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


##Look at distribution of data
feedRateDf%>%
  ggplot(aes(x = numEaten)) +
  geom_histogram() +
  theme_classic()


#Now to run the glm's on the predator data. I'll try separate models for the 3 sites. First create subset data for each site, then create a function to run all possible models to compare AIC's within site for model selection
#Subset data

MIdata <- feedRateDf%>%
  filter(site == "MI", round(numEaten) > 0, pond != "Control", temp < 35)

MOdata <- feedRateDf%>%
  filter(site == "MO", round(numEaten) > 0, pond != "Control", temp < 35)

TXdata <- feedRateDf%>%
  filter(site == "TX", round(numEaten) > 0, pond != "Control", temp < 35)

library(MuMIn)

selectmodel<- function(stateDf){
  
  ##Create each of the models
  model1 <- glmer(cbind(round(numEaten), round(100-numEaten)) ~ poly(temp, 2) + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = stateDf,
                  control = glmerControl(optimizer = "bobyqa"))
  model2 <- glmer(cbind(round(numEaten), round(100-numEaten)) ~ poly(temp, 2) + (1|bath) + scale(predmass), family = "binomial", data = stateDf,
                  control = glmerControl(optimizer = "bobyqa"))
  model3 <- glmer(cbind(round(numEaten), round(100-numEaten)) ~ temp + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = stateDf,
                  control = glmerControl(optimizer = "bobyqa"))
  model4 <- glmer(cbind(round(numEaten), round(100-numEaten)) ~ temp + (1|bath) + scale(predmass), family = "binomial", data = stateDf,
                  control = glmerControl(optimizer = "bobyqa"))
  model5 <- glmer(cbind(round(numEaten), round(100-numEaten)) ~ poly(temp, 2) + (1|pond) + (1|bath), family = "binomial", data = stateDf,
                  control = glmerControl(optimizer = "bobyqa"))
  model6 <- glmer(cbind(round(numEaten), round(100-numEaten)) ~ poly(temp, 2) + (1|bath), family = "binomial", data = stateDf,
                  control = glmerControl(optimizer = "bobyqa"))
  model7 <- glmer(cbind(round(numEaten), round(100-numEaten)) ~ temp + (1|pond) + (1|bath), family = "binomial", data = stateDf,
                  control = glmerControl(optimizer = "bobyqa"))
  model8 <- glmer(cbind(round(numEaten), round(100-numEaten)) ~ temp + (1|bath), family = "binomial", data = stateDf,
                  control = glmerControl(optimizer = "bobyqa"))
  
  
  data.frame(model = c("model1", "model2", "model3", "model4", "model5", "model6", "model7", "model8"),
             AICc = c(AICc(model1), AICc(model2), AICc(model3), AICc(model4),
                     AICc(model5), AICc(model6), AICc(model7), AICc(model8)))%>%
    arrange(AICc)%>%
    mutate(delta = AICc - .[1,2])%>%
    return()
  
  
}

##Run the function to compare models
selectmodel(MIdata)

#MI models 3 and 1 have ∆AICc of <2. Running a likelihood ratio test

library(lmtest)

lrtest(glmer(cbind(round(numEaten), round(100-numEaten)) ~ temp + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = MIdata,
               control = glmerControl(optimizer = "bobyqa")),
         glmer(cbind(round(numEaten), round(100-numEaten)) ~ poly(temp, 2) + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = MIdata,
               control = glmerControl(optimizer = "bobyqa")))


selectmodel(MOdata)

selectmodel(TXdata)

##Six of the texas models are giving singularity warnings. Running each one here to figure out which ones

#Model 1
TXdata%>%
  glmer(cbind(round(numEaten), round(100-numEaten)) ~ poly(temp, 2) + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))%>%str()

#Model 2
TXdata%>%
  glmer(cbind(round(numEaten), round(100-numEaten)) ~ poly(temp, 2) + (1|bath) + scale(predmass), family = "binomial", data =.,
        control = glmerControl(optimizer = "bobyqa"))

#Model 3
model3 <- TXdata%>%
  glmer(cbind(round(numEaten), round(100-numEaten)) ~ temp + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Model 4
TXdata%>%
  glmer(cbind(round(numEaten), round(100-numEaten)) ~ temp + (1|bath) + scale(predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Model 5
TXdata%>%
  glmer(cbind(round(numEaten), round(100-numEaten)) ~ poly(temp, 2) + (1|pond) + (1|bath), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Model 6
feedRateDf%>%
  glmer(cbind(round(numEaten), round(100-numEaten)) ~ poly(temp, 2) + (1|bath), family = "binomial", data =.,
        control = glmerControl(optimizer = "bobyqa"))

#Model 7
model7 <- TXdata%>%
  glmer(cbind(round(numEaten), round(100-numEaten)) ~ temp + (1|pond) + (1|bath), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Model 8
TXdata%>%
  glmer(cbind(round(numEaten), round(100-numEaten)) ~ temp + (1|bath), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Models 3 and 7 are the ones that are converging. 7 had the lowest AIC of the eight models, and 3 had the third lowest with delta AIC = 2.79. I'll go with 7

##The delta AIC for those two is <2, so performing a likelihood ratio test

lrtest(model3, model7)

#They aren't significantly different, so use model 7 since it's simpler

##Assign the best fit glm's for each species

MIglmer <- feedRateDf%>%
  filter(site == "MI", round(numEaten) >0, pond != "Control", temp < 35)%>%
  glmer(cbind(round(numEaten), round(100-numEaten)) ~ temp + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

MOglmer <- feedRateDf%>%
  filter(site == "MO", round(numEaten) > 0, pond != "Control", temp < 35)%>%
  glmer(cbind(round(numEaten), round(100-numEaten)) ~ temp + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

TXglmer <- feedRateDf%>%
  filter(site == "TX", round(numEaten) > 0, pond != "Control", temp < 35)%>%
  glmer(cbind(round(numEaten), round(100-numEaten)) ~ temp + (1|pond) + (1|bath), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))


#Get model results
library(car)
MIglmer
Anova(MIglmer, 3)

MOglmer
Anova(MOglmer, 3)

TXglmer
Anova(TXglmer, 3)


#Get predicted values and add to data frame
##MI


###Make the predator masses equal average value to smooth the prediction curve
MIdata$predmass <- mean(MIdata$predmass) #Make the predator masses equal average value to smooth the prediction curve


###Calculate the predicted values with CI
library(bootpredictlme4)

MIfit <- predict(MIglmer, newdata = MIdata, re.form = NA, se.fit = T, nsim = 1000, type = "response")

###Create df with both the raw data and fit values
MIfit2 <- MIdata%>%
  mutate(fit = unname(data.frame(MIfit$fit)$MIfit.fit),
         lwr = unname(t(data.frame(MIfit$ci.fit))[,1]),
         upr = unname(t(data.frame(MIfit$ci.fit))[,2]))



##MO


###Make the predator masses equal average value to smooth the prediction curve
MOdata$predmass <- mean(MOdata$predmass) #Make the predator masses equal average value to smooth the prediction curve

###Calculate predicted values and CI

MOfit <- predict(MOglmer, newdata = MOdata, re.form = NA, se.fit = T, nsim = 1000, type = "response")

###Create df with both the raw data and fit values
MOfit2 <- MOdata%>%
  mutate(fit = unname(data.frame(MOfit$fit)$MOfit.fit),
         lwr = unname(t(data.frame(MOfit$ci.fit))[,1]),
         upr = unname(t(data.frame(MOfit$ci.fit))[,2]))

##TX

###Make the predator masses equal average value to smooth the prediction curve
TXdata$predmass <- mean(TXdata$predmass) 

###Calxulate predicted values with CI

TXfit <- predict(TXglmer, newdata = TXdata, re.form = NA, se.fit = T, nsim = 1000, type = "response")

###Create df with both the raw data and fit values
TXfit2 <- TXdata%>%
  mutate(fit = unname(data.frame(TXfit$fit)$TXfit.fit),
         lwr = unname(t(data.frame(TXfit$ci.fit))[,1]),
         upr = unname(t(data.frame(TXfit$ci.fit))[,2]))

###Combine each of these to a final data set
final.data <- rbind(MIfit2, MOfit2, TXfit2)




#############################
#Calculate Topt for MI
#########
#analysis of Topt differences

##For loop to bootstrap

#Create empty data frame that the Topt's will be inserted into
MIToptDf <- data.frame(Topt = NULL, bath = NULL, pond = NULL, predmass = NULL, sim = NULL)

#For loop to resample data, run model with resampled data, extract Topt, and add Topt to data frame
for (i in 1:500){
  
  #Resample data
  rand.df <- MIdata%>%
    group_by(pond)%>%
    slice_sample(prop = 1, replace = T)
  
  ##Run model with resampled data. Set it to continue the loop even if the model generates an error, which occasionally happens with some resampled data
  model <- rand.df%>%
    glmer(cbind(round(numEaten), round(100-numEaten)) ~ poly(temp, 2) + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = .,
          control = glmerControl(optimizer = "bobyqa"))
  
  #Create data set of dummy variables to get precise estimate of Topt, then extract predicted Topt.
 
  ndat <- data.frame(temp = seq(min(feedRateDf$temp), max(feedRateDf$temp), by = 0.01), 
                     predmass = mean(rand.df$predmass))
  ##Extract topt and add to df
  MIToptDf <- data.frame(Topt = ndat[which.max(predict(model, newdata = ndat, type = "response", re.form = ~0 )),1] , sim = i)%>%
    rbind(MIToptDf)
  
}

#Uncomment to save the randomized data set
#write.csv(MIToptDf, "Data/MITopt.csv", row.names = F)

#Uncomment if uploading a previously saved randomized df
#MIToptDf <- read.csv("Data/MITopt.csv")


#Create df with average Topt's and confidence intervals for plotting

MIavgToptDf <- MIToptDf%>%
  summarize(temp = mean(Topt),
            lwrCI = quantile(Topt, 0.025),
            uprCI = quantile(Topt, 0.975))%>%
  mutate(predmass = mean(MIdata$predmass))%>%
  mutate(numEaten = round(predict(MIglmer, newdata = ., type = "response", re.form = ~ 0 )))











###############################
#Graph it

##colorblind friendly palette

cbPalette <- c("#CC79A7", "#78C1EA", "#009E73", "#E69F00", "#D55E00", "#0072B2")

##Faceted by site
final.data%>%
  ggplot(aes(x = temp, y = fit*100, color = site, fill = site)) +
  facet_wrap(~site, labeller = labeller(site = c("MI" = "Michigan", "MO" = "Missouri", "TX" = "Texas"))) +
  geom_point(aes(y = numEaten)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(alpha = 0.25, aes(ymin = lwr*100, ymax = upr*100), linetype = 0) +
  theme_classic() +
  labs(x = "Temperature (°C)", y = "Feeding rate ± CI") +   
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size=20),
        strip.text = element_text(size=16),
        legend.position = 0) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) 
#geom_line(data = controlFit2, aes(x = temp, y = backMort), color = "Black") +
#geom_ribbon(data = controlFit2, aes(ymin = 100 - lwr*100, ymax = 100 - upr*100), 
#color = "black", alpha = 0.25, fill = "black", linetype = 0)




ggsave("Figures/Latresults1.pdf", width = 11.95, height = 5.37)

##All species in one panel WITHOUT error bars


final.data%>%
  ggplot(aes(x = temp, y = fit*100)) +
  #geom_point(aes(y = numEaten)) +
  geom_line(aes(color = site), size = 2) +
  theme_classic() +
  labs(x = "Temperature (°C)", y = "# Eaten") +   
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  scale_color_manual(values = cbPalette, name = "site") +
  scale_fill_manual(values = cbPalette) 



ggsave("Figures/Latresults2.pdf", width = 6.5, height = 6.0)


##All species in one panel WITH error bars

final.data%>%
  ggplot(aes(x = temp, y = fit*100)) +
  geom_line(aes(color = site), size = 2) +
  geom_ribbon(alpha = 0.25, aes(ymin = lwr*100, ymax = upr*100, fill = site)) +
  theme_classic() +
  labs(x = "Temperature (°C)", y = "# Eaten +/- SE") +   
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  scale_color_manual(values = cbPalette, name = "site") +
  scale_fill_manual(values = cbPalette) 

ggsave("Figures/Latresults3.pdf", width = 9.5, height = 5.9)




final.data%>%
  filter(site == "TX")%>%
  ggplot(aes(x = temp, y = fit*100, color = pond, fill = pond)) +
  # facet_wrap(~site) +
  geom_point(aes(y = numEaten)) +
  geom_line(size = 1.5) +
  geom_ribbon(alpha = 0.25, aes(ymin = lwr*100, ymax = upr*100), linetype = 0) +
  theme_classic() +
  labs(x = "Temperature (C)", y = "# Eaten +/- CI") +   
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size=20),
        strip.text = element_text(size=16, face = "italic"),
        legend.position = 0) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  ylim(0,50)






################################
#Full model with all Sites
library()
feedRateDf%>%
  glmer(cbind(round(numEaten), round(100-numEaten)) ~ temp*site + I(temp^2) + (1|pond) + (1|bath), data = .,
        family = "binomial", control = glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 1000)))

#Mixed effects model won't converge. Trying main effects model.
#Compare AICc's to select best fit model


feedRateDf%>%
  filter(pond != "Control")%>%
  glm(cbind(round(numEaten), round(100-numEaten)) ~ temp*site + I(temp^2) + scale(predmass), 
             family = "binomial", data = .)%>%
  plot()

#Top two models are full (AICc = 2577.966) and full without predmass (2686.456). ∆AICc is small. Try likelihood ratio test


lrtest(glm(cbind(round(numEaten), round(100-numEaten)) ~ temp*site + I(temp^2) , 
      family = "binomial", data = filter(feedRateDf, pond != "Control")),
      glm(cbind(round(numEaten), round(100-numEaten)) ~ temp*site + I(temp^2) + scale(predmass), 
          family = "binomial", data = filter(feedRateDf, pond != "Control")))

#Two models fit significantly differently. Move forward with the full model

fullGlm <- feedRateDf%>%
  glm(cbind(round(numEaten), round(100-numEaten)) ~ temp*site + I(temp^2) + scale(predmass), 
      family = "binomial", data = .)


#Get predicted values and add to observed value data set
feedRateDf2 <- feedRateDf%>%
  group_by(site)%>%
  mutate(predmass = mean(predmass)) %>%
  ungroup()%>%
  mutate(fit = predict(fullGlm, type = "response", newdata = .) *100,
         se = predict(fullGlm, type = "response", se.fit = T, newdata = .)$se.fit * 100) %>%
  group_by(pond)%>%
  mutate(predmass = mean(predmass))


feedRateDf2%>%
  ggplot(aes(x = temp, y = numEaten)) + 
    geom_point() +
    geom_line(aes(y = fit)) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se, alpha = 0.5)) +
    facet_wrap(vars(site)) +
    theme_classic()
  
