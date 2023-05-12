
library(tidyverse)




#Load data

df <- read.csv("Data/lat-feed-data.csv")
  

#Calculate background mortality with binomial glmer


library(lme4)

##Create a subset of the data with just control treatments. Remove the bath over 35C because that appears to be outlier
Control.data <- df%>%
  filter(pond == "Control")%>%
  filter(temp < 35)

#Run glmer. Remove the water bath above 35C because it had very high background mortality of prey     
backmort.glmer <- Control.data%>%
  glmer(data = ., cbind(surv, dead) ~ temp  + (1|site) + (1|bath), family = "binomial")


##Obtain predicted values of background mortality
library(bootpredictlme4)

Control.predicts <- Control.data%>%
  predict(backmort.glmer, newdata = ., re.form = NA, se.fit = T, nsim = 100, type = "response")

Control.predicts2 <- Control.data%>%
  mutate(fit = unname(data.frame(Control.predicts$fit)$Control.predicts.fit),
         lwr = unname(t(data.frame(Control.predicts$ci.fit))[,1]),
         upr = unname(t(data.frame(Control.predicts$ci.fit))[,2]),
         Backmort = 100 - fit*100)

##Plot with data
Control.predicts2%>%
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
       
ggsave("Figures/backmort.jpg")

##Create a final data set 
feed.data <- df%>%
  filter(temp < 35)%>%
  left_join(dplyr::select(Control.predicts2, site, chamber, bath, Backmort), by = c("site", "chamber", "bath"))%>%
  mutate(Numeaten = dead - Backmort)%>%
  mutate(Numeaten = replace(Numeaten, Numeaten < 0, 0))





#Now to run the glm's on the predator data. I'll try separate models for the 3 sites. First I'll create a function to run all possible models to compare AIC's within site for model selection

library(MuMIn)

selectmodel<- function(State){
  
  ##Subset the data
  df <- feed.data%>%
    filter(site == State, Numeaten > 0, pond != "Control")
  
  
  
  ##Create each of the models
  model1 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model2 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|bath) + scale(predmass), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model3 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model4 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|bath) + scale(predmass), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model5 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|pond) + (1|bath), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model6 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|bath), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model7 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|pond) + (1|bath), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model8 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|bath), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  
  
  data.frame(model = c("model1", "model2", "model3", "model4", "model5", "model6", "model7", "model8"),
             AICc = c(AICc(model1), AICc(model2), AICc(model3), AICc(model4),
                     AICc(model5), AICc(model6), AICc(model7), AICc(model8)))%>%
    arrange(AIC)%>%
    mutate(delta = AIC - .[1,2])%>%
    return()
  
  
}

##Run the function to compare models
selectmodel("MI")
selectmodel("MO")
selectmodel("TX")

##Six of the texas models are giving singularity warnings. Running each one here to figure out which ones

#Model 1
feed.data%>%
  filter(site == "TX", Numeaten > 0, pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))%>%str()

#Model 2
feed.data%>%
  filter(site == "TX", Numeaten > 0, pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|bath) + scale(predmass), family = "binomial", data =.,
        control = glmerControl(optimizer = "bobyqa"))

#Model 3
feed.data%>%
  filter(site == "TX", Numeaten > 0, pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Model 4
feed.data%>%
  filter(site == "TX", Numeaten > 0, pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|bath) + scale(predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Model 5
feed.data%>%
  filter(site == "TX", Numeaten > 0, pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|pond) + (1|bath), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Model 6
feed.data%>%
  filter(site == "TX", Numeaten > 0, pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|bath), family = "binomial", data =.,
        control = glmerControl(optimizer = "bobyqa"))

#Model 7
feed.data%>%
  filter(site == "TX", Numeaten > 0, pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|pond) + (1|bath), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Model 8
feed.data%>%
  filter(site == "TX", Numeaten > 0, pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|bath), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Models 3 and 7 are the ones that are converging. 7 had the lowest AIC of the eight models, and 3 had the third lowest with delta AIC = 2.76. I'll go with 7

##The delta AIC for those two is <2, so performing a likelihood ratio test



##Running the best fit glm's for each species

MIglmer <- feed.data%>%
  filter(site == "MI", round(Numeaten) >0, pond != "Control", temp < 35)%>%
  drop_na(predmass)%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(temp, 2) + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

MOglmer <- feed.data%>%
  filter(site == "MO", round(Numeaten) > 0, pond != "Control", temp < 35)%>%
  drop_na(predmass)%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|pond) + (1|bath) + scale(predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

TXglmer <- feed.data%>%
  filter(site == "TX", round(Numeaten) > 0, pond != "Control", temp < 35)%>%
  drop_na(predmass)%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ temp + (1|pond) + (1|bath), family = "binomial", data = .,
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

###Create the subset data
Midata <- feed.data%>%
  filter(site == "MI", round(Numeaten) > 0, pond != "Control", temp < 35)%>%
  drop_na(predmass) #I'll probably need to find a different solution to NA's

###Make the predator masses equal average value to smooth the prediction curve
Midata$predmass <- mean(Midata$predmass) #Make the predator masses equal average value to smooth the prediction curve


###Calculate the predicted values with CI
library(bootpredictlme4)

MI.predict <- predict(MIglmer, newdata = Midata, re.form = NA, se.fit = T, nsim = 1000, type = "response")

###Create df with both the raw data and fit values
MI.predict2 <- Midata%>%
  mutate(fit = unname(data.frame(MI.predict$fit)$MI.predict.fit),
         lwr = unname(t(data.frame(MI.predict$ci.fit))[,1]),
         upr = unname(t(data.frame(MI.predict$ci.fit))[,2]))



##MO
###Create the subset data
MOdata <- feed.data%>%
  filter(site == "MO", round(Numeaten) > 0, pond != "Control", temp < 35)%>%
  drop_na(predmass) #I'll need to find a better thing to do with NA'a

###Make the predator masses equal average value to smooth the prediction curve
MOdata$predmass <- mean(MOdata$predmass) #Make the predator masses equal average value to smooth the prediction curve

###Calculate predicted values and CI

MO.predict <- predict(MOglmer, newdata = MOdata, re.form = NA, se.fit = T, nsim = 1000, type = "response")

###Create df with both the raw data and fit values
MO.predict2 <- MOdata%>%
  mutate(fit = unname(data.frame(MO.predict$fit)$MO.predict.fit),
         lwr = unname(t(data.frame(MO.predict$ci.fit))[,1]),
         upr = unname(t(data.frame(MO.predict$ci.fit))[,2]))

##TX
###Create the subset data
TXdata <- feed.data%>%
  filter(site == "TX", round(Numeaten) > 0, pond != "Control", temp < 35)

###Make the predator masses equal average value to smooth the prediction curve
TXdata$predmass <- mean(TXdata$predmass) 

###Calxulate predicted values with CI

TX.predict <- predict(TXglmer, newdata = TXdata, re.form = NA, se.fit = T, nsim = 1000, type = "response")

###Create df with both the raw data and fit values
TX.predict2 <- TXdata%>%
  mutate(fit = unname(data.frame(TX.predict$fit)$TX.predict.fit),
         lwr = unname(t(data.frame(TX.predict$ci.fit))[,1]),
         upr = unname(t(data.frame(TX.predict$ci.fit))[,2]))

###Combine each of these to a final data set
final.data <- rbind(MI.predict2, MO.predict2, TX.predict2)




#############################
#Calculate Topt for MI
#########
#analysis of Topt differences

##For loop to bootstrap

#Create empty data frame that the Topt's will be inserted into
MIToptDf <- data.frame(Topt = NULL, bath = NULL, pond = NULL, sim = NULL)

#For loop to resample data, run model with resampled data, extract Topt, and add Topt to data frame
for (i in 1:1000){
  
  #Resample data
  rand.df <- pulexDf%>%
    filter(pulexsurv > 10)%>%
    group_by(pred, comp)%>%
    slice_sample(prop = 1, replace = T)
  
  ##Run model with resampled data. Set it to continue the loop even if the model generates an error, which occasionally happens with some resampled data
  try(model <- rand.df%>%
        glm.nb(pulexsurv ~ temp*pred*comp + I(temp^2), data = ., control = glm.control(maxit = 10000)), silent = TRUE)
  
  #Create data set of dummy variables to get precise estimate of Topt, then extract predited Topt. This will be done for each comp x pred combination
  
  ##Create no pred, no comp dummy data
  ndat <- data.frame(temp = seq(10,31, by = 0.01), pred = "N", comp = "N")
  ##Extract no pred, no comp topt and add to df
  pulexToptDf <- data.frame(Topt = ndat[which.max(predict(model, newdata = ndat, type = "response")),1],
                            pred = "N", comp = "N", sim = i)%>%
    rbind(pulexToptDf)
  
  ##Create no pred, yes comp dummy data
  ndat <- data.frame(temp = seq(10,31, by = 0.01), pred = "N", comp = "Y")
  ##Extract no pred, no comp topt and add to df
  pulexToptDf <- data.frame(Topt = ndat[which.max(predict(model, newdata = ndat, type = "response")),1],
                            pred = "N", comp = "Y", sim = i)%>%
    rbind(pulexToptDf)
  
  ##Create yes pred, no comp dummy data
  ndat <- data.frame(temp = seq(10,31, by = 0.01), pred = "Y", comp = "N")
  ##Extract yes pred, no comp topt and add to df
  pulexToptDf <- data.frame(Topt = ndat[which.max(predict(model, newdata = ndat, type = "response")),1],
                            pred = "Y", comp = "N", sim = i)%>%
    rbind(pulexToptDf)
  
  ##Create yes pred, no comp dummy data
  ndat <- data.frame(temp = seq(10,31, by = 0.01), pred = "Y", comp = "Y")
  ##Extract no pred, no comp topt and add to df
  pulexToptDf <- data.frame(Topt = ndat[which.max(predict(model, newdata = ndat, type = "response")),1],
                            pred = "Y", comp = "Y", sim = i)%>%
    rbind(pulexToptDf)
}

#Uncomment to save the randomized data set
#write.csv(pulexToptDf, "Data/pulexTopt.csv", row.names = F)

#Uncomment if uploading a previously saved randomized df
#pulexToptDf <- read.csv("Data/pulexTopt.csv")


#Create df with average Topt's and confidence intervals for plotting

pulexavgToptDf <- pulexToptDf%>%
  group_by(pred, comp)%>%
  summarize(temp = mean(Topt),
            lwrCI = quantile(Topt, 0.025),
            uprCI = quantile(Topt, 0.975))%>%
  ungroup()%>%
  mutate(pulexsurv = round(predict(pulexNegBin, newdata = ., type = "response")))










###############################
#Graph it

##colorblind friendly palette

cbPalette <- c("#CC79A7", "#78C1EA", "#009E73", "#E69F00", "#D55E00", "#0072B2")

##Faceted by site
final.data%>%
  ggplot(aes(x = temp, y = fit*100, color = site, fill = site)) +
  facet_wrap(~site, labeller = labeller(site = c("MI" = "Michigan", "MO" = "Missouri", "TX" = "Texas"))) +
  geom_point(aes(y = Numeaten)) +
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
#geom_line(data = Control.predicts2, aes(x = temp, y = Backmort), color = "Black") +
#geom_ribbon(data = Control.predicts2, aes(ymin = 100 - lwr*100, ymax = 100 - upr*100), 
#color = "black", alpha = 0.25, fill = "black", linetype = 0)




ggsave("Figures/Latresults1.pdf", width = 11.95, height = 5.37)

##All species in one panel WITHOUT error bars


final.data%>%
  ggplot(aes(x = temp, y = fit*100)) +
  #geom_point(aes(y = Numeaten)) +
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
  geom_point(aes(y = Numeaten)) +
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

