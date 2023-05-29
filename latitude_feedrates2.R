#Load data

library(tidyverse)

df <- read.csv("Data/lat-feed-data.csv")%>%
  mutate(across(c(site, chamber, bath, pond, id), as.factor))

#Calculate background mortality with binomial glmer

##Create a subset of the data with just control treatments. Remove the bath over 35C because that appears to be outlier
controlDf <- df%>%
  filter(pond == "Control")%>%
  filter(temp < 35)

#Run glmer. Remove the water bath above 35C because it had very high background mortality of prey    
library(lme4)
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
         backMort = round(100 - fit*100))

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
  filter(numEaten > 0)%>%
  filter(pond != "Control")



################################
#Full model with all Sites

#Try mixed effects first
feedRateDf%>%
  glmer(cbind(numEaten, 100-numEaten) ~ temp*site + I(temp^2) + (1|pond) + (1|bath), data = .,
        family = "binomial", control = glmerControl(optimizer = "bobyqa",optCtrl = list(maxfun = 1000)))

#Mixed effects model won't converge. Trying main effects model.

#Compare AICc's to select best fit model

library(MuMIn)
feedRateDf%>%
  glm(cbind(numEaten, 100-numEaten) ~ temp*site + I(temp^2)*site + scale(predmass), 
      family = "binomial", data = .)%>%
  AICc()

#Top two models are full (2119.671) and full without predmass (2124.186). ∆AICc = 4.515
#Run likelihood ratio test

lrtest(glm(cbind(numEaten, 100-numEaten) ~ temp*site + I(temp^2)*site , 
           family = "binomial", data = feedRateDf),
       glm(cbind(numEaten, 100-numEaten) ~ temp*site + I(temp^2)*site + scale(predmass), 
           family = "binomial", data = feedRateDf))


#Two models fit significantly differently. Move forward with the full model

fullGlm <- feedRateDf%>%
  glm(cbind(numEaten, 100-numEaten) ~ temp*site + site*I(temp^2) + scale(predmass), 
      family = "binomial", data = .)

cooks.distance(fullGlm)

as.numeric(names(cooks.distance(fullGlm))[(cooks.distance(fullGlm) > (4/nrow(feedRateDf)))])


# #Create df of the site-specific average predator masses for generating predicted values
avgPredmassDf <- data.frame(site = c("MI", "MO", "TX"),
                            predmass = c(mean(filter(feedRateDf, site == "MI")$predmass),
                                         mean(filter(feedRateDf, site == "MO")$predmass),
                                         mean(filter(feedRateDf, site == "TX")$predmass)))
#Create df of predicted values
fitFeedDf <- data.frame(temp = rep(seq(min(feedRateDf$temp), max(feedRateDf$temp), by = 0.01), times = 3),
           site = as.factor(rep(c("MI", "MO", "TX"), each = 2729)))%>%
  left_join(avgPredmassDf, by = "site")%>%
  mutate(fit = predict(fullGlm, type = "response", newdata = .) *100,
         se = predict(fullGlm,  se.fit = T, type = "response", newdata = .)$se.fit * 100)



#############################
#Calculate Topt's with bootstrap analysis

#Create empty data frame that the Topt's will be inserted into
ToptDf <- data.frame(Topt = NULL, site = NULL, sim = NULL)

##Create data set of dummy variables to get more precise estimate of Topt in predict function

ndat <- data.frame(temp = rep(seq(min(feedRateDf$temp), max(feedRateDf$temp), by = 0.01), each = 3))%>%
  mutate(site = as.factor(rep(c("MI", "MO", "TX"), each = 1, length.out = length(temp))))%>%
  left_join(avgPredmassDf, by = "site")


##For loop to resample data, run model with resampled data, extract Topt, and add Topt to data frame
for (i in 1:10){

  #Resample data
  rand.df <- feedRateDf%>%
    group_by(site)%>%
    slice_sample(prop = 1, replace = T)

  ##Run model with resampled data. Set it to continue the loop even if the model generates an error, which occasionally happens with some resampled data
  model <- rand.df%>%
    glm(cbind(numEaten, 100 - numEaten) ~ temp*site + I(temp^2)*site + scale(predmass), 
        family = "binomial", data = .)

  ##Extract Michigan Topt and add to ToptDf
  ToptDf <- data.frame(Topt = ndat[which.max(predict(model, newdata = filter(ndat, site =="MI"), 
                                                     type = "response")),1] , sim = i, site = (c("MI")))%>%
    rbind(ToptDf)

  ##Extract Missouri Topt and add to ToptDf
  ToptDf <- data.frame(Topt = ndat[which.max(predict(model, newdata = filter(ndat, site == "MO"), 
                                                     type = "response")),1] , sim = i, site = (c("MO")))%>%
    rbind(ToptDf)

  ##Extract Missouri Topt and add to ToptDf
  ToptDf <- data.frame(Topt = ndat[which.max(predict(model, newdata = filter(ndat, site == "TX"), 
                                                     type = "response")),1] , sim = i, site = (c("TX")))%>%
    rbind(ToptDf)
}

#Uncomment to save the randomized data set
#write.csv(ToptDf, "Data/Topt.csv", row.names = F)

#Uncomment if uploading a previously saved randomized df
ToptDf <- read.csv("Data/Topt.csv")


#Create df with average Topt's and confidence intervals for plotting

avgToptDf <- ToptDf%>%
  group_by(site)%>%
  summarize(temp = mean(Topt),
            lwrCI = quantile(Topt, 0.025),
            uprCI = quantile(Topt, 0.975)) %>%
  mutate(predmass = avgPredmassDf$predmass)%>%
  mutate(numEaten = predict(fullGlm, newdata = ., type = "response")*100)

#Calculate differences between site Topts in each simulated analysis
ToptDiffDf <- ToptDf%>%
  pivot_wider(names_from = site, values_from = Topt)%>%
  rowwise()%>%
  summarize("TX-MO" = TX - MO,
            "TX-MI" = TX - MI,
            "MO-MI" = MO - MI)%>%
  summarize(across(.cols = everything(), c(mean = mean,
                                           uprCI = ~ quantile(., probs = c(0.975)),
                                           lwrCI = ~ quantile(., probs = c(0.025)))))%>%
  pivot_longer(everything())%>%
  separate(name, c("sites", "statistic"), "_")%>%
  pivot_wider(id_cols = sites, names_from = statistic)



###############################
#Graphs

##colorblind friendly palette

cbPalette <- c("#CC79A7", "#78C1EA", "#009E73", "#E69F00", "#D55E00", "#0072B2")

##Faceted by site
fitFeedDf%>%
  ggplot(aes(x = temp, y = fit, color = site, fill = site)) +
  facet_wrap(~site, labeller = labeller(site = c("MI" = "Michigan", "MO" = "Missouri", "TX" = "Texas"))) +
  geom_point(data = feedRateDf, aes(y = numEaten)) +
  geom_line(linewidth = 0.3) +
  geom_ribbon(alpha = 0.3, aes(ymin = fit - se, ymax = fit + se), linetype = 0) +
  theme_classic() +
  labs(x = "Temperature (°C)", y = "Number of Prey Eaten ± CI") +   
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size=20),
        strip.text = element_text(size=16),
        legend.position = 0) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  geom_point(data = avgToptDf, aes(x = temp, y = numEaten), color = "Black", size = 3) +
  geom_linerange(data = avgToptDf, aes(xmin = lwrCI, xmax = uprCI, y = numEaten), color = "Black")
#geom_line(data = controlFit2, aes(x = temp, y = backMort), color = "Black") +
#geom_ribbon(data = controlFit2, aes(ymin = 100 - lwr*100, ymax = 100 - upr*100), 
#color = "black", alpha = 0.25, fill = "black", linetype = 0)


#ggsave("Figures/Latresults1.pdf", width = 11.95, height = 5.37)


