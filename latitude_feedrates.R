#Load data

library(tidyverse)

#Load data
#Feed trial data
df <- read.csv("Data/lat-feed-data.csv")%>%
  mutate(across(c(site, chamber, bath, pond, id), as.factor))

#Climate data
climate <- read.csv("Data/climates.csv")%>%
  mutate(across(site, as.factor),
         Tavg = (Tavg - 32)*5/9, #Convert temperature to celcius
         Tmin = (Tmin - 32)*5/9,
         Tmax = (Tmax - 32)*5/9)
         


#There's an outlier in the predator mass data. Replace with NA

df[73,6] <- NA

#Assign missing predator mass values with pond-level averages

df <- df%>%
  group_by(site, pond)%>%
  mutate_at(vars(predmass), ~replace_na(., mean(., na.rm = T)))

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
  ungroup()%>%
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
  labs(y = expression(paste(italic("Daphnia pulex"), " background mortality")),
       x = "Temperature (°C)") +
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

library(MuMIn)

#With pond, mixed effects
library(lme4)


fullGlmer <- feedRateDf%>%
  mutate(predmass = scale(predmass))%>%
  glmer(cbind(numEaten, 100-numEaten) ~ poly(temp, 2)*site + (1|predmass) + (1|pond), 
        family = "binomial", data = .)


library(car)
Anova(fullGlmer, type = 3)


#Extract model coefficients
coefDf <- summary(fullGlmer)$coefficients%>%
  data.frame()%>%
  rownames_to_column(var = "Term")

#Create table with site-specific coefficients on link and response scales. MO and TX coefficients need to be added to the MI coefficient, which is the baseline
data.frame(site = c(rep(c("MI", "MO", "TX"), 2)),
           term = c(rep("Temperature", 3), rep("Temperature2", 3)),
           coef = c(coefDf[2,2], #MI temp
                         coefDf[2,2] + coefDf[6,2], #MO temp
                         coefDf[2,2] + coefDf[8,2], #TX temp
                         coefDf[3,2], #MI temp2
                         coefDf[3,2] + coefDf[7,2], #MO temp2
                         coefDf[3,2] + coefDf[9,2]), #TX temp2
           se = c(coefDf[2,3], #MI temp
                       coefDf[2,3] + coefDf[6,3], #MO temp
                       coefDf[2,3] + coefDf[8,3], #TX temp
                       coefDf[3,3], #MI temp2
                       coefDf[3,3] + coefDf[7,3], #MO temp2
                       coefDf[3,3] + coefDf[9,3]))%>%#TX temp2
  view()
  

# #Create df of the site-specific average predator masses for generating smooth predicted curves
avgPredmassDf <- data.frame(site = c("MI", "MO", "TX"),
                            predmass = c(mean(filter(feedRateDf, site == "MI")$predmass),
                                         mean(filter(feedRateDf, site == "MO")$predmass),
                                         mean(filter(feedRateDf, site == "TX")$predmass)))



#Estimate fit values for plotting
library(AICcmodavg) #Package for generating standard errors from glmer objects


#glmer predictions

#Get reverse transformation function for generating predicted values on response scale
ilink <- fullGlmer@resp$family$linkinv

#Create df of fit values
fitFeedDf <- data.frame(temp = rep(seq(min(feedRateDf$temp), max(feedRateDf$temp), by = 0.1), times =9))%>%
  mutate(site = as.factor(rep(c("MI", "MO", "TX"), each = 1/3*length(temp))),
         pond = as.factor(rep(c("29", "36", "64", "Blue", "Grey", "Twin", "108", "17", "19"), each = 1/9*length(temp))))%>%
  left_join(avgPredmassDf, by = "site")%>%
  mutate(fit_link = predictSE(fullGlmer, type = "link", newdata = .,  re.form = NA)$fit)%>%
  mutate(se_link = predictSE(fullGlmer, se.fit = T, type = "link", newdata = ., re.form = NA)$se.fit)%>%
  mutate(fit_resp = (ilink(fit_link))*100,
         lwr = (ilink(fit_link - 1.96*se_link))*100,
         upr = (ilink(fit_link + 1.96*se_link))*100)
#############################
#Calculate Topt's with bootstrap analysis

#Create empty data frame that the Topt's will be inserted into
ToptDf <- data.frame(Topt = NULL, site = NULL, sim = NULL)

##Create data sets of dummy variables for each site to get more precise estimate of Topt in predict function

MIndat <- data.frame(temp = rep(seq(min(feedRateDf$temp), max(feedRateDf$temp), by = 0.1), times = 3),
                     site = "MI")%>%
  mutate(pond = rep(unique(filter(feedRateDf, site == "MI")$pond), times = 1/3*nrow(.)))%>%
  left_join(avgPredmassDf, by = "site")

MOndat <- data.frame(temp = rep(seq(min(feedRateDf$temp), max(feedRateDf$temp), by = 0.11), times = 3),
                     site = "MO")%>%
  mutate(pond = rep(unique(filter(feedRateDf, site == "MO")$pond), times = 1/3*nrow(.)))%>%
  left_join(avgPredmassDf, by = "site")

TXndat <- data.frame(temp = rep(seq(min(feedRateDf$temp), max(feedRateDf$temp), by = 0.1), times = 3),
                     site = "TX")%>%
  mutate(pond = rep(unique(filter(feedRateDf, site == "TX")$pond), times = 1/3*nrow(.)))%>%
  left_join(avgPredmassDf, by = "site")

#Create empty df to insert the randomized data used in each simulation
randomized <- data.frame(site = NULL, predmass = NULL, temp = NULL, numEaten = NULL, sim = NULL)

##For loop to resample data, run model with resampled data, extract Topt, and add Topt to data frame
for (i in 1:1000){

  #Resample data
  rand.df <- feedRateDf%>%
    dplyr::select(site, predmass, temp, numEaten, pond)%>%
    group_by(site)%>%
    slice_sample(prop = 1, replace = T)%>%
    ungroup()
  
  #Insert the resampled data into randomized df
  randomized <- rand.df%>%
    mutate(sim = i)%>%
    rbind(randomized)

  ##Run model with resampled data. Set it to continue the loop even if the model generates an error, which occasionally happens with some resampled data
  model <- rand.df%>%
    glmer(cbind(numEaten, 100 - numEaten) ~ poly(temp, 2)*site + (1|predmass) + (1|pond), 
        family = "binomial", data = .)

  ##Extract Michigan Topt and add to ToptDf
  ToptDf <- data.frame(Topt = MIndat[which.max(predict(model, newdata = MIndat, re.form = NA, 
                                                     type = "response", allow.new.levels = T)),1] , sim = i, site = (c("MI")))%>%
    rbind(ToptDf)

  ##Extract Missouri Topt and add to ToptDf
  ToptDf <- data.frame(Topt = MOndat[which.max(predict(model, newdata = MOndat, re.form = NA,
                                                     type = "response", allow.new.levels = T)),1] , sim = i, site = (c("MO")))%>%
    rbind(ToptDf)

  ##Extract Missouri Topt and add to ToptDf
  ToptDf <- data.frame(Topt = TXndat[which.max(predict(model, newdata = TXndat, re.form = NA,
                                                     type = "response", allow.new.levels = T)),1] , sim = i, site = (c("TX")))%>%
    rbind(ToptDf)
  
}

#Uncomment to save the randomized data set
#write.csv(ToptDf, "Data/Topt.csv", row.names = F)

#Uncomment if uploading a previously saved randomized df
#ToptDf <- read.csv("Data/Topt.csv")


#Create df with average Topt's and confidence intervals for plotting

avgToptDf <- ToptDf%>%
  group_by(site)%>%
  summarize(temp = mean(Topt),
            lwrCI = quantile(Topt, 0.025),
            uprCI = quantile(Topt, 0.975)) %>%
  mutate(predmass = avgPredmassDf$predmass)%>%  
  mutate(numEaten = predict(fullGlmer, newdata = ., type = "response", re.form = NA)*100)

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

#Climate plot

climate%>%
  ggplot(aes(x = month, y = Tavg)) +
  geom_smooth(color = "black", se = F) +
  geom_smooth(aes(y = Tmax), color = "black", se = F, linetype = "dashed") +
  geom_smooth(aes(y = Tmin), color = "black", se = F, linetype = "dashed") + 
  facet_grid(~site, labeller = labeller(site = c( "MI" = "Michigan", "MO" = "Missouri", "TX" = "Texas" )))+
  labs(x = "Month", y = "Temperature (°C)") +
  scale_x_continuous(breaks = seq(1,12),
                     labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 24, vjust = 0),
        axis.title.y = element_text(size = 24),
        axis.text.x = element_text(size=14, angle = 90, vjust = 0.3),
        axis.text.y = element_text(size = 20),
        strip.text = element_text(size=20),
        title = element_text(size = 20)) 
 
ggsave("Figures/climate.jpg", width = 11.95, height = 5.37)



 ##Feeding rate plot

#Assign ponds values of 1, 2, and 3 within sites for assigning shapes of points
feedRateDf <- feedRateDf%>%
  mutate(pondNum = as.factor(case_when(pond == "29" ~ 1, 
                                       pond == "36" ~ 2,
                                       pond == "64" ~ 3,
                                       pond == "Blue" ~ 1,
                                       pond == "Grey" ~ 2,
                                       pond == "Twin" ~ 3,
                                       pond == "108" ~ 1,
                                       pond == "17" ~ 2,
                                       pond == "19" ~ 3)))
##Df of historic mean temperatures of month of collection
siteTemp <- data.frame(site = c("MI", "MO", "TX"), temp = c(filter(climate, site == "MI", month == 7)$Tavg,
                                                            filter(climate, site == "MO", month == 7)$Tavg,
                                                            filter(climate, site == "TX", month == 8)$Tavg))

##Faceted by site
fitFeedDf%>%
  ggplot(aes(x = temp, y = fit_resp, color = site, fill = site)) +
  facet_wrap(~site, labeller = labeller(site = c("MI" = "Michigan", "MO" = "Missouri", "TX" = "Texas"))) +
  geom_point(data = feedRateDf, aes(y = numEaten, shape = pondNum)) +
  geom_line(linewidth = 1) +
  geom_ribbon(alpha = 0.2, aes(ymin = lwr, ymax = upr), linetype = 0) +
  geom_vline(data = siteTemp, aes(xintercept = temp), linetype = 3, linewidth = 0.7) +
  theme_classic() +
  labs(x = "Temperature (°C)", y = "Number of prey eaten") +   
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size=20),
        strip.text = element_text(size=20),
        legend.position = 0,
        title = element_text(size = 20)) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) +
  geom_point(data = avgToptDf, aes(x = temp, y = numEaten), color = "Black", size = 3) +
  geom_linerange(data = avgToptDf, aes(xmin = lwrCI, xmax = uprCI, y = numEaten), color = "Black")
#geom_line(data = controlFit2, aes(x = temp, y = backMort), color = "Black") + #These lines add the background mortality fit values
#geom_ribbon(data = controlFit2, aes(ymin = 100 - lwr*100, ymax = 100 - upr*100), 
#color = "black", alpha = 0.25, fill = "black", linetype = 0)


ggsave("Figures/latfig.jpeg", width = 11.95, height = 5.37)


