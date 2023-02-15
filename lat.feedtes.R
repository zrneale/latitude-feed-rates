
library(tidyverse)




#Load data

##Survivorship data
survdf <- read.csv("Data/Lat.data.csv")

##Predator Mass data
massdf <- read.csv("Data/Lat.pred.mass.csv")


##Combine the two
df <- massdf%>%
  dplyr::select(Site, Chamber, Bath, Pond, ID, Predmass)%>%
  right_join(survdf, by = c("Site", "Chamber", "Bath", "Pond", "ID"))%>%
  mutate_at(vars(Site, Chamber, Bath, Pond, ID), as.factor)



##A bit of data tidying. I'll set any survivor values above 100 to equal 100 to avoid negative mortality. Also, I need to add a column for number dead for the binomial glm, create a column with temperature in C, and drop rows missing survivor values.


df <- df %>%
  mutate(Surv = replace(Surv, Surv > 100, 100), Dead = 100 - Surv, Temp2C = (Temp2 - 32)*5/9)%>%
  drop_na(Surv)

##Visualize control data
df%>%
  filter(Pond == "Control")%>%
  ggplot(aes(x = Temp2C, y = Surv)) +
  geom_point()

##There's an outlier at the hottest temperature. This seems to track with what I've seen before that daphnia survivorship crashes as soon as the temp gets above 35C. Maybe it would be best to exclude data points from that temp? It would mean losing 3 points from the Michigan data. I'll have to talk to Volker about that


##Visualize predator data
df%>%
  filter(Pond != "Control")%>%
  ggplot(aes(x = Temp2C, y = Surv, color = Site)) + 
  geom_point()


#Calculate background mortality with binomial glmer


library(lme4)

##Create a subset of the data with just control treatments
Control.data <- df%>%
  filter(Pond == "Control")

##Run to check resid plot
backmort.glmer1 <- Control.data%>%
  glmer(data = ., cbind(Surv, Dead) ~ Temp2C  + (1|Site) + (1|Bath), family = "binomial")

plot(backmort.glmer1)

##Residual plot doesn't look good. There's a point way off to the left that's probably that point above 35C. #The rest of it also looks off, but it's hard to diagnose since it's so zoomed out. Running it again without #that point.

backmort.glmer2 <- Control.data%>%
  filter(Temp2C < 35)%>%
  glmer(data = ., cbind(Surv, Dead) ~ Temp2C  + (1|Site) + (1|Bath), family = "binomial")

plot(backmort.glmer2)



##That doesn't look as bad. I'll move forward without that point, but it means I'll need to remove the predator data points in that bath too.


##Trying another method for getting predicted values.
library(bootpredictlme4)

Control.predicts <- Control.data%>%
  filter(Temp2C < 35)%>%
  predict(backmort.glmer2, newdata = ., re.form = NA, se.fit = T, nsim = 100, type = "response")

Control.predicts2 <- Control.data%>%
  filter(Temp2C<35)%>%
  mutate(fit = unname(data.frame(Control.predicts$fit)$Control.predicts.fit),
         lwr = unname(t(data.frame(Control.predicts$ci.fit))[,1]),
         upr = unname(t(data.frame(Control.predicts$ci.fit))[,2]),
         Backmort = 100 - fit*100)

##Plot with data
Control.predicts2%>%
  ggplot(aes(x = Temp2C, y = Dead)) +
  geom_line(aes(y = 100 - fit*100 )) +
  geom_point() +
  geom_ribbon(aes(ymax = 100 - upr*100, ymin = 100 - lwr*100), alpha = 0.25, linetype = 0) +
  theme_classic() +
  geom_point(data = filter(df, Pond == "Control", Temp2C >35), shape = 1) +
  labs(y = "Background mortality ± CI",
       x = "Temperature (C)") +
  theme(axis.title = element_text(size = 18),
        axis.text = element_text(size=16))
       
#ggsave("Figures/backmort.jpg")

##Create a final data set 
feed.data <- df%>%
  filter(Temp2C < 35)%>%
  left_join(dplyr::select(Control.predicts2, Site, Chamber, Bath, Backmort), by = c("Site", "Chamber", "Bath"))%>%
  mutate(Numeaten = Dead - Backmort)%>%
  mutate(Numeaten = replace(Numeaten, Numeaten < 0, 0))





#Now to run the glm's on the predator data. I'll try separate models for the 3 sites. First I'll create a function to run all possible models to compare AIC's within site for model selection

library(MuMIn)

selectmodel<- function(State){
  
  ##Subset the data
  df <- feed.data%>%
    filter(Site == State, Numeaten > 0, Pond != "Control")
  
  
  
  ##Create each of the models
  model1 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp2C, 2) + (1|Pond) + (1|Bath) + scale(Predmass), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model2 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp2C, 2) + (1|Bath) + scale(Predmass), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model3 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp2C + (1|Pond) + (1|Bath) + scale(Predmass), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model4 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp2C + (1|Bath) + scale(Predmass), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model5 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp2C, 2) + (1|Pond) + (1|Bath), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model6 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp2C, 2) + (1|Bath), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model7 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp2C + (1|Pond) + (1|Bath), family = "binomial", data = df,
                  control = glmerControl(optimizer = "bobyqa"))
  model8 <- glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp2C + (1|Bath), family = "binomial", data = df,
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
  filter(Site == "TX", Numeaten > 0, Pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp2C, 2) + (1|Pond) + (1|Bath) + scale(Predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))%>%str()

#Model 2
feed.data%>%
  filter(Site == "TX", Numeaten > 0, Pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp2C, 2) + (1|Bath) + scale(Predmass), family = "binomial", data =.,
        control = glmerControl(optimizer = "bobyqa"))

#Model 3
feed.data%>%
  filter(Site == "TX", Numeaten > 0, Pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp2C + (1|Pond) + (1|Bath) + scale(Predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Model 4
feed.data%>%
  filter(Site == "TX", Numeaten > 0, Pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp2C + (1|Bath) + scale(Predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Model 5
feed.data%>%
  filter(Site == "TX", Numeaten > 0, Pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp2C, 2) + (1|Pond) + (1|Bath), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Model 6
feed.data%>%
  filter(Site == "TX", Numeaten > 0, Pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp2C, 2) + (1|Bath), family = "binomial", data =.,
        control = glmerControl(optimizer = "bobyqa"))

#Model 7
feed.data%>%
  filter(Site == "TX", Numeaten > 0, Pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp2C + (1|Pond) + (1|Bath), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Model 8
feed.data%>%
  filter(Site == "TX", Numeaten > 0, Pond != "Control")%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp2C + (1|Bath), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

#Models 3 and 7 are the ones that are converging. 7 had the lowest AIC of the eight models, and 3 had the third lowest with delta AIC = 2.76. I'll go with 7

##The delta AIC for those two is <2, so performing a likelihood ratio test



##Running the best fit glm's for each species

MIglmer <- feed.data%>%
  filter(Site == "MI", round(Numeaten) >0, Pond != "Control", Temp2C < 35)%>%
  drop_na(Predmass)%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ poly(Temp2C, 2) + (1|Pond) + (1|Bath) + scale(Predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

MOglmer <- feed.data%>%
  filter(Site == "MO", round(Numeaten) > 0, Pond != "Control", Temp2C < 35)%>%
  drop_na(Predmass)%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp2C + (1|Pond) + (1|Bath) + scale(Predmass), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))

TXglmer <- feed.data%>%
  filter(Site == "TX", round(Numeaten) > 0, Pond != "Control", Temp2C < 35)%>%
  drop_na(Predmass)%>%
  glmer(cbind(round(Numeaten), round(100-Numeaten)) ~ Temp2C + (1|Pond) + (1|Bath), family = "binomial", data = .,
        control = glmerControl(optimizer = "bobyqa"))



##MI

###Create the subset data
MIdata <- feed.data%>%
  filter(Site == "MI", round(Numeaten) > 0, Pond != "Control", Temp2C < 35)%>%
  drop_na(Predmass) #I'll probably need to find a different solution to NA's

###Make the predator masses equal average value to smooth the prediction curve
MIdata$Predmass <- mean(MIdata$Predmass) #Make the predator masses equal average value to smooth the prediction curve


###Calculate the predicted values with CI
library(bootpredictlme4)

MI.predict <- predict(MIglmer, newdata = MIdata, re.form = NA, se.fit = T, nsim = 1000, type = "response")

###Create df with both the raw data and fit values
MI.predict2 <- MIdata%>%
  mutate(fit = unname(data.frame(MI.predict$fit)$MI.predict.fit),
         lwr = unname(t(data.frame(MI.predict$ci.fit))[,1]),
         upr = unname(t(data.frame(MI.predict$ci.fit))[,2]))



##MO
###Create the subset data
MOdata <- feed.data%>%
  filter(Site == "MO", round(Numeaten) > 0, Pond != "Control", Temp2C < 35)%>%
  drop_na(Predmass) #I'll need to find a better thing to do with NA'a

###Make the predator masses equal average value to smooth the prediction curve
MOdata$Predmass <- mean(MOdata$Predmass) #Make the predator masses equal average value to smooth the prediction curve

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
  filter(Site == "TX", round(Numeaten) > 0, Pond != "Control", Temp2C < 35)

###Make the predator masses equal average value to smooth the prediction curve
TXdata$Predmass <- mean(TXdata$Predmass) 

###Calxulate predicted values with CI

TX.predict <- predict(TXglmer, newdata = TXdata, re.form = NA, se.fit = T, nsim = 1000, type = "response")

###Create df with both the raw data and fit values
TX.predict2 <- TXdata%>%
  mutate(fit = unname(data.frame(TX.predict$fit)$TX.predict.fit),
         lwr = unname(t(data.frame(TX.predict$ci.fit))[,1]),
         upr = unname(t(data.frame(TX.predict$ci.fit))[,2]))

###Combine each of these to a final data set
final.data <- rbind(MI.predict2, MO.predict2, TX.predict2)


#Graph it

##colorblind friendly palette

cbPalette <- c("#CC79A7", "#78C1EA", "#009E73", "#E69F00", "#D55E00", "#0072B2")

##Faceted by Site
final.data%>%
  ggplot(aes(x = Temp2C, y = fit*100, color = Site, fill = Site)) +
  facet_wrap(~Site, labeller = labeller(Site = c("MI" = "Michigan", "MO" = "Missouri", "TX" = "Texas"))) +
  geom_point(aes(y = Numeaten)) +
  geom_line(linewidth = 1.5) +
  geom_ribbon(alpha = 0.25, aes(ymin = lwr*100, ymax = upr*100), linetype = 0) +
  theme_classic() +
  labs(x = "Temperature (°C)", y = "Number Eaten ± CI") +   
  theme(axis.title = element_text(size = 24),
        axis.text = element_text(size=20),
        strip.text = element_text(size=16),
        legend.position = 0) +
  scale_color_manual(values = cbPalette) +
  scale_fill_manual(values = cbPalette) 
#geom_line(data = Control.predicts2, aes(x = Temp2C, y = Backmort), color = "Black") +
#geom_ribbon(data = Control.predicts2, aes(ymin = 100 - lwr*100, ymax = 100 - upr*100), 
#color = "black", alpha = 0.25, fill = "black", linetype = 0)




ggsave("Figures/Latresults1.pdf", width = 11.95, height = 5.37)

##All species in one panel WITHOUT error bars


final.data%>%
  ggplot(aes(x = Temp2C, y = fit*100)) +
  #geom_point(aes(y = Numeaten)) +
  geom_line(aes(color = Site), size = 2) +
  theme_classic() +
  labs(x = "Temperature (°C)", y = "# Eaten") +   
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  scale_color_manual(values = cbPalette, name = "Site") +
  scale_fill_manual(values = cbPalette) 



ggsave("Figures/Latresults2.pdf", width = 6.5, height = 6.0)


##All species in one panel WITH error bars

final.data%>%
  ggplot(aes(x = Temp2C, y = fit*100)) +
  geom_line(aes(color = Site), size = 2) +
  geom_ribbon(alpha = 0.25, aes(ymin = lwr*100, ymax = upr*100, fill = Site)) +
  theme_classic() +
  labs(x = "Temperature (°C)", y = "# Eaten +/- SE") +   
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size=18),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16)) +
  scale_color_manual(values = cbPalette, name = "Site") +
  scale_fill_manual(values = cbPalette) 

ggsave("Figures/Latresults3.pdf", width = 9.5, height = 5.9)




final.data%>%
  filter(Site == "TX")%>%
  ggplot(aes(x = Temp2C, y = fit*100, color = Pond, fill = Pond)) +
  # facet_wrap(~Site) +
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

