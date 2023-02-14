

library(tidyverse)

#Load data
Tempdf <- read.csv("Data/Monthly.temps.csv", header = T)%>%
  mutate_at(vars(Month, State), as.factor)%>%
  mutate(TempC = (TempF - 32)*5/9) #Convert temperature to Celsius

#Calculate summary stats

Avgtemps <- Tempdf%>%
  group_by(State)%>%
  summarize(meantemp = mean(TempC), stdv = sd(TempC))


#Plot them
##Create the text labels and fonts to standardize across plots



textsize <- theme(axis.title = element_blank(),
                  axis.text.y = element_text(size = 14),
                  axis.text.x = element_blank())


#All three states together
Avgtemps%>%
  ggplot(aes(x = State, y = meantemp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = meantemp - stdv, ymax = meantemp + stdv)) +
  theme_classic() +
  labels +
  textsize

#Michigan
Avgtemps%>%
  filter(State == "MI")%>%
  ggplot(aes(x = State, y = meantemp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = meantemp - stdv, ymax = meantemp + stdv)) +
  theme_classic() +
  labels +
  textsize +
  ylim(-1,27)

ggsave("Figures/Michigan.temp.pdf", height = 1.61, width = 1.12)

#Missouri
Avgtemps%>%
  filter(State == "MO")%>%
  ggplot(aes(x = State, y = meantemp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = meantemp - stdv, ymax = meantemp + stdv)) +
  theme_classic() +
  labels +
  textsize +
  ylim(-1,27)

ggsave("Figures/Missouri.temp.pdf", height = 1.61, width = 1.12)

#Texas
#Missouri
Avgtemps%>%
  filter(State == "TX")%>%
  ggplot(aes(x = State, y = meantemp)) +
  geom_point(size = 4) +
  geom_linerange(aes(ymin = meantemp - stdv, ymax = meantemp + stdv)) +
  theme_classic() +
  labels +
  textsize +
  ylim(-1,27)

ggsave("Figures/Texas.temp.pdf", height = 1.61, width = 1.12)
