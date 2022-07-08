rm(list=ls())


##--------------
## LOAD PACKAGES
##--------------
library(colormap)
library(GGally)
library(ggplot2)
library(gtsummary)
library(reshape2)
library(splitstackshape)
library(survBootOutliers)
library(survival)
library(survminer)
library(tidyverse)


##-------------
## LOAD DATASET
##-------------
data <- read.csv("data/insect_survival_experiment_data.txt", sep = "")

data_death <- data %>% select(Island, Year, Condition, Species, Order, Status, Days, D_SL)
data_death <- data_death %>% filter(D_SL > 0)  
data_death <- as.data.frame(expandRows(data_death, "D_SL"))
rownames(data_death) <- NULL
data_death$event <- 1

data_alive <- data %>% select(Island, Year, Condition, Species, Order, Status, Days, A)
data_alive <- data_alive %>% filter(A > 0)  
data_alive <- as.data.frame(expandRows(data_alive, "A"))
rownames(data_alive) <- NULL
data_alive$event <- 0

data_survival <- rbind(data_death, data_alive)
data_survival$Order_Status <- paste0(data_survival$Order, "_", data_survival$Status)
data_survival <- data_survival %>% select(Island, Year, Condition, Species, Order, Status, Order_Status, Days, event)


##------
## LT50
##------
tab.lt50 <- matrix(nrow = 0, ncol = 7)


for (i in 1:length(unique(data_survival$Condition))){
  sub1_data_survival <- data_survival %>% filter(Condition == unique(data_survival$Condition)[i])
  
  for (j in 1:length(unique(data_survival$Island))){
    sub2_data_survival <- sub1_data_survival %>% filter(Island == unique(sub1_data_survival$Island)[j])
    
    surv.species <- survreg(Surv(Days, event) ~ Species, data = sub2_data_survival)
    mortalite = 0.50
    
    for (k in 1:length(unique(sub2_data_survival$Species))){
      sub3_data_survival <- sub2_data_survival %>% filter(Species == unique(sub2_data_survival$Species)[k]) 
      
      sub.tab.lt50 <- matrix(nrow = 1, ncol = 7)
    
      p.surv.species <- predict(surv.species, newdata = list(Species = unique(sub2_data_survival$Species)[k]), type = "quantile", p = mortalite, se = TRUE)
      
      sub.tab.lt50[1, 1] <- unique(sub3_data_survival$Island)
      sub.tab.lt50[1, 2] <- unique(sub3_data_survival$Condition)
      sub.tab.lt50[1, 3] <- unique(sub3_data_survival$Species)
      sub.tab.lt50[1, 4] <- unique(sub3_data_survival$Order)
      sub.tab.lt50[1, 5] <- unique(sub3_data_survival$Status)
      sub.tab.lt50[1, 6] <- p.surv.species$fit[1]
      sub.tab.lt50[1, 7] <- 1.96*p.surv.species$se.fit[1]
    
      tab.lt50 <- rbind(tab.lt50, sub.tab.lt50)
    
    }    
  }
}

rm(i, j, k, sub1_data_survival, sub2_data_survival, sub3_data_survival, sub.tab.lt50, surv.species, mortalite, p.surv.species)
colnames(tab.lt50) <- c("Island","Condition", "Cd.species", "Order", "Status", "lt50","95ci.lt50")
write.table(tab.lt50,"outputs/LT50.txt")
