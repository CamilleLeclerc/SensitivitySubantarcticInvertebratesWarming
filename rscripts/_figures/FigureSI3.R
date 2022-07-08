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


species.color <- as.data.frame(sort(unique(data_survival$Species)))
#scales::show_col(colormap(colormap = colormaps$inferno, nshades = 7))
#scales::show_col(colormap(colormap = colormaps$viridis, nshades = 4))
species.color$color <- c("#6b176cff", "#400154ff", "#23978aff", "#79d051ff", "#df5636ff", "#f9940fff", "#2a0d4fff", "#ab305bff", "#000004ff", "#39558bff")
colnames(species.color)[1] <- "species"

order.color <- as.data.frame(sort(unique(data_survival$Order_Status)))
order.color$color <- c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854")
colnames(order.color)[1] <- "order.status"


##------------------------------------------------
## FIT SURVIVAL DATA USING THE KAPLAN-MEIER METHOD
##------------------------------------------------
##Pairwise survdiff
res.order <- pairwise_survdiff(Surv(Days, event) ~ Island + Order + Status + Condition, data = data_survival)
res.order

res.order.pvalue <- res.order$p.value
res.order.pvalue[0.05 < res.order.pvalue] <- 1
res.order.pvalue[0.01 < res.order.pvalue & res.order.pvalue <= 0.05] <- 0.05
res.order.pvalue[0.001 < res.order.pvalue & res.order.pvalue <= 0.01] <- 0.01
res.order.pvalue[res.order.pvalue <= 0.001] <- 0.001
#Signif. codes:  0 ’***’ 0.001 ’**’ 0.01 ’*’ 0.05 ’.’ 0.1 ’ ’ 1

melted.order <- melt(res.order.pvalue, na.rm = TRUE)
melted.order$value <- as.factor(melted.order$value)

melted.order <- data.frame(lapply(melted.order, function(x) { gsub("Island=", "", x)}))
melted.order <- data.frame(lapply(melted.order, function(x) { gsub("Order=", "", x)}))
melted.order <- data.frame(lapply(melted.order, function(x) { gsub("Status=", "", x)}))
melted.order <- data.frame(lapply(melted.order, function(x) { gsub("Condition=", "", x)}))
melted.order <- data.frame(lapply(melted.order, function(x) { gsub(",", "-", x)}))
melted.order <- data.frame(lapply(melted.order, function(x) { gsub(" ", "", x)}))


ggplot(data = melted.order, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("#d73027", "#fc8d59", "#fee090", "white")) +
  coord_fixed() +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))+
  
  geom_segment(aes(x = 0.5, y = 0.5, xend = 0.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 1.5, y = 0.5, xend = 1.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 2.5, y = 1.5, xend = 2.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 3.5, y = 2.5, xend = 3.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 4.5, y = 3.5, xend = 4.5, yend = 23.5), data = melted.order)+ 
  geom_segment(aes(x = 5.5, y = 4.5, xend = 5.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 6.5, y = 5.5, xend = 6.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 7.5, y = 6.5, xend = 7.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 8.5, y = 7.5, xend = 8.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 9.5, y = 8.5, xend = 9.5, yend = 23.5), data = melted.order)+ 
  geom_segment(aes(x = 10.5, y = 9.5, xend = 10.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 11.5, y = 10.5, xend = 11.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 12.5, y = 11.5, xend = 12.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 13.5, y = 12.5, xend = 13.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 14.5, y = 13.5, xend = 14.5, yend = 23.5), data = melted.order)+ 
  geom_segment(aes(x = 15.5, y = 14.5, xend = 15.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 16.5, y = 15.5, xend = 16.5, yend = 23.5), data = melted.order)+ 
  geom_segment(aes(x = 17.5, y = 16.5, xend = 17.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 18.5, y = 17.5, xend = 18.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 19.5, y = 18.5, xend = 19.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 20.5, y = 19.5, xend = 20.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 21.5, y = 20.5, xend = 21.5, yend = 23.5), data = melted.order)+ 
  geom_segment(aes(x = 22.5, y = 21.5, xend = 22.5, yend = 23.5), data = melted.order)+
  geom_segment(aes(x = 23.5, y = 22.5, xend = 23.5, yend = 23.5), data = melted.order)+
  
  geom_segment(aes(x = 0.5, y = 0.5, xend = 1.5, yend = 0.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 1.5, xend = 2.5, yend = 1.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 2.5, xend = 3.5, yend = 2.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 3.5, xend = 4.5, yend = 3.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 4.5, xend = 5.5, yend = 4.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 5.5, xend = 6.5, yend = 5.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 6.5, xend = 7.5, yend = 6.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 7.5, xend = 8.5, yend = 7.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 8.5, xend = 9.5, yend = 8.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 9.5, xend = 10.5, yend = 9.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 10.5, xend = 11.5, yend = 10.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 11.5, xend = 12.5, yend = 11.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 12.5, xend = 13.5, yend = 12.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 13.5, xend = 14.5, yend = 13.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 14.5, xend = 15.5, yend = 14.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 15.5, xend = 16.5, yend = 15.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 16.5, xend = 17.5, yend = 16.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 17.5, xend = 18.5, yend = 17.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 18.5, xend = 19.5, yend = 18.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 19.5, xend = 20.5, yend = 19.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 20.5, xend = 21.5, yend = 20.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 21.5, xend = 22.5, yend = 21.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 22.5, xend = 23.5, yend = 22.5), data = melted.order)+
  geom_segment(aes(x = 0.5, y = 23.5, xend = 23.5, yend = 23.5), data = melted.order)
#save pdf 8x8  
