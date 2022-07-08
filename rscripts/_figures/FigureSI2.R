rm(list=ls())


##--------------
## LOAD PACKAGES
##--------------
library(colormap)
library(GGally)
library(gtsummary)
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
##Model fitting and plot
for (i in 1:length(unique(data_survival$Condition))){
  sub1_data_survival <- data_survival %>% filter(Condition == unique(data_survival$Condition)[i])
  
  for (j in 1:length(unique(data_survival$Island))){
    sub2_data_survival <- sub1_data_survival %>% filter(Island == unique(sub1_data_survival$Island)[j])   
    
    surv.order <- survreg(Surv(Days, event) ~ Order_Status, data = sub2_data_survival)
    summary(surv.order)
    
    
    #pdf(file = paste0("pdf/Order_", unique(sub2_data_survival$Condition), "_", unique(sub2_data_survival$Island),".pdf"),   # The directory you want to save the file in
    #    width = 4, # The width of the plot in inches
    #    height = 4)
    
    plot(1, 1, xlim = c(0,70), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "",
         main = paste("Condition", unique(sub2_data_survival$Condition), unique(sub2_data_survival$Island), sep = " "),
         cex.lab = 1, cex.main = 1, cex.sub = 1)
    axis(1, mgp = c(3, 2, 0), cex.axis = 1)
    axis(2, las = 2, mgp = c(3, 1.2, 0), cex.axis = 1)
    
    
    for (k in 1:length(unique(sub2_data_survival$Order_Status))){
      p.surv.order <- predict(surv.order, newdata = list(Order_Status = unique(sub2_data_survival$Order_Status)[k]), type = "quantile", p = seq(0.01, 1, by = 0.01), se = TRUE)
      lines(predict(surv.order, newdata = list(Order_Status = unique(sub2_data_survival$Order_Status)[k]), type = "quantile", p = seq(0,1,0.01)), 1-seq(0,1,0.01), lwd = 3, col = unlist(order.color %>% filter(order.status == unique(sub2_data_survival$Order_Status)[k]) %>% select(color)))
      lines(predict(surv.order, newdata = list(Order_Status = unique(sub2_data_survival$Order_Status)[k]), type = "quantile", p = seq(0,1,0.01))- 1.96*p.surv.order$se.fit, 1-seq(0,1,0.01), lty = 3, col = unlist(order.color %>% filter(order.status == unique(sub2_data_survival$Order_Status)[k]) %>% select(color))) 
      lines(predict(surv.order, newdata = list(Order_Status = unique(sub2_data_survival$Order_Status)[k]), type = "quantile", p = seq(0,1,0.01))+ 1.96*p.surv.order$se.fit, 1-seq(0,1,0.01), lty = 3, col = unlist(order.color %>% filter(order.status == unique(sub2_data_survival$Order_Status)[k]) %>% select(color)))
      legend(60, 1-as.numeric(paste0(0,".", k)), 
             legend = unique(sub2_data_survival$Order_Status)[k],
             col = unlist(order.color %>% filter(order.status == unique(sub2_data_survival$Order_Status)[k]) %>% select(color)),
             lty = 1,
             lwd = 2)
    }
    #dev.off()
  }
}
rm(i, j, k, sub1_data_survival, sub2_data_survival, surv.order, p.surv.order)
