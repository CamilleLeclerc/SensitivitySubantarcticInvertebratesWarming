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
##The Kaplan-Meier estimator - Species // Island + Condition
#The Kaplan-Meier estimator is a non-parametric statistic that allows us to estimate the survival function.
km.species <- survfit(Surv(Days, event) ~ Species + Island + Condition, data = data_survival)
km.species
summary(km.species)
survdiff(Surv(Days, event) ~ Species + Island + Condition, data = data_survival)
ggsurvplot(km.species)
ggsurv <- ggsurvplot(km.species, conf.int = TRUE, ggtheme = theme_bw())
ggsurv$plot +
  theme_bw() + 
  theme (legend.position = "right") +
  facet_grid(Island ~ Condition)


##Model fitting and plot
for (i in 1:length(unique(data_survival$Condition))){
  sub1_data_survival <- data_survival %>% filter(Condition == unique(data_survival$Condition)[i])
  
  for (j in 1:length(unique(data_survival$Island))){
    sub2_data_survival <- sub1_data_survival %>% filter(Island == unique(sub1_data_survival$Island)[j])    
    
    surv.species <- survreg(Surv(Days, event) ~ Species, data = sub2_data_survival)
    summary(surv.species)
    
    #pdf(file = paste0("pdf/Species_", unique(sub2_data_survival$Condition), "_", unique(sub2_data_survival$Island),".pdf"),   # The directory you want to save the file in
    #    width = 4, # The width of the plot in inches
    #    height = 4)
    
    plot(1, 1, xlim = c(0,70), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "",
          main = paste("Condition", unique(sub2_data_survival$Condition), unique(sub2_data_survival$Island), sep = " "),
          cex.lab = 1, cex.main = 1, cex.sub = 1)
    axis(1, mgp = c(3, 2, 0), cex.axis = 1)
    axis(2, las = 2, mgp = c(3, 1.2, 0), cex.axis = 1)

    
    for (k in 1:length(unique(sub2_data_survival$Species))){
      p.surv.species <- predict(surv.species, newdata = list(Species = unique(sub2_data_survival$Species)[k]), type = "quantile", p = seq(0.01, 1, by = 0.01), se = TRUE)
      lines(predict(surv.species, newdata = list(Species = unique(sub2_data_survival$Species)[k]), type = "quantile", p = seq(0,1,0.01)), 1-seq(0,1,0.01), lwd = 3, col = unlist(species.color %>% filter(species == unique(sub2_data_survival$Species)[k]) %>% select(color)))
      lines(predict(surv.species, newdata = list(Species = unique(sub2_data_survival$Species)[k]), type = "quantile", p = seq(0,1,0.01))- 1.96*p.surv.species$se.fit, 1-seq(0,1,0.01), lty = 3, col = unlist(species.color %>% filter(species == unique(sub2_data_survival$Species)[k]) %>% select(color))) 
      lines(predict(surv.species, newdata = list(Species = unique(sub2_data_survival$Species)[k]), type = "quantile", p = seq(0,1,0.01))+ 1.96*p.surv.species$se.fit, 1-seq(0,1,0.01), lty = 3, col = unlist(species.color %>% filter(species == unique(sub2_data_survival$Species)[k]) %>% select(color)))
      legend(60, 1-as.numeric(paste0(0,".", k)), 
             legend = unique(sub2_data_survival$Species)[k],
             col = unlist(species.color %>% filter(species == unique(sub2_data_survival$Species)[k]) %>% select(color)),
             lty = 1,
             lwd = 2)
          }
    #dev.off()
  }
}
rm(i, j, k, sub1_data_survival, sub2_data_survival, surv.species, p.surv.species)


##Residuals checking
for (i in 1:length(unique(data_survival$Condition))){
  sub1_data_survival <- data_survival %>% filter(Condition == unique(data_survival$Condition)[i])
  
  for (j in 1:length(unique(data_survival$Island))){
    sub2_data_survival <- sub1_data_survival %>% filter(Island == unique(sub1_data_survival$Island)[j])    
    
    for (k in 1:length(unique(sub2_data_survival$Species))){
      sub3_data_survival <- sub2_data_survival %>% filter(Species == unique(sub2_data_survival$Species)[k])    
    
      surv.species <- survreg(Surv(Days, event) ~ 1, data = sub3_data_survival)
    
    fitted_values <- surv.species$linear.predictors
    resids <- (log(surv.species$y[, 1]) - fitted_values) / surv.species$scale
    resKM <- survfit(Surv(resids, event) ~ 1, data = sub3_data_survival)
    plot(resKM, mark.time = FALSE, xlab = "AFT Residuals", ylab = "Survival Probability")
    xx <- seq(min(resids), max(resids), length.out = 35)
    yy <- exp(- exp(xx))
    lines(xx, yy, col = "red", lwd = 2)
    legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                       "Survival function of Extreme Value distribution"), 
            lty = c(1,2,1), col = c(1,1,2), bty = "n")
    }
  }
}
rm(i, j, k, sub1_data_survival, sub2_data_survival, sub3_data_survival, surv.species, fitted_values, resids, resKM, xx, yy)


##Pairwise survdiff
res.species <- pairwise_survdiff(Surv(Days, event) ~  Island + Species + Condition, data = data_survival)
res.species

res.species.pvalue <- res.species$p.value
res.species.pvalue[0.05 < res.species.pvalue] <- 1
res.species.pvalue[0.01 < res.species.pvalue & res.species.pvalue <= 0.05] <- 0.05
res.species.pvalue[0.001 < res.species.pvalue & res.species.pvalue <= 0.01] <- 0.01
res.species.pvalue[res.species.pvalue <= 0.001] <- 0.001
#Signif. codes:  0 ’***’ 0.001 ’**’ 0.01 ’*’ 0.05 ’.’ 0.1 ’ ’ 1

melted.species <- melt(res.species.pvalue, na.rm = TRUE)
melted.species$value <- as.factor(melted.species$value)

melted.species <- data.frame(lapply(melted.species, function(x) { gsub("Island=", "", x)}))
melted.species <- data.frame(lapply(melted.species, function(x) { gsub("Species=", "", x)}))
melted.species <- data.frame(lapply(melted.species, function(x) { gsub("Condition=", "", x)}))
melted.species <- data.frame(lapply(melted.species, function(x) { gsub(",", "-", x)}))
melted.species <- data.frame(lapply(melted.species, function(x) { gsub(" ", "", x)}))


ggplot(data = melted.species, aes(Var2, Var1, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("#d73027", "#fc8d59", "#fee090", "white")) +
  coord_fixed() +
  theme(panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill = NA))+
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), legend.title = element_blank()) +
  theme(axis.text.x = element_text(angle = 90))+
 
  geom_segment(aes(x = 0.5, y = 0.5, xend = 0.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 1.5, y = 0.5, xend = 1.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 2.5, y = 1.5, xend = 2.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 3.5, y = 2.5, xend = 3.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 4.5, y = 3.5, xend = 4.5, yend = 39.5), data = melted.species)+ 
  geom_segment(aes(x = 5.5, y = 4.5, xend = 5.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 6.5, y = 5.5, xend = 6.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 7.5, y = 6.5, xend = 7.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 8.5, y = 7.5, xend = 8.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 9.5, y = 8.5, xend = 9.5, yend = 39.5), data = melted.species)+ 
  geom_segment(aes(x = 10.5, y = 9.5, xend = 10.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 11.5, y = 10.5, xend = 11.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 12.5, y = 11.5, xend = 12.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 13.5, y = 12.5, xend = 13.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 14.5, y = 13.5, xend = 14.5, yend = 39.5), data = melted.species)+ 
  geom_segment(aes(x = 15.5, y = 14.5, xend = 15.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 16.5, y = 15.5, xend = 16.5, yend = 39.5), data = melted.species)+ 
  geom_segment(aes(x = 17.5, y = 16.5, xend = 17.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 18.5, y = 17.5, xend = 18.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 19.5, y = 18.5, xend = 19.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 20.5, y = 19.5, xend = 20.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 21.5, y = 20.5, xend = 21.5, yend = 39.5), data = melted.species)+ 
  geom_segment(aes(x = 22.5, y = 21.5, xend = 22.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 23.5, y = 22.5, xend = 23.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 24.5, y = 23.5, xend = 24.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 25.5, y = 24.5, xend = 25.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 26.5, y = 25.5, xend = 26.5, yend = 39.5), data = melted.species)+ 
  geom_segment(aes(x = 27.5, y = 26.5, xend = 27.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 28.5, y = 27.5, xend = 28.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 29.5, y = 28.5, xend = 29.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 30.5, y = 29.5, xend = 30.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 31.5, y = 30.5, xend = 31.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 32.5, y = 31.5, xend = 32.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 33.5, y = 32.5, xend = 33.5, yend = 39.5), data = melted.species)+ 
  geom_segment(aes(x = 34.5, y = 33.5, xend = 34.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 35.5, y = 34.5, xend = 35.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 36.5, y = 35.5, xend = 36.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 37.5, y = 36.5, xend = 37.5, yend = 39.5), data = melted.species)+
  geom_segment(aes(x = 38.5, y = 37.5, xend = 38.5, yend = 39.5), data = melted.species)+ 
  geom_segment(aes(x = 39.5, y = 38.5, xend = 39.5, yend = 39.5), data = melted.species)+
  
  geom_segment(aes(x = 0.5, y = 0.5, xend = 1.5, yend = 0.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 1.5, xend = 2.5, yend = 1.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 2.5, xend = 3.5, yend = 2.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 3.5, xend = 4.5, yend = 3.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 4.5, xend = 5.5, yend = 4.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 5.5, xend = 6.5, yend = 5.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 6.5, xend = 7.5, yend = 6.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 7.5, xend = 8.5, yend = 7.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 8.5, xend = 9.5, yend = 8.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 9.5, xend = 10.5, yend = 9.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 10.5, xend = 11.5, yend = 10.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 11.5, xend = 12.5, yend = 11.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 12.5, xend = 13.5, yend = 12.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 13.5, xend = 14.5, yend = 13.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 14.5, xend = 15.5, yend = 14.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 15.5, xend = 16.5, yend = 15.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 16.5, xend = 17.5, yend = 16.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 17.5, xend = 18.5, yend = 17.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 18.5, xend = 19.5, yend = 18.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 19.5, xend = 20.5, yend = 19.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 20.5, xend = 21.5, yend = 20.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 21.5, xend = 22.5, yend = 21.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 22.5, xend = 23.5, yend = 22.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 23.5, xend = 24.5, yend = 23.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 24.5, xend = 25.5, yend = 24.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 25.5, xend = 26.5, yend = 25.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 26.5, xend = 27.5, yend = 26.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 27.5, xend = 28.5, yend = 27.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 28.5, xend = 29.5, yend = 28.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 29.5, xend = 30.5, yend = 29.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 30.5, xend = 31.5, yend = 30.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 31.5, xend = 32.5, yend = 31.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 32.5, xend = 33.5, yend = 32.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 33.5, xend = 34.5, yend = 33.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 34.5, xend = 35.5, yend = 34.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 35.5, xend = 36.5, yend = 35.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 36.5, xend = 37.5, yend = 36.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 37.5, xend = 38.5, yend = 37.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 38.5, xend = 39.5, yend = 38.5), data = melted.species)+
  geom_segment(aes(x = 0.5, y = 39.5, xend = 39.5, yend = 39.5), data = melted.species)
#save pdf 8x8  
 
 


##-----------------------------------------------------------------------------------------------
##-----------------------------------------------------------------------------------------------

##The Kaplan-Meier estimator - Order // Island + Condition
#The Kaplan-Meier estimator is a non-parametric statistic that allows us to estimate the survival function.
km.order <- survfit(Surv(Days, event) ~ Order + Island + Condition, data = data_survival)
km.order
summary(km.order)
survdiff(Surv(Days, event) ~ Order + Island + Condition, data = data_survival)
ggsurvplot(km.order)
ggsurv <- ggsurvplot(km.order, conf.int = TRUE, ggtheme = theme_bw())
ggsurv$plot +
  theme_bw() + 
  theme (legend.position = "right") +
  facet_grid(Island ~ Condition)


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


##Residuals checking
for (i in 1:length(unique(data_survival$Condition))){
  sub1_data_survival <- data_survival %>% filter(Condition == unique(data_survival$Condition)[i])
  
  for (j in 1:length(unique(data_survival$Island))){
    sub2_data_survival <- sub1_data_survival %>% filter(Island == unique(sub1_data_survival$Island)[j])    
    
    for (k in 1:length(unique(sub2_data_survival$Order))){
      sub3_data_survival <- sub2_data_survival %>% filter(Order == unique(sub2_data_survival$Order)[k])    
      
      surv.order <- survreg(Surv(Days, event) ~ 1, data = sub3_data_survival)
      
      fitted_values <- surv.order$linear.predictors
      resids <- (log(surv.order$y[, 1]) - fitted_values) / surv.order$scale
      resKM <- survfit(Surv(resids, event) ~ 1, data = sub3_data_survival)
      plot(resKM, mark.time = FALSE, xlab = "AFT Residuals", ylab = "Survival Probability")
      xx <- seq(min(resids), max(resids), length.out = 35)
      yy <- exp(- exp(xx))
      lines(xx, yy, col = "red", lwd = 2)
      legend("bottomleft", c("KM estimate", "95% CI KM estimate", 
                             "Survival function of Extreme Value distribution"), 
             lty = c(1,2,1), col = c(1,1,2), bty = "n")
    }
  }
}
rm(i, j, k, sub1_data_survival, sub2_data_survival, sub3_data_survival, surv.order, fitted_values, resids, resKM, xx, yy)


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