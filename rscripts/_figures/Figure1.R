rm(list=ls())


##--------------
## LOAD PACKAGES
##--------------
library(ggplot2)
library(ggthemes)


##-----
## DATA
##-----
experimental_design <- read.delim("data/experimental_design.txt")
str(experimental_design)
experimental_design$Hour <- as.factor(experimental_design$Hour)
sort(experimental_design$Hour)
experimental_design$Hour <- factor(experimental_design$Hour, levels = c("0:00",
                                                                        "1:00", "2:00", "3:00", "4:00", 
                                                                        "5:00", "6:00", "7:00", "8:00",
                                                                        "9:00", "10:00", "11:00", "12:00", 
                                                                        "13:00", "14:00", "15:00", "16:00", 
                                                                        "17:00", "18:00", "19:00", "20:00", 
                                                                        "21:00", "22:00", "23:00", "24:00"))


##---------
## FIGURE 1
##---------
ggplot(experimental_design, aes(x = Hour, y = Temperature, group=1)) +
  geom_line(size = 1.5)+
  geom_point(size = 1.5) +
  facet_wrap(~ Condition) +
  scale_y_continuous(breaks = seq(4, 28, len = 4)) +
  scale_x_discrete(labels = c("0:00",
                              "", "", "", "4:00", 
                              "", "", "", "8:00",
                              "", "", "", "12:00", 
                              "", "", "", "16:00", 
                              "", "", "", "20:00", 
                              "", "", "", "24:00")) +
  xlab("Hours") + ylab("Temperature (°C)") +
  theme_bw() +
  theme(axis.text = element_text(size = 12, colour="black"),
        axis.title = element_text(size = 16, colour="black", face="bold"),
        axis.line.x = element_line(color="black"),
        axis.line.y = element_line(color="black"),
        strip.text.x = element_text(size = 16, colour="black", face="bold"))
