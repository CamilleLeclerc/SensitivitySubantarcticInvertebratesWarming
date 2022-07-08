rm(list=ls())




##--------------
## LOAD PACKAGES
##--------------
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




## ANATK - C1 - KER
anatk_c1 <- data_survival %>% filter(Condition == "C1" & Island == "KER" & Species == "ANATK" & Days <= 10)

## ANATK - C2 - KER
anatk_c2 <- data_survival %>% filter(Condition == "C2" & Island == "KER" & Species == "ANATK" & Days <= 10)

## ANATK - C3 - KER
anatk_c3 <- data_survival %>% filter(Condition == "C3" & Island == "KER" & Species == "ANATK")

## ANATK - C4 - KER
anatk_c4 <- data_survival %>% filter(Condition == "C4" & Island == "KER" & Species == "ANATK" & Days <= 10)

## AMALO - C3 - KER
amalo_c3 <- data_survival %>% filter(Condition == "C3" & Island == "KER" & Species == "AMALO" & Days <= 10)

## CACHX - C1 - KER
cachx_c1 <- data_survival %>% filter(Condition == "C1" & Island == "KER" & Species == "CACHX" & Days <= 10)

## CACHX - C2 - KER
cachx_c2 <- data_survival %>% filter(Condition == "C2" & Island == "KER" & Species == "CACHX" & Days <= 10)

## CACHX - C3 - KER
cachx_c3 <- data_survival %>% filter(Condition == "C3" & Island == "KER" & Species == "CACHX" & Days <= 10)

## CACHX - C4 - KER
cachx_c4 <- data_survival %>% filter(Condition == "C4" & Island == "KER" & Species == "CACHX" & Days <= 10)

## CALAI - C1 - KER
calai_c1 <- data_survival %>% filter(Condition == "C1" & Island == "KER" & Species == "CALAI" & Days <= 10)

## CALAI - C2 - KER
calai_c2 <- data_survival %>% filter(Condition == "C2" & Island == "KER" & Species == "CALAI" & Days)

## CALAI - C3 - KER
calai_c3 <- data_survival %>% filter(Condition == "C3" & Island == "KER" & Species == "CALAI" & Days <= 10)

## CALAI - C4 - KER
calai_c4 <- data_survival %>% filter(Condition == "C4" & Island == "KER" & Species == "CALAI" & Days <= 10)

## ANATK - C1 - CRO
anatkc_c1 <- data_survival %>% filter(Condition == "C1" & Island == "CRO" & Species == "ANATK" & Days <= 10)

## ANATK - C2 - CRO
anatkc_c2 <- data_survival %>% filter(Condition == "C2" & Island == "CRO" & Species == "ANATK" & Days <= 10)

## ANATK - C3 - CRO
anatkc_c3 <- data_survival %>% filter(Condition == "C3" & Island == "CRO" & Species == "ANATK" & Days <= 10)

## ANATK - C4 - CRO
anatkc_c4 <- data_survival %>% filter(Condition == "C4" & Island == "CRO" & Species == "ANATK" & Days <= 10)

## AMBP - C4 - CRO
ambp_c4 <- data_survival %>% filter(Condition == "C4" & Island == "CRO" & Species == "AMBP")

## MYROK - C3 - CRO
myrok_c3 <- data_survival %>% filter(Condition == "C3" & Island == "CRO" & Species == "MYROKC" & Days <= 10)

## MYROK - C4 - CRO
myrok_c4 <- data_survival %>% filter(Condition == "C4" & Island == "CRO" & Species == "MYROKC" & Days <= 16)




##----------------
## STATISTIC TABLE
##----------------
statistic.table <- as.data.frame(matrix(NA, nrow = 20, ncol = 3, dimnames = list(NULL, c("Modality", "Chisq", "pvalue"))))
statistic.table$Modality <- c("anatk_c1", "anatk_c2", "anatk_c3", "anatk_c4", "amalo_c3", "cachx_c1", "cachx_c2", "cachx_c3", "cachx_c4", "calai_c1", "calai_c2", "calai_c3", "calai_c4", "anatkc_c1", "anatkc_c2", "anatkc_c3", "anatkc_c4", "ambp_c4", "myrok_c3", "myrok_c4")
statistic.table[1, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = anatk_c1)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = anatk_c1)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = anatk_c1)$n)-1))
statistic.table[2, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = anatk_c2)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = anatk_c2)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = anatk_c2)$n)-1))
statistic.table[3, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = anatk_c3)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = anatk_c3)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = anatk_c3)$n)-1))
statistic.table[4, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = anatk_c4)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = anatk_c4)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = anatk_c4)$n)-1))
statistic.table[5, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = amalo_c3)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = amalo_c3)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = amalo_c3)$n)-1))
statistic.table[6, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = cachx_c1)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = cachx_c1)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = cachx_c1)$n)-1))
statistic.table[7, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = cachx_c2)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = cachx_c2)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = cachx_c2)$n)-1))
statistic.table[8, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = cachx_c3)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = cachx_c3)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = cachx_c3)$n)-1))
statistic.table[9, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = cachx_c4)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = cachx_c4)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = cachx_c4)$n)-1))
statistic.table[10, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = calai_c1)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = calai_c1)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = calai_c1)$n)-1))
statistic.table[11, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = calai_c2)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = calai_c2)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = calai_c2)$n)-1))
statistic.table[12, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = calai_c3)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = calai_c3)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = calai_c3)$n)-1))
statistic.table[13, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = calai_c4)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = calai_c4)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = calai_c4)$n)-1))
statistic.table[14, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = anatkc_c1)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = anatkc_c1)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = anatkc_c1)$n)-1))
statistic.table[15, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = anatkc_c2)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = anatkc_c2)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = anatkc_c2)$n)-1))
statistic.table[16, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = anatkc_c3)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = anatkc_c3)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = anatkc_c3)$n)-1))
statistic.table[17, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = anatkc_c4)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = anatkc_c4)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = anatkc_c4)$n)-1))
statistic.table[18, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = ambp_c4)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = ambp_c4)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = ambp_c4)$n)-1))
statistic.table[19, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = myrok_c3)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = myrok_c3)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = myrok_c3)$n)-1))
statistic.table[20, 2:3] <- c(survdiff(Surv(Days, event) ~ Year, data = myrok_c4)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ Year, data = myrok_c4)$chisq, length(survdiff(Surv(Days, event) ~ Year, data = myrok_c4)$n)-1))


netpvalues <- as.data.frame(statistic.table$pvalue)
rownames(netpvalues) <- netpvalues$Modality
colnames(netpvalues) <- c("pvalue")
(sum(p.adjust(netpvalues$pvalue, method = "bonferroni") >= 0.05)/(length(netpvalues$pvalue)))*100
