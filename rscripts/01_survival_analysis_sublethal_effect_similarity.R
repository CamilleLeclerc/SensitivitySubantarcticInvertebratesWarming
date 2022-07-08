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

data_DSL<- data %>% select(Island, Year, Condition, Species, Order, Status, Days, D_SL)
data_DSL <- data_DSL %>% filter(D_SL > 0)  
data_DSL <- as.data.frame(expandRows(data_DSL, "D_SL"))
rownames(data_DSL) <- NULL
data_DSL$event <- 1

data_A <- data %>% select(Island, Year, Condition, Species, Order, Status, Days, A)
data_A <- data_A %>% filter(A > 0)  
data_A <- as.data.frame(expandRows(data_A, "A"))
rownames(data_A) <- NULL
data_A$event <- 0

data_survival_1 <- rbind(data_DSL, data_A)
data_survival_1$sublethal_effect <- "D_SL"

##-------------
data_D <- data %>% select(Island, Year, Condition, Species, Order, Status, Days, D)
data_D <- data_D %>% filter(D > 0)  
data_D <- as.data.frame(expandRows(data_D, "D"))
rownames(data_D) <- NULL
data_D$event <- 1

data_survival_2 <- rbind(data_D, data_A)
data_survival_2$sublethal_effect <- "D"

data_survival <- rbind(data_survival_1, data_survival_2)
rm(data_DSL, data_A, data_D, data_survival_1, data_survival_2)


##-----------------
## AMALO - C2 - KER
##-----------------
amalo_c2 <- data_survival %>% filter(Condition == "C2" & Island == "KER" & Species == "AMALO")
km_amalo_c2 <- survfit(Surv(Days, event) ~ sublethal_effect, data = amalo_c2)
km_amalo_c2
survdiff(Surv(Days, event) ~ sublethal_effect, data = amalo_c2)
ggsurvplot(km_amalo_c2)
surv_pvalue(km_amalo_c2)$pval.txt


##-----------------
## ANATK - C1 - KER
##-----------------
anatk_c1 <- data_survival %>% filter(Condition == "C1" & Island == "KER" & Species == "ANATK")
km_anatk_c1 <- survfit(Surv(Days, event) ~ sublethal_effect, data = anatk_c1)
km_anatk_c1
survdiff(Surv(Days, event) ~ sublethal_effect, data = anatk_c1)
ggsurvplot(km_anatk_c1)


##-----------------
## ANATK - C2 - KER
##-----------------
anatk_c2 <- data_survival %>% filter(Condition == "C2" & Island == "KER" & Species == "ANATK")
km_anatk_c2 <- survfit(Surv(Days, event) ~ sublethal_effect, data = anatk_c2)
km_anatk_c2
survdiff(Surv(Days, event) ~ sublethal_effect, data = anatk_c2)
ggsurvplot(km_anatk_c2)


##-----------------
## ANATK - C4 - KER
##-----------------
anatk_c4 <- data_survival %>% filter(Condition == "C4" & Island == "KER" & Species == "ANATK")
km_anatk_c4 <- survfit(Surv(Days, event) ~ sublethal_effect, data = anatk_c4)
km_anatk_c4
survdiff(Surv(Days, event) ~ sublethal_effect, data = anatk_c4)
ggsurvplot(km_anatk_c4)


##-----------------
## CACHX - C1 - KER
##-----------------
cachx_c1 <- data_survival %>% filter(Condition == "C1" & Island == "KER" & Species == "CACHX")
km_cachx_c1 <- survfit(Surv(Days, event) ~ sublethal_effect, data = cachx_c1)
km_cachx_c1
survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c1)
ggsurvplot(km_cachx_c1)


##-----------------
## CACHX - C2 - KER
##-----------------
cachx_c2 <- data_survival %>% filter(Condition == "C2" & Island == "KER" & Species == "CACHX")
km_cachx_c2 <- survfit(Surv(Days, event) ~ sublethal_effect, data = cachx_c2)
km_cachx_c2
survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c2)
ggsurvplot(km_cachx_c2)


##-----------------
## CACHX - C3 - KER
##-----------------
cachx_c3 <- data_survival %>% filter(Condition == "C3" & Island == "KER" & Species == "CACHX")
km_cachx_c3 <- survfit(Surv(Days, event) ~ sublethal_effect, data = cachx_c3)
km_cachx_c3
survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c3)
ggsurvplot(km_cachx_c3)


##-----------------
## CACHX - C4 - KER
##-----------------
cachx_c4 <- data_survival %>% filter(Condition == "C4" & Island == "KER" & Species == "CACHX")
km_cachx_c4 <- survfit(Surv(Days, event) ~ sublethal_effect, data = cachx_c4)
km_cachx_c4
survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c4)
ggsurvplot(km_cachx_c4)


##-----------------
## CALAI - C1 - KER
##-----------------
calai_c1 <- data_survival %>% filter(Condition == "C1" & Island == "KER" & Species == "CALAI")
km_calai_c1 <- survfit(Surv(Days, event) ~ sublethal_effect, data = calai_c1)
km_calai_c1
survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c1)
ggsurvplot(km_calai_c1)


##-----------------
## CALAI - C2 - KER
##-----------------
calai_c2 <- data_survival %>% filter(Condition == "C2" & Island == "KER" & Species == "CALAI")
km_calai_c2 <- survfit(Surv(Days, event) ~ sublethal_effect, data = calai_c2)
km_calai_c2
survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c2)
ggsurvplot(km_calai_c2)


##-----------------
## CALAI - C3 - KER
##-----------------
calai_c3 <- data_survival %>% filter(Condition == "C3" & Island == "KER" & Species == "CALAI")
km_calai_c3 <- survfit(Surv(Days, event) ~ sublethal_effect, data = calai_c3)
km_calai_c3
survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c3)
ggsurvplot(km_calai_c3)


##-----------------
## CALAI - C4 - KER
##-----------------
calai_c4 <- data_survival %>% filter(Condition == "C4" & Island == "KER" & Species == "CALAI")
km_calai_c4 <- survfit(Surv(Days, event) ~ sublethal_effect, data = calai_c4)
km_calai_c4
survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c4)
ggsurvplot(km_calai_c4)


##----------------
## ECVI - C4 - KER
##----------------
ecvi_c4 <- data_survival %>% filter(Condition == "C4" & Island == "KER" & Species == "ECVI")
km_ecvi_c4 <- survfit(Surv(Days, event) ~ sublethal_effect, data = ecvi_c4)
km_ecvi_c4
survdiff(Surv(Days, event) ~ sublethal_effect, data = ecvi_c4)
ggsurvplot(km_ecvi_c4)


##-----------------
## AMBP - C3 - CRO
##-----------------
ambp_c4 <- data_survival %>% filter(Condition == "C3" & Island == "CRO" & Species == "AMBP")
km_ambp_c4 <- survfit(Surv(Days, event) ~ sublethal_effect, data = ambp_c4)
km_ambp_c4
survdiff(Surv(Days, event) ~ sublethal_effect, data = ambp_c4)
ggsurvplot(km_ambp_c4)


##-----------------
## ANATC - C3 - CRO
##-----------------
anatc_c3 <- data_survival %>% filter(Condition == "C3" & Island == "CRO" & Species == "ANATC")
km_anatc_c3 <- survfit(Surv(Days, event) ~ sublethal_effect, data = anatc_c3)
km_anatc_c3
survdiff(Surv(Days, event) ~ sublethal_effect, data = anatc_c3)
ggsurvplot(km_anatc_c3)


##-----------------
## MYROK - C1 - CRO
##-----------------
myrok_c1 <- data_survival %>% filter(Condition == "C1" & Island == "CRO" & Species == "MYROKC")
km_myrok_c1 <- survfit(Surv(Days, event) ~ sublethal_effect, data = myrok_c1)
km_myrok_c1
survdiff(Surv(Days, event) ~ sublethal_effect, data = myrok_c1)
ggsurvplot(km_myrok_c1)


##-----------------
## MYROK - C3 - CRO
##-----------------
myrok_c3 <- data_survival %>% filter(Condition == "C3" & Island == "CRO" & Species == "MYROKC")
km_myrok_c3 <- survfit(Surv(Days, event) ~ sublethal_effect, data = myrok_c3)
km_myrok_c3
survdiff(Surv(Days, event) ~ sublethal_effect, data = myrok_c3)
ggsurvplot(km_myrok_c3)


##-----------------
## MYROK - C4 - CRO
##-----------------
myrok_c4 <- data_survival %>% filter(Condition == "C4" & Island == "CRO" & Species == "MYROKC")
km_myrok_c4 <- survfit(Surv(Days, event) ~ Year, data = myrok_c4)
km_myrok_c4
survdiff(Surv(Days, event) ~ sublethal_effect, data = myrok_c4)
ggsurvplot(km_myrok_c4)


##----------------
## STATISTIC TABLE
##----------------
statistic.table <- as.data.frame(matrix(NA, nrow = 18, ncol = 3, dimnames = list(NULL, c("Modality", "Chisq", "pvalue"))))
statistic.table$Modality <- c("amalo_c2", "anatk_c1", "anatk_c2", "anatk_c4", "cachx_c1", "cachx_c2", "cachx_c3", "cachx_c4", "calai_c1", "calai_c2", "calai_c3", "calai_c4", "ecvi_c4", "ambp_c3", "anatc_c3", "myrok_c1", "myrok_c3", "myrok_c4")
statistic.table[1, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = amalo_c2)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = amalo_c2)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = amalo_c2)$n)-1))
statistic.table[2, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = anatk_c1)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = anatk_c1)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = anatk_c1)$n)-1))
statistic.table[3, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = anatk_c2)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = anatk_c2)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = anatk_c2)$n)-1))
statistic.table[4, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = anatk_c4)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = anatk_c4)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = anatk_c4)$n)-1))
statistic.table[5, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c1)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c1)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c1)$n)-1))
statistic.table[6, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c2)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c2)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c2)$n)-1))
statistic.table[7, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c3)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c3)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c3)$n)-1))
statistic.table[8, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c4)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c4)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = cachx_c4)$n)-1))
statistic.table[9, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c1)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c1)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c1)$n)-1))
statistic.table[10, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c2)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c2)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c2)$n)-1))
statistic.table[11, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c3)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c3)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c3)$n)-1))
statistic.table[12, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c4)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c4)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = calai_c4)$n)-1))
statistic.table[13, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = ecvi_c4)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = ecvi_c4)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = ecvi_c4)$n)-1))
statistic.table[14, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = ambp_c4)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = ambp_c4)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = ambp_c4)$n)-1))
statistic.table[15, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = anatc_c3)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = anatc_c3)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = anatc_c3)$n)-1))
statistic.table[16, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = myrok_c1)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = myrok_c1)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = myrok_c1)$n)-1))
statistic.table[17, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = myrok_c3)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = myrok_c3)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = myrok_c3)$n)-1))
statistic.table[18, 2:3] <- c(survdiff(Surv(Days, event) ~ sublethal_effect, data = myrok_c4)$chisq , 1-pchisq(survdiff(Surv(Days, event) ~ sublethal_effect, data = myrok_c4)$chisq, length(survdiff(Surv(Days, event) ~ sublethal_effect, data = myrok_c4)$n)-1))


netpvalues <- as.data.frame(statistic.table$pvalue)
rownames(netpvalues) <- netpvalues$Modality
colnames(netpvalues) <- c("pvalue")
(sum(p.adjust(netpvalues$pvalue, method = "bonferroni") >= 0.05)/(length(netpvalues$pvalue)))*100
