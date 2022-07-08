##---------
## TUTORIAL
##---------
#https://mctavishlab.github.io/R_OpenTree_tutorials/index.html


rm(list=ls())


##--------------
## LOAD PACKAGES
##--------------
library(ape)
library(datelife)
library(datelifeplot)
library(devtools)
library(ggpubr)
library(httr)
library(mulTree)
library(nlme)
library(phytools)
library(rotl)
library(stringr)


##-----------------------------------------------
## FINDING TAXA IN THE OPEN TREE OF LIFE TAXONOMY
##-----------------------------------------------
my_taxa <- c("Amalopteryx maritima",
             "Amblystogenium pacificum",
             "Anatalanta aptera",
             "Anatalanta crozetensis",
             "Calycopteryx mosleyi",
             "Ectemnorhinus viridis",
             "Fucellia tergina",
             "Merizodus soledadinus",
             "Myro kerguelenensis")
resolved_names <- rotl::tnrs_match_names(names = my_taxa, context_name = "All life")
resolved_names

my_ott_ids <- ott_id(resolved_names)

OpenTreeOfLife <- read.tree("data/labelled_supertree/labelled_supertree.tre")
str(OpenTreeOfLife)
class(OpenTreeOfLife)
tip.labels.list <- as.data.frame(OpenTreeOfLife$tip.label, stringsAsFactors = TRUE) ; colnames(tip.labels.list) <- "tip.label"
nrow(tip.labels.list)
tip.labels.list <- as.data.frame(tip.labels.list[tip.labels.list$tip.label != "ott4424258" &
                                                         tip.labels.list$tip.label != "ott5040473" &
                                                         tip.labels.list$tip.label != "ott1081185" &
                                                         tip.labels.list$tip.label != "ott4398130" &
                                                         tip.labels.list$tip.label != "ott4410377" &
                                                         tip.labels.list$tip.label != "ott556821" &
                                                         tip.labels.list$tip.label != "ott4360222" &
                                                         tip.labels.list$tip.label != "ott679489" &
                                                         tip.labels.list$tip.label != "ott3560180" ,])
nrow(tip.labels.list)
my_tree <- drop.tip(OpenTreeOfLife, as.character(tip.labels.list$tip.label), root.edge = 0)
plot(my_tree)
my_tree$tip.label
my_tree$tip.label <- c("Fucellia_tergina",
                       "Amalopteryx_maritima",
                       "Anatalanta_aptera",
                       "Anatalanta_crozetensis",
                       "Calycopteryx_mosleyi",
                       "Ectemnorhinus_viridis",
                       "Merizodus_soledadinus",
                       "Amblystogenium_pacificum",
                       "Myro_kerguelenensis")
plot(my_tree)
my_tree <- compute.brlen(my_tree)
vcv(my_tree)
str(my_tree)



##------------------------
## DATA LETHAL TEMPERATURE
##------------------------
LT50 <- read.csv("outputs/LT50.txt", sep="")
summary(LT50)
LT50[5, 6:7] <- NA
LT50[9, 6:7] <- NA
LT50 <- LT50[complete.cases(LT50), ]
LT50 <- LT50[is.finite(LT50$lt50),]


LT50$Species[LT50$Cd.species == "AMALO"] <- "Amalopteryx_maritima"
LT50$Species[LT50$Cd.species == "AMBP"] <- "Amblystogenium_pacificum"
LT50$Species[LT50$Cd.species == "ANATC"] <- "Anatalanta_crozetensis"
LT50$Species[LT50$Cd.species == "ANATK"] <- "Anatalanta_aptera"
LT50$Species[LT50$Cd.species == "CACHX"] <- "Calycopteryx_mosleyi"
LT50$Species[LT50$Cd.species == "CALAI"] <- "Calycopteryx_mosleyi"
LT50$Species[LT50$Cd.species == "ECVI"] <- "Ectemnorhinus_viridis"
LT50$Species[LT50$Cd.species == "FUCE"] <- "Fucellia_tergina"
LT50$Species[LT50$Cd.species == "MERIZ"] <- "Merizodus_soledadinus"
LT50$Species[LT50$Cd.species == "MYROKC"] <- "Myro_kerguelenensis"


ggdensity(LT50$lt50, 
          main = "Density plot of LT50", xlab = "Time to 50% survival (LT50, in days) of species under the warming scenarios")
ggqqplot(LT50$lt50)
shapiro.test(LT50$lt50) #From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.

#https://www.r-bloggers.com/2020/01/a-guide-to-data-transformation/
LT50$log_lt50 <- scale(log10(LT50$lt50))
ggdensity(LT50$log_lt50, main = "Density plot of log(LT50)", xlab = "Log of time to 50% survival (LT50, in days) of species under the warming scenarios")
ggqqplot(LT50$log_lt50)
shapiro.test(LT50$log_lt50) #From the output, the p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution. In other words, we can assume the normality.

LT50 <- LT50[, c("Species", "Condition", "Island", "Order", "Status", "log_lt50")]




##---------
## MCMCglmm
##---------
#https://stackoverflow.com/questions/51132259/phylogenetic-model-using-multiple-entries-for-each-species
#https://github.com/TGuillerme/mulTree/blob/master/doc/Vanilla_flavoured_phylogenetic_analyses.Rmd

## Creates a mulTree object specifying species as random terms
mulTree_data <- as.mulTree(LT50, my_tree, taxa = "Species", rand.terms = ~Species)

## The glmm formula
formula <- log_lt50 ~ Condition + Island + Order + Status
#formula <- log_lt50 ~ Condition * Island * Order * Status

## The MCMC parameters (number of generations, thin/sampling, burnin)
mcmc_parameters <- c(100000, 10, 250)

## The MCMCglmm priors
mcmc_priors <- list(R = list(V = 1, nu = 0.002),
                    G = list(G1 = list(V = 1, nu = 0.002))) 

## Running MCMCglmm on multiple trees
mulTree(mulTree_data, formula = formula, parameters = mcmc_parameters,
        priors = mcmc_priors, chains = 2, output = "log_lt50")


model <- read.mulTree("log_lt50-tree1_chain1", model = TRUE)
str(model)
summary(model, use.hdr = FALSE)
#summary(model, use.hdr = FALSE, cent.tend = mean, prob = c(5,95))

model.final <- read.mulTree("log_lt50-tree1_chain1", model = FALSE)
str(model.final)
#summary(model.final, use.hdr = FALSE)
summary(model.final, use.hdr = FALSE, cent.tend = mean, prob = c(5,95))

plot(summary(model.final, use.hdr = FALSE))
plot(summary(model.final, use.hdr = FALSE),
     horizontal = TRUE,
     ylab = "",
     cex.coeff = 0.8,
     main = "Posterior distributions", 
     ylim = c(-3,3),
     cex.terms = 0.5,
     terms = c("Intercept",
               "Condition - C2", "Condition - C3", "Condition - C4",
               "Island - Kerguelen",
               "Order - Coleoptera", "Order - Diptera",
               "Status - Non-native",
               "Phylogeny", "Residuals"),
     cex.main = 0.8)
abline(v = 0, lty = 3)

