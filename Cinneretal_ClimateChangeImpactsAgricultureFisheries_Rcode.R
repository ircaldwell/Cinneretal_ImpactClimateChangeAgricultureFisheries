#### Project: Cinner et al. The potential impacts of climate change on agriculture and fisheries production in 73 tropical coastal communities ####
#### Objective : Analysis and creating figures ####
#### Authors : Lauric Thiault & Iain R. Caldwell ####
#### Last update : 2022/02/01 ####
####  Contents:
####  1) Analysis - calculate statistics for the main text
####    a) Model agreement
####    b) Fisheries and agriculture exposure comparison
####    c) Fisheries and agriculture sensitivity comparison 
####    d) Within country exposure and sensitivity variability
####    e) Quantifying double burden (sites that experience losses in fisheries and agriculture) 
####    f) Potential impact (Euclidean distance of sensitivity and exposure) vs. Material style of life (MSL)
####    g) Sensitivity (combined agriculture and fisheries) vs. Material style of life (MSL)
####    h) Exposure (combined agriculture and fisheries) vs. MSL
####  2) Main figures
####    a) Figure 1: Exposure vs. sensitivity with model agreement for agriculture and fisheries with map of study sites
####    b) Figure 2: Agriculture vs. fisheries exposure and comparison with random locations
####    c) Figure 3: Agriculture vs. fisheries exposure with sector dependency and across SSP's
####    d) Figure 4: Potential impacts vs. MSL across SSP's
####    e) Figure 5: Sensitivity and material style of life over time in two PNG locations
####  3) Supplementary figures
####    a) Figure S1: Projected agricultural changes by crop
####    b) Figure S2: Agriculture vs. fisheries exposure with random sites compared between SSP's
####    c) Figure S3: Combined exposure vs combined sensitivity across SSPs
####    d) Figure S4: Combined sensitivity and combined exposure vs. MSL
####    e) Figure S5: Map of locations from which agriculture and fisheries data were extracted 
####    f) Figure S9: Finer scale maps showing locations of study sites in each country with average model agreement

# Remove variables from memory
rm(list = ls())
dev.off()

# Load packages
library(tidyverse)
library(ggplot2)
library(lme4) #mixed models
library(lmerTest) #extend results of mixed models to include t and p-values
library(effsize) #for Cohen's D
library(scales) #for rescaling
library(car) #for Anova function
library(MuMIn) #for r.squaredGLMM function
#remotes::install_github("RemkoDuursma/bootpredictlme4")
library(bootpredictlme4)
library(RColorBrewer)
library(egg) #for combining plots
library(lemon) # for facet_rep_wrap() function
library(ggtext) #for geom_richtext
library(ggpubr) #to combine plots within plots
library(rgdal)

#### Set directories ####
dataDir <- "Data/"
plotDir <- "Plots/"

#### Import data ####
agFishSiteTBL <- read_csv(file = paste0(dataDir, "Cinneretal_ImpactClimateChangeAgricultureFisheries_SiteNumData.csv")) #Site level data
agFishRandTBL <- read_csv(file = paste0(dataDir, "Cinneretal_ImpactClimateChangeAgricultureFisheries_RandomSiteData.csv")) #Random test data
pngTimeTBL <- read_csv(file = paste0(dataDir, "Cinneretal_ImpactClimateChangeAgricultureFisheries_TemporalPNGdata.csv")) #Sensitivity and MSL data over time for two sites in PNG
mslPcaResTBL <- read_csv(file = paste0(dataDir, "Cinneretal_ImpactClimateChangeAgricultureFisheries_MSL_PCAresults.csv")) #Results from principal component analysis for MSL
agFishExpLocsTBL <- read_csv(file = paste0(dataDir, "Cinneretal_ImpactClimateChangeAgricultureFisheries_ExpCells.csv")) #Locations neighboring communities from which exposure data was extracted

#Create folder for plots if it does not exist
if(!file.exists(plotDir)){
  dir.create(plotDir)
}

####  1) Analysis - calculate statistics for the main text ####
# Convert exposure data to show magnitude of LOSS (i.e. change the sign of exposure)
agFishSiteTBL <- agFishSiteTBL %>% 
  mutate_at(vars(matches("E_")), ~.*-1) %>% 
  mutate_at(vars(matches("_PercModelAgree")), ~.*-1) #Convert back the % model agreement

agFishRandTBL <- agFishRandTBL %>% 
  mutate_at(vars(matches("E_")), ~.*-1) %>% 
  mutate_at(vars(matches("_PercModelAgree")), ~.*-1) #Convert back the % model agreement

# Calculate stats for main text
std_errFn <- function(x) sd(x)/sqrt(length(x))

####    a) Model agreement ####
fishModelAgreeSsp585Lmer <- lmer(data = agFishSiteTBL, formula = E_fish_Ssp585_PercModelAgree ~ 1 + (1 | Country), REML = F)
summary(fishModelAgreeSsp585Lmer) #Average model agreement = 84.7 +/- 4.50 after accounting for country

fishModelAgreeSsp126Lmer <- lmer(data = agFishSiteTBL, formula = E_fish_Ssp126_PercModelAgree ~ 1 + (1 | Country), REML = F)
summary(fishModelAgreeSsp126Lmer) #Average model agreement = 89.2 +/- 4.06 after accounting for country

agModelAgreeSsp585Lmer <- lmer(data = agFishSiteTBL, formula = E_agr_Ssp585_PercModelAgree ~ 1 + (1 | Country), REML = F)
summary(agModelAgreeSsp585Lmer) #Average model agreement = 69.1 +/- 4.82 after accounting for country

agModelAgreeSsp126Lmer <- lmer(data = agFishSiteTBL, formula = E_agr_Ssp126_PercModelAgree ~ 1 + (1 | Country), REML = F)
summary(agModelAgreeSsp126Lmer) #Average model agreement = 70.4 +/- 3.27 after accounting for country

####    b) Fisheries and agriculture exposure comparison ####
###   Study sites
##  Average and SE of fisheries losses - study data - SSP5 - 8.5
fishLossesSsp585Lmer <- lmer(data = agFishSiteTBL, formula = E_fish_Ssp585 ~ 1 + (1 | Country), REML = F)
fishLossesSsp585MeanEst <- signif(coef(summary(fishLossesSsp585Lmer))["(Intercept)", "Estimate"], digits = 3) #14.7% loss
fishLossesSsp585SE <- signif(coef(summary(fishLossesSsp585Lmer))["(Intercept)", "Std. Error"], digits = 3) #+/- 4.27%

##  Average and SE of agriculture losses - study data - SSP5 - 8.5
agLossesSsp585Lmer <- lmer(data = agFishSiteTBL, formula = E_agr_Ssp585_Unweighted ~ 1 + (1 | Country), REML = F)
agLossesSsp585MeanEst <- signif(coef(summary(agLossesSsp585Lmer))["(Intercept)", "Estimate"], digits = 3) #-1.17% - gains
agLossesSsp585SE <- signif(coef(summary(agLossesSsp585Lmer))["(Intercept)", "Std. Error"], digits = 3) #+/- 1.47%

##  Average and SE of fisheries losses - study data - SSP1 - 2.6
fishLossesSsp126Lmer <- lmer(data = agFishSiteTBL, formula = E_fish_Ssp126 ~ 1 + (1 | Country), REML = F)
fishLossesSsp126MeanEst <- signif(coef(summary(fishLossesSsp126Lmer))["(Intercept)", "Estimate"], digits = 3) #11.2% loss
fishLossesSsp126SE <- signif(coef(summary(fishLossesSsp126Lmer))["(Intercept)", "Std. Error"], digits = 3) #+/- 2.71%

##  Average and SE of agriculture losses - study data - SSP1 - 2.6
agLossesSsp126Lmer <- lmer(data = agFishSiteTBL, formula = E_agr_Ssp126_Unweighted ~ 1 + (1 | Country), REML = F)
agLossesSsp126MeanEst <- signif(coef(summary(agLossesSsp126Lmer))["(Intercept)", "Estimate"], digits = 3) #-1.06%
agLossesSsp126SE <- signif(coef(summary(agLossesSsp126Lmer))["(Intercept)", "Std. Error"], digits = 3) #+/- 0.736%

###   Random sites
##  Average and SE of fisheries losses - random sites - SSP5 - 8.5
fishLossesRandSsp585Lmer <- lmer(data = agFishRandTBL, formula = E_fish_Ssp585 ~ 1 + (1 | Country), REML = F)
fishLossesRandSsp585MeanEst <- signif(coef(summary(fishLossesRandSsp585Lmer))["(Intercept)", "Estimate"], digits = 3) #14.6% loss
fishLossesRandSsp585SE <- signif(coef(summary(fishLossesRandSsp585Lmer))["(Intercept)", "Std. Error"], digits = 3) #+/- 3.85%

##  Average and SE of agriculture losses - random sites - SSP5 - 8.5
agLossesRandSsp585Lmer <- lmer(data = agFishRandTBL, formula = E_agr_Ssp585_Unweighted ~ 1 + (1 | Country), REML = F)
agLossesRandSsp585MeanEst <- signif(coef(summary(agLossesRandSsp585Lmer))["(Intercept)", "Estimate"], digits = 3) #-1.01%
agLossesRandSsp585SE <- signif(coef(summary(agLossesRandSsp585Lmer))["(Intercept)", "Std. Error"], digits = 3) #+/- 1.38%

##  Average and SE of fisheries losses - random sites - SSP1 - 2.6
fishLossesRandSsp126Lmer <- lmer(data = agFishRandTBL, formula = E_fish_Ssp126 ~ 1 + (1 | Country), REML = F)
fishLossesRandSsp126MeanEst <- signif(coef(summary(fishLossesRandSsp126Lmer))["(Intercept)", "Estimate"], digits = 3) #11.7% loss
fishLossesRandSsp126SE <- signif(coef(summary(fishLossesRandSsp126Lmer))["(Intercept)", "Std. Error"], digits = 3) #+/- 2.29%

##  Average and SE of agriculture losses - random sites - SSP1 - 2.6
agLossesRandSsp126Lmer <- lmer(data = agFishRandTBL, formula = E_agr_Ssp126_Unweighted ~ 1 + (1 | Country), REML = F)
agLossesRandSsp126MeanEst <- signif(coef(summary(agLossesRandSsp126Lmer))["(Intercept)", "Estimate"], digits = 3) #-0.895%
agLossesRandSsp126SE <- signif(coef(summary(agLossesRandSsp126Lmer))["(Intercept)", "Std. Error"], digits = 3) #+/- 0.666%

### Production-weighted agriculture losses
##  Average and SE of agriculture losses - weighted study data - SSP5 - 8.5
weightedAgLossesSsp585Lmer <- lmer(data = agFishSiteTBL, formula = E_agr_Ssp585_Weighted ~ 1 + (1 | Country), REML = F)
summary(weightedAgLossesSsp585Lmer) #-0.741 +/- 1.39% 

##Comparing weighted and unweighted data
agFishSiteTBL <- agFishSiteTBL %>% 
  mutate(WeightedMinusUnweighted_E_agr_Ssp585 = E_agr_Ssp585_Weighted - E_agr_Ssp585_Unweighted)

weightVsUnweightDiffLmer <- lmer(agFishSiteTBL,
                                 formula = WeightedMinusUnweighted_E_agr_Ssp585 ~ 1 + (1 | Country), REML = F)
summary(weightVsUnweightDiffLmer) #0.433 +/- 

##  Comparing the fisheries and agriculture losses
agFishSiteTBL <- agFishSiteTBL %>%
  mutate(E_fishMinAgr_Ssp585 = E_fish_Ssp585 - E_agr_Ssp585_Unweighted)

fishVsAgrDiffLmer <- lmer(data = agFishSiteTBL, formula = E_fishMinAgr_Ssp585 ~ 1 + (1 | Country), REML = F)
summary(fishVsAgrDiffLmer) #Fisheries losses are 15.9 +/- 5.65% greater than agricultural losses (t = 2.81; p = 0.0379)

##  Comparing random sites - fisheries and agriculture losses
agFishRandTBL <- agFishRandTBL %>% 
  mutate(E_fishMinAgr_Ssp585 = E_fish_Ssp585 - E_agr_Ssp585_Unweighted)

randFishVsAgrDiffLmer <- lmer(data = agFishRandTBL, formula = E_fishMinAgr_Ssp585 ~ 1 + (1 | Country), REML = F)
summary(randFishVsAgrDiffLmer) #Fisheries losses are 15.6 +/- 5.12% greater than agricultural losses (t = 3.06; p = 0.0282)

##  Cohen's D comparing the study sites and random sites
#SSP5 - 8.5
cohen.d(agFishSiteTBL$E_agr_Ssp585_Unweighted, agFishRandTBL$E_agr_Ssp585_Unweighted) #d = -0.308 ("small" difference in agriculture exposure)
cohen.d(agFishSiteTBL$E_fish_Ssp585, agFishRandTBL$E_fish_Ssp585) #d = -0.0216 ("negligible" difference in fisheries exposure)

#SSP1 - 2.6
cohen.d(agFishSiteTBL$E_agr_Ssp126_Unweighted, agFishRandTBL$E_agr_Ssp126_Unweighted) #d = -0.352 ("small" difference in agriculture exposure)
cohen.d(agFishSiteTBL$E_fish_Ssp126, agFishRandTBL$E_fish_Ssp126) #d = -0.0296 ("negligible" difference in fisheries exposure)

####    c) Fisheries and agriculture sensitivity comparison #### 
##  Average and SE of fisheries sensitivity
fishSensLmer <- lmer(data = agFishSiteTBL, formula = S_fsh ~ 1 + (1 | Country), REML = F)
summary(fishSensLmer) #Estimate = 0.0772; SE = 0.00748

##  Average and SE of agriculture sensitivity
agSensLmer <- lmer(data = agFishSiteTBL, formula = S_agr ~ 1 + (1 | Country), REML = F)
summary(agSensLmer) #Estimate = 0.0402; SE = 0.00964

##  Comparing the fisheries and agriculture sensitivity
agFishSiteTBL <- agFishSiteTBL %>% 
  mutate(S_fishMinAgr = S_fsh - S_agr)

fishVsAgrSensDiffLmer <- lmer(data = agFishSiteTBL, formula = S_fishMinAgr ~ 1 + (1 | Country), REML = F)
summary(fishVsAgrSensDiffLmer) #Fisheries sensitivity is 0.0410 +/- 0.0136 greater than agricultural sensitivity (t = 3.01; p = 0.0815)

####    d) Within country exposure, sensitivity, and model agreement variability #### 
## Average and variability of fisheries exposure for Indonesian sites
mean(agFishSiteTBL$E_fish_Ssp585[agFishSiteTBL$Country == "indonesia"]) #16.5% loss
std_errFn(agFishSiteTBL$E_fish_Ssp585[agFishSiteTBL$Country == "indonesia"]) #2.21%
range(agFishSiteTBL$E_fish_Ssp585[agFishSiteTBL$Country == "indonesia"]) #6.54 to 32.0%

## Variability of fisheries exposure and sensitivity for Philippines sites
range(agFishSiteTBL$E_fish_Ssp585[agFishSiteTBL$Country == "philippines"]) #8.90 - 12.6% loss
range(agFishSiteTBL$S_fsh[agFishSiteTBL$Country == "philippines"]) #0.000580 - 0.319

## Variability of model agreement for Indonesian sites
range(agFishSiteTBL$E_agr_Ssp585_PercModelAgree[agFishSiteTBL$Country == "indonesia"]) #50 - 85% agreement
range(agFishSiteTBL$E_fish_Ssp585_PercModelAgree[agFishSiteTBL$Country == "indonesia"]) #56 - 100% agreement
range(agFishSiteTBL$E_agr_Ssp126_PercModelAgree[agFishSiteTBL$Country == "indonesia"]) #50 - 80% agreement
range(agFishSiteTBL$E_fish_Ssp126_PercModelAgree[agFishSiteTBL$Country == "indonesia"]) #50 - 94% agreement


####    e) Quantifying double burden (sites that experience losses in fisheries and agriculture) #### 
agFishSiteTBL <- agFishSiteTBL %>% 
  mutate(DoubleBurden_Ssp585 = ifelse(test = E_agr_Ssp585_Unweighted > 0 & E_fish_Ssp585 > 0, yes = "yes", no = "no"),
         DoubleBurden_Ssp126 = ifelse(test = E_agr_Ssp126_Unweighted > 0 & E_fish_Ssp126 > 0, yes = "yes", no = "no"))

agFishRandTBL <- agFishRandTBL %>% 
  mutate(DoubleBurden_Ssp585 = ifelse(test = E_agr_Ssp585_Unweighted > 0 & E_fish_Ssp585 > 0, yes = "yes", no = "no"),
         DoubleBurden_Ssp126 = ifelse(test = E_agr_Ssp126_Unweighted > 0 & E_fish_Ssp126 > 0, yes = "yes", no = "no"))

## Percentage of study sites that will lose in both fisheries and agriculture under SSP5 - 8.5
sum(agFishSiteTBL$DoubleBurden_Ssp585 == "yes")/nrow(agFishSiteTBL)*100 #63.9%

## Percentage of study sites that will lose in both fisheries and agriculture under SSP1 - 2.6
sum(agFishSiteTBL$DoubleBurden_Ssp126 == "yes")/nrow(agFishSiteTBL)*100 #37.5%

## Percentage of random sites that will lose in both fisheries and agriculture under SSP5 - 8.5
sum(agFishRandTBL$DoubleBurden_Ssp585 == "yes")/nrow(agFishRandTBL)*100 #70.4%

## Percentage of random sites that will lose in both fisheries and agriculture under SSP1 - 2.6
sum(agFishRandTBL$DoubleBurden_Ssp126 == "yes")/nrow(agFishRandTBL)*100 #47.6%

## Range of all cumulative sensitivity values vs. those for sites with double burden
range(agFishSiteTBL$S_agrfsh) #0.00332 to 0.403
range(agFishSiteTBL$S_agrfsh[agFishSiteTBL$DoubleBurden_Ssp585 == "yes"]) #0.00332 to 0.339

## Model agreement in those sites with double burden vs. other sites
agrModelAgreeDoubleBurdLmer <- lmer(data = agFishSiteTBL, formula = E_agr_Ssp585_PercModelAgree ~ DoubleBurden_Ssp585 + (1 | Country), REML = F)
summary(agrModelAgreeDoubleBurdLmer) #those with double burden have 25.8 +/- 2.78% less agriculture model agreement (p < 0.001)

fishModelAgreeDoubleBurdLmer <- lmer(data = agFishSiteTBL, formula = E_fish_Ssp585_PercModelAgree ~ DoubleBurden_Ssp585 + (1 | Country), REML = F)
summary(fishModelAgreeDoubleBurdLmer) #those with double burden have 20.3 +/- 3.93% more fisheries model agreement (p < 0.001)

#Calculate average model agreement (between agriculture and fisheries)
agFishSiteTBL <- agFishSiteTBL %>% 
  mutate(E_agrfsh_Ssp126_AvgPercModelAgreement = (E_agr_Ssp126_PercModelAgree + E_fish_Ssp126_PercModelAgree)/2,
         E_agrfsh_Ssp585_AvgPercModelAgreement = (E_agr_Ssp585_PercModelAgree + E_fish_Ssp585_PercModelAgree)/2)

agfishModelAgreeDoubleBurdLmer <- lmer(data = agFishSiteTBL, formula = E_agrfsh_Ssp585_AvgPercModelAgreement ~ DoubleBurden_Ssp585 + (1 | Country), REML = F)
summary(agfishModelAgreeDoubleBurdLmer) #those with double burden have 5.15 +/- 1.97% less average model agreement (p = 0.0149)

####    f) Potential impact (Euclidean distance of sensitivity and exposure) vs. Material style of life (MSL) #### 
###     Calculate "potential impact"
##        Calculate the average loss (exposure from agr and fish) and the direction of the loss
agFishSiteTBL <- agFishSiteTBL %>% 
  mutate(E_agrfsh_Ssp126 = (E_agr_Ssp126_Unweighted + E_fish_Ssp126)/2,
         E_agrfsh_Ssp585 = (E_agr_Ssp585_Unweighted + E_fish_Ssp585)/2)

##Make the data longer by turning the combined exposure from SSPs into one column
##  Then rescale the exposure, sensitivity, and MSL data
##  Then get the Euclidean distance from the origin for S and E to get "potential impact"
agFishSiteLongTBL <- agFishSiteTBL %>% 
  dplyr::select(SiteNum, Country, S_agrfsh, MSL1, E_agrfsh_Ssp126, E_agrfsh_Ssp585) %>% 
  gather(key = "SSP", value = "E_agrfsh", "E_agrfsh_Ssp126":"E_agrfsh_Ssp585") %>% 
  mutate(SSP = dplyr::recode(SSP,
                             E_agrfsh_Ssp126 = "SSP1 - 2.6",
                             E_agrfsh_Ssp585 = "SSP5 - 8.5"),
         S_agrfsh_rescale = scales::rescale(S_agrfsh),
         MSL_rescale = scales::rescale(MSL1),
         E_agrfsh_rescale = scales::rescale(abs(E_agrfsh)),
         PotImp = scales::rescale(sqrt(S_agrfsh_rescale^2 + E_agrfsh_rescale^2)))

####Fit mixed models to get the slope, p, and R2 values
fit26 <- lmer(data = agFishSiteLongTBL %>% filter(SSP == "SSP1 - 2.6"),
              formula = PotImp ~ MSL_rescale + (1|Country))
summary(fit26) #Estimate (slope) = -0.530 +/- 0.107

fit85 <- lmer(data = agFishSiteLongTBL %>% filter(SSP == "SSP5 - 8.5"),
              formula = PotImp ~ MSL_rescale + (1|Country))
summary(fit85) #Estimate (slope) = -0.485 +/- 0.138

# extract p-values
pvalue26 <- coef(summary(fit26))["MSL_rescale", "Pr(>|t|)"] #0.00764
pvalue85 <- coef(summary(fit85))["MSL_rescale", "Pr(>|t|)"] #0.0352

# extract R squared
r2_26 <- r.squaredGLMM(fit26) #R2m = 0.385; R2c = 0.442
r2_85 <- r.squaredGLMM(fit85) #R2m = 0.261; R2c = 0.352

# extract 95% confidence intervals for each value
pred26 <- predict(fit26, newdata = data.frame(MSL_rescale = agFishSiteLongTBL$MSL_rescale[agFishSiteLongTBL$SSP == "SSP1 - 2.6"]),
                  re.form = NA, se.fit = TRUE, nsim = 1000)
pred85 <- predict(fit85, newdata = data.frame(MSL_rescale = agFishSiteLongTBL$MSL_rescale[agFishSiteLongTBL$SSP == "SSP5 - 8.5"]),
                  re.form = NA, se.fit = TRUE, nsim = 1000)

#Assign values back to the data
agFishSiteLongTBL <- agFishSiteLongTBL %>% 
  mutate(Pvalue = ifelse(test = SSP == "SSP1 - 2.6", yes = pvalue26, no = pvalue85),
         R2m = ifelse(test = SSP == "SSP1 - 2.6", yes = r2_26[1], no = r2_85[1]),
         R2c = ifelse(test = SSP == "SSP1 - 2.6", yes = r2_26[2], no = r2_85[2]),
         fit = as.numeric(NA),
         seBoot = as.numeric(NA))

agFishSiteLongTBL$fit[agFishSiteLongTBL$SSP == "SSP1 - 2.6"] <- pred26$fit
agFishSiteLongTBL$fit[agFishSiteLongTBL$SSP == "SSP5 - 8.5"] <- pred85$fit

agFishSiteLongTBL$seBoot[agFishSiteLongTBL$SSP == "SSP1 - 2.6"] <- pred26$se.boot
agFishSiteLongTBL$seBoot[agFishSiteLongTBL$SSP == "SSP5 - 8.5"] <- pred85$se.boot

####    g) Sensitivity (combined agriculture and fisheries) vs. Material style of life (MSL) #### 
####      Run lmer models for the sensitvity vs. MSL
agFishSiteTBL <- agFishSiteTBL %>% 
  mutate(S_agrfsh_rescale = rescale(S_agrfsh),
         MSL_rescale = rescale(MSL1))

####Fit mixed models to get the slope, R2 and p values for S vs. MSL
fitSvsMSL <- lmer(data = agFishSiteTBL,
                  formula = S_agrfsh_rescale ~ MSL_rescale + (1|Country))
summary(fitSvsMSL) #Estimate = -0.715 +/- 0.174

# extract p-values
pvalueSvsMSL <- coef(summary(fitSvsMSL))["MSL_rescale", "Pr(>|t|)"] #0.00107

# extract R squared
r2_SvsMSL <- r.squaredGLMM(fitSvsMSL) #R2m = 0.407; R2c = 0.683

# extract 95%confidence intervals
predSvsMSL <- predict(fitSvsMSL,
                      newdata = data.frame(MSL_rescale = agFishSiteTBL$MSL_rescale),
                      re.form = NA, se.fit = TRUE, nsim = 1000)

#Assign these back to the data
agFishSiteTBL <- agFishSiteTBL %>% 
  mutate(fitSvsMSL = predSvsMSL$fit,
         seBootSvsMSL = predSvsMSL$se.boot)

####    h) Exposure (combined agriculture and fisheries) vs. Material style of life (MSL) #### 
####      Fit mixed models to get the slope, R2 and p values for Exposure vs. MSL for RCP 2.6 and 8.5 - with different models for gains and losses
fit26_ExpVsMSL <- lmer(data = agFishSiteLongTBL %>% filter(SSP == "SSP1 - 2.6"),
                       formula = E_agrfsh_rescale ~ MSL_rescale + (1|Country))
summary(fit26_ExpVsMSL) #Estimate (slope) = -0.349 +/- 0.163

fit85_ExpVsMSL <- lmer(data = agFishSiteLongTBL %>% filter(SSP == "SSP5 - 8.5"),
                       formula = E_agrfsh_rescale ~ MSL_rescale + (1|Country))
summary(fit85_ExpVsMSL) #Estimate (slope) = -0.480 +/- 0.236

# extract p-values
pvalue26_ExpVsMSL <- coef(summary(fit26_ExpVsMSL))["MSL_rescale", "Pr(>|t|)"] #0.0486
pvalue85_ExpVsMSL <- coef(summary(fit85_ExpVsMSL))["MSL_rescale", "Pr(>|t|)"] #0.0491

# extract R squared
r2_26_ExpVsMSL <- r.squaredGLMM(fit26_ExpVsMSL) #R2m = 0.157; R2c = 0.538
r2_85_ExpVsMSL <- r.squaredGLMM(fit85_ExpVsMSL) #R2m = 0.122; R2c = 0.697

# extract 95%confidence intervals
pred26_ExpVsMSL <- predict(fit26_ExpVsMSL,
                           newdata = data.frame(MSL_rescale = agFishSiteLongTBL$MSL_rescale[agFishSiteLongTBL$SSP == "SSP1 - 2.6"]),
                           re.form = NA, se.fit = TRUE, nsim = 1000)

pred85_ExpVsMSL <- predict(fit85_ExpVsMSL,
                           newdata = data.frame(MSL_rescale = agFishSiteLongTBL$MSL_rescale[agFishSiteLongTBL$SSP == "SSP5 - 8.5"]),
                           re.form = NA, se.fit = TRUE, nsim = 1000)

#Assign these back to the data
agFishSiteLongTBL <- agFishSiteLongTBL %>% 
  mutate(Pvalue_ExpVsMSL = ifelse(test = SSP == "SSP1 - 2.6", yes = pvalue26_ExpVsMSL, no = pvalue85_ExpVsMSL),
         Rm2_ExpVsMSL = ifelse(test = SSP == "SSP1 - 2.6", yes = r2_26_ExpVsMSL[1], no = r2_85_ExpVsMSL[1]),
         Rc2_ExpVsMSL = ifelse(test = SSP == "SSP1 - 2.6", yes = r2_26_ExpVsMSL[2], no = r2_85_ExpVsMSL[2]),
         fit_ExpVsMSL = NA,
         seBoot_ExpVsMSL = NA)

agFishSiteLongTBL$fit_ExpVsMSL[agFishSiteLongTBL$SSP == "SSP1 - 2.6"] <- pred26_ExpVsMSL$fit
agFishSiteLongTBL$fit_ExpVsMSL[agFishSiteLongTBL$SSP == "SSP5 - 8.5"] <- pred85_ExpVsMSL$fit

agFishSiteLongTBL$seBoot_ExpVsMSL[agFishSiteLongTBL$SSP == "SSP1 - 2.6"] <- pred26_ExpVsMSL$se.boot
agFishSiteLongTBL$seBoot_ExpVsMSL[agFishSiteLongTBL$SSP == "SSP5 - 8.5"] <- pred85_ExpVsMSL$se.boot


####  2) Main figures ####
### Set shapes and colors for countries
countryShapes <- c(21, 22, 23, 24, 25)
names(countryShapes) <- c("indonesia", "papua new guinea", "madagascar", "philippines", "tanzania")

countryColors <- c("#00264D", "#309292", "#86FAF2", "#93A42A", "#ECBE13")
names(countryColors) <- names(countryShapes)

countryLabels = c("Indonesia", "Papua New Guinea", "Madagascar", "Philippines", "Tanzania")
names(countryLabels) <- names(countryShapes)

# Change the order of countries
agFishSiteTBL <- agFishSiteTBL %>% 
  mutate(Country = fct_relevel(Country, c("tanzania", "madagascar", "indonesia", "philippines", "papua new guinea")))

####    a) Figure 1: Exposure vs. sensitivity with model agreement for agriculture and fisheries with map of study sites ####
#Set limits for plots
sensMax <- 0.32
sensBreaks <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3)
expLims <- c(-15,45)
expBreaks <- c(-10, 0, 10, 20, 30, 40)
modelAgreeLims <- c(45, 100)

####        i) Figure 1a - Plot of agriculture exposure (SSP5 - 8.5) vs. sensitivity ####
fig1a_AgExpVsSensPlot <- ggplot(data = agFishSiteTBL,
                                aes(x = S_agr,
                                    y = E_agr_Ssp585_Unweighted,
                                    shape = Country,
                                    fill = E_agr_Ssp585_PercModelAgree,
                                    ymin = E_agr_Ssp585_Unweighted25perc,
                                    ymax = E_agr_Ssp585_Unweighted75perc)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_segment(aes(x = 0.2, y = 0.5, xend = 0.2, yend = 10), arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_text(aes(x = 0.21, y = 5), label = "Losses", hjust = 0) +
  geom_segment(aes(x = 0.2, y = -0.5, xend = 0.2, yend = -10), arrow = arrow(length = unit(0.02, "npc"), type = "closed")) +
  geom_text(aes(x = 0.21, y = -5), label = "Gains", hjust = 0) +
  geom_errorbar() +
  geom_point(size = 3) +
  scale_x_continuous(name = "Agriculture sensitivity",
                     limits = c(0, sensMax),
                     breaks = sensBreaks) +
  scale_y_continuous(name = "Agriculture exposure\n(% loss by mid-century)",
                     limits = expLims,
                     breaks = expBreaks) +
  scale_shape_manual(values = countryShapes,
                     labels = countryLabels,
                     name = NULL) +
  scale_fill_distiller(palette = "YlGnBu", name = "% model agreement", direction = 1, limits = modelAgreeLims) +
  guides(shape = "none",
         fill = guide_colorbar(title.position = "top")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.75,
        legend.direction = "horizontal",
        legend.position = c(0.7, 0.7),
        legend.title = element_text(hjust = 1))

####        ii) Figure 1b - Plot of Fisheries exposure (SSP5 - 8.5) vs. sensitivity ####
fig1b_FishExpVsSensPlot <- ggplot(data = agFishSiteTBL,
                                  aes(x = S_fsh,
                                      y = E_fish_Ssp585,
                                      shape = Country,
                                      fill = E_fish_Ssp585_PercModelAgree,
                                      ymin = E_fish_Ssp585_25perc,
                                      ymax = E_fish_Ssp585_75perc)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbar() +
  geom_point(size = 3) +
  scale_x_continuous(name = "Fisheries sensitivity",
                     limits = c(0, sensMax),
                     breaks = sensBreaks) +
  scale_y_continuous(name = "Fisheries exposure\n(% loss by mid-century)",
                     limits = expLims,
                     breaks = expBreaks) +
  scale_shape_manual(values = countryShapes,
                     name = NULL) +
  scale_fill_distiller(palette = "YlGnBu", name = "% model agreement", direction = 1, limits = modelAgreeLims) +
  guides(shape = "none",
         fill = "none") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 0.75)

####        iii) Figure 1c - Map of the sites ####
###Set up map
mapWorld <- map_data('world2')

fig1c_StudySiteMap <- ggplot() + 
  geom_polygon(data = mapWorld, aes(x = long, y = lat, group = group), fill = "lightgrey", color = "lightgrey") +
  coord_fixed(xlim = c(35, 155), ylim = c(-25, 20)) +
  geom_point(data = agFishSiteTBL, 
             aes(x = lonDD, y = latDD, shape = Country, fill = E_agrfsh_Ssp585_AvgPercModelAgreement),
             size = 3) +
  scale_shape_manual(values = countryShapes,
                     labels = countryLabels,
                     name = NULL) +
  guides(shape = guide_legend(override.aes = list(fill = "grey", size = 3))) +
  scale_x_continuous("", expand = c(0,0)) +
  scale_y_continuous("", expand = c(0,0)) + 
  scale_fill_distiller(palette = "YlGnBu", name = "% model agreement", direction = 1, limits = modelAgreeLims, guide = NULL) +
  annotate("text",
           x = 154,
           y = 3,
           label = "Papua New Guinea\n(n=10)",
           hjust = 1,
           #color = countryColors["papua new guinea"],
           fontface = 2,
           size = 3) + 
  annotate("text",
           x = 136,
           y = 15,
           label = "Philippines\n(n=25)",
           hjust = 1,
           #color = countryColors["philippines"],
           fontface = 2,
           size = 3) + 
  annotate("text",
           x = 110,
           y = 0,
           label = "Indonesia\n(n=25)",
           hjust = 0.5,
           #color = countryColors["indonesia"],
           fontface = 2,
           size = 3) + 
  annotate("text",
           x = 42,
           y = -5,
           label = "Tanzania\n(n=6)",
           hjust = 0,
           #color = countryColors["tanzania"],
           fontface = 2,
           size = 3) + 
  annotate("text",
           x = 52,
           y = -15,
           label = "Madagascar\n(n=6)",
           hjust = 0,
           #color = countryColors["madagascar"],
           fontface = 2,
           size = 3) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "top",
        legend.key = element_rect(fill = NA),
        aspect.ratio = 0.4)

####        iv) Assemble figure 1 ####
fig1_AgFishExposureVsSensitivityWithMap <- ggpubr::ggarrange(
  #First row with plots of exposure vs sensitivity
  egg::ggarrange(fig1a_AgExpVsSensPlot,
                 fig1b_FishExpVsSensPlot,
                 ncol = 2, labels = c("A", "B"),
                 label.args = list(gp = grid::gpar(font = 2))),
  #Second row with map
  egg::ggarrange(fig1c_StudySiteMap,
                 ncol = 1, labels = "C",
                 label.args = list(gp = grid::gpar(font = 2))),
  nrow = 2, heights = c(1,1.4)
)

ggsave(plot = fig1_AgFishExposureVsSensitivityWithMap, 
       filename = paste0(plotDir, "AgFish_Fig1_AgFishExposureVsSensitivityWithStudySites.tiff"),
       width = 8,
       height = 6)

####    b) Figure 2: Agriculture vs. fisheries exposure and comparison with random locations ####
#Combine the exposure data for the study sites and random sites
randomSiteExposureTBL <- agFishRandTBL %>% 
  dplyr::select_at(vars(matches("E_"))) %>% 
  dplyr::mutate(Type = "Random")

studySiteExposureTBL <- agFishSiteTBL %>% 
  dplyr::select(colnames(agFishSiteTBL)[colnames(agFishSiteTBL) %in% colnames(randomSiteExposureTBL)]) %>% 
  dplyr::select_at(vars(matches("E_"))) %>% 
  dplyr::mutate(Type = "Study")

allSitesExposureTBL <- bind_rows(randomSiteExposureTBL, studySiteExposureTBL)

studyRandomColors <- c("black", "grey")
names(studyRandomColors) <- c("Study", "Random")

studyRandomSizes <- c(2, 1)
names(studyRandomSizes) <- names(studyRandomColors)

#Create data for average lines (from the mixed effects models including country as a random effect)
fig2_AverageLinesTBL <- tibble(Type = c("Random", "Study"),
                               Avg_E_fish_Ssp585 = c(fishLossesRandSsp585MeanEst, fishLossesSsp585MeanEst),
                               SE_E_fish_Ssp585 = c(fishLossesRandSsp585SE, fishLossesSsp585SE),
                               Avg_E_agr_Ssp585 = c(agLossesRandSsp585MeanEst, agLossesSsp585MeanEst),
                               SE_E_agr_Ssp585 = c(agLossesRandSsp585SE, agLossesSsp585SE))

####      i) Figure 2 - main plot - Agriculture vs. fisheries exposure for study and random sites (RCP 8.5) ####
fig2_agExpLims <- c(-15, 4)
fig2_agExpBreaks <- c(-14, -12, -10, -8, -6, -4, -2, 0, 2, 4)
fig2_fishExpLims <- c(-5, 45) 

fig2main_AgVsFishExpWithRandomPlot <- ggplot() + 
  geom_point(data = allSitesExposureTBL,
             aes(x = E_fish_Ssp585,
                 y = E_agr_Ssp585_Unweighted,
                 color = Type),
             size = 2) +
  geom_vline(data = fig2_AverageLinesTBL,
             aes(xintercept = Avg_E_fish_Ssp585,
                 color = Type), linetype = "dashed") +
  geom_hline(data = fig2_AverageLinesTBL,
             aes(yintercept = Avg_E_agr_Ssp585,
                 color = Type), linetype = "dashed") +
  scale_color_manual(values = studyRandomColors,
                     name = NULL) +
  guides(color = "none", size = "none") +
  scale_x_continuous(name = "Fisheries exposure\n(% loss by mid-century)",
                     limits = fig2_fishExpLims) +
  scale_y_continuous(name = "Agriculture exposure\n(% loss by mid-century)",
                     limits = fig2_agExpLims, breaks = fig2_agExpBreaks) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

####      ii) Figure 2 - horizontal histogram - fisheries exposure density plot ####
fig2horizhist_FishExpCountDensityPlot <- ggplot() + 
  geom_density(data = allSitesExposureTBL %>% filter(Type == "Random"),
               aes(E_fish_Ssp585, after_stat(count), fill = Type)) +
  geom_density(data = allSitesExposureTBL %>% filter(Type == "Study"),
               aes(E_fish_Ssp585, after_stat(count)*20, fill = Type)) +
  scale_fill_manual(values = studyRandomColors) +
  scale_x_continuous(limits = fig2_fishExpLims) +
  guides(fill = "none") +
  theme_void() +
  theme(aspect.ratio = 1/4)

####      iii) Figure 2 - vertical histogram - agriculture exposure density ####
fig2verthist_AgExpCountDensityPlot <- ggplot() + 
  geom_density(data = allSitesExposureTBL %>% filter(Type == "Random"),
               aes(E_agr_Ssp585_Unweighted, after_stat(count), fill = Type)) +
  geom_density(data = allSitesExposureTBL %>% filter(Type == "Study"),
               aes(E_agr_Ssp585_Unweighted, after_stat(count)*20, fill = Type)) +
  scale_fill_manual(values = studyRandomColors) +
  scale_x_continuous(limits = fig2_agExpLims) +
  guides(fill = "none") +
  theme_void() +
  coord_flip() +
  theme(aspect.ratio = 4/1)

blankPlot <- ggplot()+geom_blank()+ theme_void() #Add blank plot to make space

######        iv) Assemble figure 2 ####
fig2_AgVsFishExpWithMarginalDensPlot <- egg::ggarrange(fig2horizhist_FishExpCountDensityPlot,
                                                          blankPlot,
                                                          fig2main_AgVsFishExpWithRandomPlot,
                                                          fig2verthist_AgExpCountDensityPlot,
                                                          ncol = 2,
                                                          widths = c(4, 1),
                                                          heights = c(1, 4))

ggsave(plot = fig2_AgVsFishExpWithMarginalDensPlot, 
       filename = paste0(plotDir, "AgFish_Fig2_AgVsFishExposureWithMarginalDensities.tiff"),
       width = 8,
       height = 8)

####    c) Figure 3: Agriculture vs. fisheries exposure with sector dependency and across SSP's ####
# Calculate a "relative sector dependency" (i.e. difference between fisheries and agriculture sensitivity)
agFishSiteTBL <- agFishSiteTBL %>%
  dplyr::mutate(S_index = S_fsh - S_agr)

agExpLim <- c(-15, 13)
fishExpLim <- c(-5, 45)

####      i) Figure 3a - plot of agriculture vs. fisheries exposure with bubble size and color related to sensitivity ####
fig3a_AgVsFishExpRelCumDepPlot <- ggplot(data = agFishSiteTBL %>% arrange(desc(S_agrfsh)),
                                         aes(x = E_fish_Ssp585,
                                             y = E_agr_Ssp585_Unweighted)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_linerange(aes(xmin = E_fish_Ssp585_25perc, xmax = E_fish_Ssp585_75perc, color = S_index)) +
  geom_linerange(aes(ymin = E_agr_Ssp585_Unweighted25perc, ymax = E_agr_Ssp585_Unweighted75perc, color = S_index)) +
  geom_point(aes(fill = S_index, size = S_agrfsh), color = "black", pch = 21) +
  scale_size_continuous(breaks = c(0.1, 0.25, 0.4), range = c(0.01,10)) +
  scale_fill_gradient2(
    low = "#eec71a",
    mid = "lightgrey",
    high = "#007ecc",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "fill",
    limits = c(-0.2, 0.4),
    breaks = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4)
  ) +
  scale_colour_gradient2(
    low = "#eec71a",
    mid = "lightgrey",
    high = "#007ecc",
    midpoint = 0,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour",
    limits = c(-0.2, 0.4),
    breaks = c(-0.2, -0.1, 0, 0.1, 0.2, 0.3, 0.4)
  ) +
  guides(color = "none",
         size = guide_legend(title = "Cumulative sector\ndependency", title.hjust = 0.5, override.aes = list(fill = "grey")),
         fill = guide_colorbar(title = "Relative sector\ndependency\n(Fisheries-Agriculture)", title.hjust = 0.5)
  ) +
  scale_x_continuous(name = "Fisheries exposure\n(% loss by mid-century)",
                     limits = fishExpLim) +
  scale_y_continuous(name = "Agriculture exposure\n(% loss by mid-century)",
                     limits = agExpLim) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

####      ii) Figure 3b - Effects of climate mitigation (comparing RCP 8.5 and 2.6) ####
fig3b_ClimateMitigationPlot <- ggplot(data = agFishSiteTBL) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_segment(aes(x = E_fish_Ssp585, y = E_agr_Ssp585_Unweighted,
                   xend = E_fish_Ssp126, yend = E_agr_Ssp126_Unweighted),
               color = "#301818") +
  geom_point(aes(x = E_fish_Ssp585,
                 y = E_agr_Ssp585_Unweighted),
             shape = 21,
             size = 3,
             fill = "#301818") +
  geom_point(aes(x = E_fish_Ssp126,
                 y = E_agr_Ssp126_Unweighted),
             shape = 21,
             size = 3,
             fill = "gold3") +
  geom_text(aes(x = 35, y = 9.5), vjust = 0, label = "SSP5 - 8.5", color = "#301818") +
  geom_text(aes(x = 35, y = 9), vjust = 1, label = "SSP1 - 2.6", color = "gold3") +
  scale_x_continuous(name = "Fisheries exposure\n(% loss by mid-century)",
                     limits = fishExpLim) +
  scale_y_continuous(name = "Agriculture exposure\n(% loss by mid-century)",
                     limits = agExpLim) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

######      iii) Assemble figure 3 ####
figure3_AgVsFishExpByDepRcpPlot <- egg::ggarrange(fig3a_AgVsFishExpRelCumDepPlot,
                                                  fig3b_ClimateMitigationPlot,
                                                  ncol = 2, labels = c("A", "B"),
                                                  label.args = list(gp = grid::gpar(font = 2)))

ggsave(plot = figure3_AgVsFishExpByDepRcpPlot, 
       filename = paste0(plotDir, "AgFish_Fig3_AgVsFishExposureByDependencyRCP.tiff"),
       width = 8,
       height = 6)

####    d) Figure 4: Potential impacts vs. MSL across SSP's ####
##Create the labels for the plots
agFishSiteLongTBL <- agFishSiteLongTBL %>% 
  mutate(pvalLabel = ifelse(test = Pvalue < 0.001, yes = "p < 0.001", no = paste0("p = ", round(Pvalue, digits = 3))),
         r2label = paste0("R<sup>2</sup> = ", round(x = R2m, digits = 3), " (m); ", round(x = R2c, digits = 3), " (c)"),
         fullLabel = paste0(pvalLabel, "<br>", r2label))

labelsTBL <- agFishSiteLongTBL %>% 
  dplyr::select(SSP, pvalLabel, r2label, fullLabel) %>% 
  distinct() 

fig4_PotImpVsMslColByCountryFacetByRcpPlot <- ggplot(data = agFishSiteLongTBL) +
  facet_rep_wrap(~ SSP) +
  geom_line(aes(x = MSL_rescale, y = fit)) +
  geom_ribbon(aes(x = MSL_rescale, ymin = fit-seBoot, ymax = fit+seBoot), alpha = 0.25) +
  geom_point(aes(x = MSL_rescale, y = PotImp, fill = Country, shape = Country), size = 2) +
  geom_richtext(data = labelsTBL, aes(x = 0,  y = 0.05, label = fullLabel), hjust = 0, fill = NA, label.color = NA) +
  scale_fill_manual(values = countryColors,
                    labels = countryLabels,
                    name = NULL) +
  scale_shape_manual(values = countryShapes,
                     labels = countryLabels,
                     name = NULL) +
  coord_cartesian(xlim = c(0,1), ylim = c(0,1)) +
  scale_x_continuous(name = "Material style of life",
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_y_continuous(name = "Potential impact",
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0.5, face = "bold"),
        legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

ggsave(plot = fig4_PotImpVsMslColByCountryFacetByRcpPlot, 
       filename = paste0(plotDir, "AgFish_Fig4_PotImpVsMslByCountryRCP.tiff"),
       width = 8,
       height = 6)

####    e) Figure 5: Sensitivity and material style of life over time in two PNG locations ####
####      i) Figure 5a - plot of agriculture vs. fisheries sensitivity over time in two PNG sites ####
fig5a_AgVsFishSensByTimePNG <- ggplot() +
  geom_abline(slope = 1, color = "grey", linetype = "dashed", size = 1) +
  geom_path(data = pngTimeTBL %>% filter(Community == "AH") %>% arrange(Year),
            aes(x = S_fsh, y = S_agr),
            color = "#0072B2") +
  geom_path(data = pngTimeTBL %>% filter(Community == "MU") %>% arrange(Year),
            aes(x = S_fsh, y = S_agr),
            color = "#D55E00") +
  geom_point(data = pngTimeTBL %>% filter(Community == "AH"),
             aes(x = S_fsh, y = S_agr),
             shape = 21, size = 2, fill = "#56B4E9", color = "#0072B2") +
  geom_point(data = pngTimeTBL %>% filter(Community == "MU"),
             aes(x = S_fsh, y = S_agr),
             shape = 21, size = 2, fill = "#E69F00", color = "#D55E00") +
  geom_text(data = pngTimeTBL %>% filter(Community == "AH" & Year == 2001),
            aes(x = S_fsh, y = S_agr, label = Year),
            color = "#56B4E9", nudge_x = -0.002, nudge_y = 0.01) +
  geom_text(data = pngTimeTBL %>% filter(Community == "AH" & Year == 2009),
            aes(x = S_fsh, y = S_agr, label = Year),
            color = "#56B4E9", nudge_x = 0.01) +
  geom_text(data = pngTimeTBL %>% filter(Community == "AH" & Year == 2012),
            aes(x = S_fsh, y = S_agr, label = Year),
            color = "#56B4E9", nudge_y = 0.005, nudge_x = 0.01) +
  geom_text(data = pngTimeTBL %>% filter(Community == "AH" & Year == 2016),
            aes(x = S_fsh, y = S_agr, label = Year),
            color = "#56B4E9", nudge_y = 0.006, nudge_x = -0.003) +
  geom_text(data = pngTimeTBL %>% filter(Community == "AH" & Year == 2018),
            aes(x = S_fsh, y = S_agr, label = Year),
            color = "#56B4E9", nudge_y = 0, nudge_x = -0.01) +
  geom_text(data = pngTimeTBL %>% filter(Community == "MU" & Year == 2001),
            aes(x = S_fsh, y = S_agr, label = Year),
            color = "#E69F00", nudge_x = -0.01, nudge_y = -0.005) +
  geom_text(data = pngTimeTBL %>% filter(Community == "MU" & Year == 2009),
            aes(x = S_fsh, y = S_agr, label = Year),
            color = "#E69F00", nudge_x = 0.01) +
  geom_text(data = pngTimeTBL %>% filter(Community == "MU" & Year == 2012),
            aes(x = S_fsh, y = S_agr, label = Year),
            color = "#E69F00", nudge_y = 0, nudge_x = 0.012) +
  geom_text(data = pngTimeTBL %>% filter(Community == "MU" & Year == 2016),
            aes(x = S_fsh, y = S_agr, label = Year),
            color = "#E69F00", nudge_y = -0.005, nudge_x = 0) +
  geom_point(data = pngTimeTBL, aes(x = 0.09, y = 0.02), color = "#56B4E9", size = 10) +
  geom_text(data = pngTimeTBL, aes(x = 0.09, y = 0.02), color = "#0072B2", label = "Ahus") +
  geom_point(data = pngTimeTBL, aes(x = 0.04, y = 0.11), color = "#E69F00", size = 10) +
  geom_text(data = pngTimeTBL, aes(x = 0.04, y = 0.11), color = "#D55E00", label = "Muluk") +
  scale_x_continuous(name = "Fisheries sensitivity",
                     limits = c(0, 0.14),
                     breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14)) +
  scale_y_continuous(name = "Agriculture sensitivity",
                     limits = c(0, 0.14),
                     breaks = c(0, 0.02, 0.04, 0.06, 0.08, 0.1, 0.12, 0.14)) +
  geom_richtext(data = pngTimeTBL, aes(x = 0.09, y = 0.12),
                label = "More sensitive to<br>agricultural impacts", angle = 45, color = "#D55E00",
                fill = NA, label.color = NA) +
  geom_segment(data = pngTimeTBL, aes(x = 0.105, y = 0.105, xend = 0.095, yend = 0.115),
               arrow = arrow(length = unit(0.2, "cm")), color = "#D55E00") +
  geom_segment(data = pngTimeTBL, aes(x = 0.105, y = 0.105, xend = 0.115, yend = 0.095),
               arrow = arrow(length = unit(0.2, "cm")), color = "#0072B2") +
  geom_richtext(data = pngTimeTBL, aes(x = 0.12, y = 0.09),
                label = "More sensitive to<br>fisheries impacts", angle = 45, color = "#0072B2",
                fill = NA, label.color = NA) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

####      ii) Figure 5b - plot of MSL PCA axes over time for two PNG sites ####
fig5b_MSLPC1vs2_ByTimePNG <- ggplot() +
  geom_vline(xintercept = 0, color = "grey", linetype = "dashed", size = 1) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed", size = 1) +
  geom_segment(data = mslPcaResTBL %>% filter(!is.na(Label)), aes(x = 0, y = 0, xend = RC1, yend = RC2), color = "grey") +
  geom_text(data = mslPcaResTBL %>% filter(!is.na(Label)), aes(x = RC1, y = RC2, label = Label), color = "black", size = 2,
            position = position_jitter(width = 0.01, height = 0.01)) +
  geom_path(data = pngTimeTBL %>% filter(Community == "AH") %>% arrange(Year),
            aes(x = MSL1, y = MSL2),
            color = "#0072B2") +
  geom_path(data = pngTimeTBL %>% filter(Community == "MU") %>% arrange(Year),
            aes(x = MSL1, y = MSL2),
            color = "#D55E00") +
  geom_point(data = pngTimeTBL %>% filter(Community == "AH"),
             aes(x = MSL1, y = MSL2),
             shape = 21, size = 2, fill = "#56B4E9", color = "#0072B2") +
  geom_point(data = pngTimeTBL %>% filter(Community == "MU"),
             aes(x = MSL1, y = MSL2),
             shape = 21, size = 2, fill = "#E69F00", color = "#D55E00") +
  geom_text(data = pngTimeTBL %>% filter(Community == "AH" & Year == 2001),
            aes(x = MSL1, y = MSL2, label = Year),
            color = "#56B4E9", nudge_x = 0.16, nudge_y = 0.1) +
  geom_text(data = pngTimeTBL %>% filter(Community == "AH" & Year == 2009),
            aes(x = MSL1, y = MSL2, label = Year),
            color = "#56B4E9", nudge_y = 0.07) +
  geom_text(data = pngTimeTBL %>% filter(Community == "AH" & Year == 2012),
            aes(x = MSL1, y = MSL2, label = Year),
            color = "#56B4E9", nudge_y = 0, nudge_x = 0.2) +
  geom_text(data = pngTimeTBL %>% filter(Community == "AH" & Year == 2016),
            aes(x = MSL1, y = MSL2, label = Year),
            color = "#56B4E9", nudge_y = 0, nudge_x = 0.2) +
  geom_text(data = pngTimeTBL %>% filter(Community == "AH" & Year == 2018),
            aes(x = MSL1, y = MSL2, label = Year),
            color = "#56B4E9", nudge_y = 0, nudge_x = -0.2) +
  geom_text(data = pngTimeTBL %>% filter(Community == "MU" & Year == 2001),
            aes(x = MSL1, y = MSL2, label = Year),
            color = "#E69F00", nudge_x = -0.05, nudge_y = -0.06) +
  geom_text(data = pngTimeTBL %>% filter(Community == "MU" & Year == 2009),
            aes(x = MSL1, y = MSL2, label = Year),
            color = "#E69F00", nudge_x = 0.22) +
  geom_text(data = pngTimeTBL %>% filter(Community == "MU" & Year == 2012),
            aes(x = MSL1, y = MSL2, label = Year),
            color = "#E69F00", nudge_y = 0.05, nudge_x = 0.2) +
  geom_text(data = pngTimeTBL %>% filter(Community == "MU" & Year == 2016),
            aes(x = MSL1, y = MSL2, label = Year),
            color = "#E69F00", nudge_y = -0.06, nudge_x = 0) +
  geom_point(data = pngTimeTBL, aes(x = -1.4, y = -0.9), color = "#56B4E9", size = 10) +
  geom_text(data = pngTimeTBL, aes(x = -1.4, y = -0.9), color = "#0072B2", label = "Ahus") +
  geom_point(data = pngTimeTBL, aes(x = -1, y = 0.9), color = "#E69F00", size = 10) +
  geom_text(data = pngTimeTBL, aes(x = -1, y = 0.9), color = "#D55E00", label = "Muluk") +
  scale_x_continuous(name = paste0("Material Style of Life PC1 (", round(mslPcaResTBL$RC1[mslPcaResTBL$Item == "Proportion Var"]*100, digits = 1), "%)"),
                     limits = c(-1.8, 1),
                     breaks = c(-1.5, -1, -0.5, 0, 0.5, 1)) +
  scale_y_continuous(name = paste0("Material Style of Life PC2 (", round(mslPcaResTBL$RC2[mslPcaResTBL$Item == "Proportion Var"]*100, digits = 1), "%)"),
                     limits = c(-1, 1),
                     breaks = c(-1, -0.5, 0, 0.5, 1)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)
  
######        iii) Assemble figure 5 ####
figure5_SensMslByTimePngPlot <- egg::ggarrange(fig5a_AgVsFishSensByTimePNG,
                                                  fig5b_MSLPC1vs2_ByTimePNG,
                                                  ncol = 2, labels = c("A", "B"),
                                                  label.args = list(gp = grid::gpar(font = 2)))

ggsave(plot = figure5_SensMslByTimePngPlot, 
       filename = paste0(plotDir, "AgFish_Fig5_SensMslOverTimePng.tiff"),
       width = 8,
       height = 6)


####  3) Supplementary figures ####
####    a) Figure S1: Projected agricultural changes by crop ####
#Create a single dataset with the random and study sites combined for agriculture exposure
agExpRandLongTBL <- agFishRandTBL %>% 
  dplyr::select(pointid, Country, Lat_DD, Long_DD, E_maize_Ssp126, E_rice_Ssp126, E_cassava_Ssp126,
                E_agr_Ssp126_Weighted, E_agr_Ssp126_Unweighted, E_maize_Ssp585, E_rice_Ssp585, E_cassava_Ssp585,
                E_agr_Ssp585_Weighted, E_agr_Ssp585_Unweighted) %>% 
  gather(key = "Metric_SSP", value = "Exposure", "E_maize_Ssp126":"E_agr_Ssp585_Unweighted") %>% 
  mutate(Metric_SSP = dplyr::recode(Metric_SSP,
                             E_maize_Ssp126 = "Maize_Ssp126",
                             E_rice_Ssp126 = "Rice_Ssp126",
                             E_cassava_Ssp126 = "Cassava_Ssp126",
                             E_agr_Ssp126_Unweighted = "Total_Ssp126",
                             E_agr_Ssp126_Weighted = "Weighted_Ssp126",
                             E_maize_Ssp585 = "Maize_Ssp585",
                             E_rice_Ssp585 = "Rice_Ssp585",
                             E_cassava_Ssp585 = "Cassava_Ssp585",
                             E_agr_Ssp585_Unweighted = "Total_Ssp585",
                             E_agr_Ssp585_Weighted = "Weighted_Ssp585"),
         Crop = sapply(strsplit(Metric_SSP,"_"), `[`, 1),
         SSP = sapply(strsplit(Metric_SSP,"_"), `[`, 2),
         Data = "Random sites") 
         
agExpLongTBL <- agFishSiteTBL %>% 
  dplyr::select(SiteNum, Country, latDD, lonDD, E_maize_Ssp126, E_rice_Ssp126, E_cassava_Ssp126,
                E_agr_Ssp126_Weighted, E_agr_Ssp126_Unweighted, E_maize_Ssp585, E_rice_Ssp585, E_cassava_Ssp585,
                E_agr_Ssp585_Weighted, E_agr_Ssp585_Unweighted) %>% 
  rename(pointid = SiteNum, Lat_DD = latDD, Long_DD = lonDD) %>% 
  gather(key = "Metric_SSP", value = "Exposure", "E_maize_Ssp126":"E_agr_Ssp585_Unweighted") %>% 
  mutate(Metric_SSP = dplyr::recode(Metric_SSP,
                                    E_maize_Ssp126 = "Maize_Ssp126",
                                    E_rice_Ssp126 = "Rice_Ssp126",
                                    E_cassava_Ssp126 = "Cassava_Ssp126",
                                    E_agr_Ssp126_Unweighted = "Total_Ssp126",
                                    E_agr_Ssp126_Weighted = "Weighted_Ssp126",
                                    E_maize_Ssp585 = "Maize_Ssp585",
                                    E_rice_Ssp585 = "Rice_Ssp585",
                                    E_cassava_Ssp585 = "Cassava_Ssp585",
                                    E_agr_Ssp585_Unweighted = "Total_Ssp585",
                                    E_agr_Ssp585_Weighted = "Weighted_Ssp585"),
         Crop = sapply(strsplit(Metric_SSP,"_"), `[`, 1),
         SSP = sapply(strsplit(Metric_SSP,"_"), `[`, 2),
         Data = "Study sites")

agExpAllLongTBL <- bind_rows(agExpLongTBL, agExpRandLongTBL) %>% 
  mutate(SSP = dplyr::recode(SSP,
                             Ssp126 = "SSP1 - 2.6",
                             Ssp585 = "SSP5 - 8.5"),
         Data = factor(Data, levels = c("Study sites", "Random sites")))

#Set colors for the figure
cropColors <- c("#C4961A", "#FFDB6D", "#F4EDCA", "#52854C", "#C3D7A4")
names(cropColors) <- c("Cassava", "Maize", "Rice", "Total", "Weighted")
  
figS1_ExpByCropSspRandStudyFacetsBoxplots <- ggplot(data = agExpAllLongTBL,
                                                    aes(x = Crop,
                                                        y = Exposure)) +
  facet_grid(vars(Data), vars(SSP)) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_boxplot(aes(fill = Crop)) +
  scale_fill_manual(values = cropColors) +
  scale_y_continuous(limits = c(-60, 15),
                     breaks = c(-50, -40, -30, -30, -20, -10, 0, 10)) +
  ylab(label = "Production loss (%)") +
  xlab(label = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none",
        strip.background = element_blank())
  
ggsave(plot = figS1_ExpByCropSspRandStudyFacetsBoxplots, 
       filename = paste0(plotDir, "AgFish_FigS1_AgExpByCropSspRandomStudyFacetsBoxplots.tiff"),
       width = 8,
       height = 8)

####    b) Figure S2: Agriculture vs. fisheries exposure with random sites compared between SSP's ####
#Create data for average lines (from the mixed effects models including country as a random effect)
figS2_AverageLinesTBL <- tibble(Type = c("Random", "Study"),
                               Avg_E_fish_Ssp126 = c(fishLossesRandSsp126MeanEst, fishLossesSsp126MeanEst),
                               SE_E_fish_Ssp126 = c(fishLossesRandSsp126SE, fishLossesSsp126SE),
                               Avg_E_agr_Ssp126 = c(agLossesRandSsp126MeanEst, agLossesSsp126MeanEst),
                               SE_E_agr_Ssp126 = c(agLossesRandSsp126SE, agLossesSsp126SE))


####      i) Figure S2 - main plot - Agriculture vs. fisheries exposure for study and random sites (RCP 8.5) ####
figS2_main_AgVsFishExpWithRandom126Plot <- ggplot() + 
  geom_point(data = allSitesExposureTBL,
             aes(x = E_fish_Ssp126,
                 y = E_agr_Ssp126_Unweighted,
                 color = Type),
             size = 2) +
  geom_vline(data = figS2_AverageLinesTBL,
             aes(xintercept = Avg_E_fish_Ssp126,
                 color = Type), linetype = "dashed") +
  geom_hline(data = figS2_AverageLinesTBL,
             aes(yintercept = Avg_E_agr_Ssp126,
                 color = Type), linetype = "dashed") +
  scale_color_manual(values = studyRandomColors,
                     name = NULL) +
  guides(color = "none", size = "none") +
  scale_x_continuous(name = "Fisheries exposure\n(% loss by mid-century)",
                     limits = fig2_fishExpLims) +
  scale_y_continuous(name = "Agriculture exposure\n(% loss by mid-century)",
                     limits = fig2_agExpLims, breaks = fig2_agExpBreaks) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

####      ii) Figure S2 - horizontal histogram - fisheries exposure density plot ####
figS2horizhist_FishExpCountDensityPlot <- ggplot() + 
  geom_density(data = allSitesExposureTBL %>% filter(Type == "Random"),
               aes(E_fish_Ssp126, after_stat(count), fill = Type)) +
  geom_density(data = allSitesExposureTBL %>% filter(Type == "Study"),
               aes(E_fish_Ssp126, after_stat(count)*20, fill = Type)) +
  scale_fill_manual(values = studyRandomColors) +
  scale_x_continuous(limits = fig2_fishExpLims) +
  guides(fill = "none") +
  theme_void() +
  theme(aspect.ratio = 1/4)

####      iii) Figure S2 - vertical histogram - agriculture exposure density ####
figS2verthist_AgExpCountDensityPlot <- ggplot() + 
  geom_density(data = allSitesExposureTBL %>% filter(Type == "Random"),
               aes(E_agr_Ssp126_Unweighted, after_stat(count), fill = Type)) +
  geom_density(data = allSitesExposureTBL %>% filter(Type == "Study"),
               aes(E_agr_Ssp126_Unweighted, after_stat(count)*20, fill = Type)) +
  scale_fill_manual(values = studyRandomColors) +
  scale_x_continuous(limits = fig2_agExpLims) +
  guides(fill = "none") +
  theme_void() +
  coord_flip() +
  theme(aspect.ratio = 4/1)

######        iv) Assemble figure S2 ####
figS2_AgVsFishExpWithMarginalDensPlot <- egg::ggarrange(figS2horizhist_FishExpCountDensityPlot,
                                                        blankPlot,
                                                        figS2_main_AgVsFishExpWithRandom126Plot,
                                                        figS2verthist_AgExpCountDensityPlot,
                                                        ncol = 2,
                                                        widths = c(4, 1),
                                                        heights = c(1, 4))

###         Now combine that with figure 2 (SSP5 - 8.5)
figS2_AgVsFishExpWithMargDensPlotsSsps126and585 <- ggarrange(figS2_AgVsFishExpWithMarginalDensPlot,
                                                             fig2_AgVsFishExpWithMarginalDensPlot,
                                                             labels = c("A.        SSP1 - 2.6", "B.        SSP5 - 8.5"),
                                                             ncol = 2, label.y = 0.9, label.x = 0.5, hjust = 1)

ggsave(plot = figS2_AgVsFishExpWithMargDensPlotsSsps126and585, 
       filename = paste0(plotDir, "AgFish_FigS2_AgVsFishExposureWithMarginalDensitiesBothSSP.tiff"),
       width = 8,
       height = 6)

####    c) Figure S3: Combined exposure vs combined sensitivity across SSPs ####
####      i) Fig S3a - Agriculture-fisheries exposure vs. sensitivity for RCP 2.6 and 8.5 as two facets ####
figS3a_AgFishExpVsSensByCountryPotImpSspFacetScatterplots <- ggplot(data = agFishSiteLongTBL) +
  facet_rep_wrap(~ SSP) +
  geom_point(aes(x = S_agrfsh_rescale, y = E_agrfsh_rescale, size = PotImp, fill = Country, shape = Country)) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                     name = "Agriculture and fisheries sensitivity") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                     name = "Agriculture and fisheries exposure") +
  scale_shape_manual(values = countryShapes,
                     labels = countryLabels,
                     name = NULL) +
  scale_fill_manual(values = countryColors,
                    labels = countryLabels,
                    name = NULL) +
  scale_size_continuous(name = "Magnitude of potential impact:",
                        limits = c(0,1),
                        breaks = c(0.25, 0.5, 0.75, 1)) +
  guides(fill = guide_legend(override.aes = list(size = 4)),
         size = guide_legend(override.aes = list(fill = "grey", shape = 21))) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

figS3_legend <- get_legend(figS3a_AgFishExpVsSensByCountryPotImpSspFacetScatterplots)

figS3a_AgFishExpVsSensByCountryPotImpSspFacetScatterplotsNoLegend <- ggplot(data = agFishSiteLongTBL) +
  facet_rep_wrap(~ SSP) +
  geom_point(aes(x = S_agrfsh_rescale, y = E_agrfsh_rescale, size = PotImp, fill = Country, shape = Country)) +
  scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                     limits = c(0,1),
                     name = "Agriculture and fisheries sensitivity") +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                     limits = c(0,1),
                     name = "Agriculture and fisheries exposure") +
  scale_shape_manual(values = countryShapes,
                     labels = countryLabels,
                     name = NULL) +
  scale_fill_manual(values = countryColors,
                    labels = countryLabels,
                    name = NULL) +
  scale_size_continuous(name = "Magnitude of potential impact:",
                        limits = c(0,1),
                        breaks = c(0.25, 0.5, 0.75, 1)) +
  guides(fill = FALSE,
         size = FALSE,
         shape = FALSE) +
  theme_bw() +
  theme(strip.background = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.box = "vertical",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

####      ii) Fig S3b - Plot of the difference between RCP 2.6 and 8.5 ####
#Split the long dataset into 2.6 and 8.5 then subtract to get the difference
agFishSiteLong126TBL <- agFishSiteLongTBL %>% 
  filter(SSP == "SSP1 - 2.6") %>% 
  dplyr::select(SiteNum, Country, S_agrfsh_rescale, E_agrfsh_rescale, PotImp) %>% 
  rename(S_agrfsh_rescale126 = S_agrfsh_rescale,
         E_agrfsh_rescale126 = E_agrfsh_rescale,
         PotImp126 = PotImp)

agFishSiteLong585TBL <- agFishSiteLongTBL %>% 
  filter(SSP == "SSP5 - 8.5") %>% 
  dplyr::select(SiteNum, Country, S_agrfsh_rescale, E_agrfsh_rescale, PotImp) %>% 
  rename(S_agrfsh_rescale585 = S_agrfsh_rescale,
         E_agrfsh_rescale585 = E_agrfsh_rescale,
         PotImp585 = PotImp)

diffAgFishSiteLongTBL <- left_join(agFishSiteLong126TBL, agFishSiteLong585TBL)  %>%
  mutate(E_agrfsh_rescaleDiff = E_agrfsh_rescale585 - E_agrfsh_rescale126,
         PotImpDiff = abs(PotImp585 - PotImp126))

figS3b_DiffAgFishExpVsSensByCountryPotImpScatterplot <- ggplot(data = diffAgFishSiteLongTBL) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
  geom_point(aes(x = S_agrfsh_rescale585, y = E_agrfsh_rescaleDiff, size = PotImpDiff, fill = Country, shape = Country)) +
  scale_x_continuous(name = "Agriculture and fisheries sensitivity",
                     limits = c(0,1),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_y_continuous(name = "\u0394 Agriculture and fisheries exposure",
                     limits = c(-0.2,1),
                     breaks = c(-0.2, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_shape_manual(values = countryShapes,
                     labels = countryLabels,
                     name = NULL) +
  scale_fill_manual(values = countryColors,
                    labels = countryLabels,
                    name = NULL) +
  scale_size_continuous(name = "Magnitude of potential impact:",
                        limits = c(0,1),
                        breaks = c(0.25, 0.5, 0.75, 1)) +
  guides(fill = FALSE,
         size = FALSE,
         shape = FALSE) +
  labs(subtitle = "SSP1 - 2.6 to SSP5 - 8.5") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        plot.subtitle = element_text(hjust = 0.5),
        aspect.ratio = 1)

######      iii) Assemble figure S3 ####
figS3_AgFishExpVsSensCompSsp126and585Plots <- ggarrange(figS3a_AgFishExpVsSensByCountryPotImpSspFacetScatterplotsNoLegend,
                                                        figS3_legend,
                                                        figS3b_DiffAgFishExpVsSensByCountryPotImpScatterplot,
                                                        nrow = 3, labels = c("A", "", "B"),
                                                        heights = c(1, 0.3, 1))

ggsave(plot = figS3_AgFishExpVsSensCompSsp126and585Plots, 
       filename = paste0(plotDir, "AgFish_FigS3_AgFishExposureVsSensCompSsp126and585.tiff"),
       width = 8,
       height = 10)



####    d) Figure S4: Combined sensitivity and combined exposure vs. MSL ####
####      i) Fig S4a - Rescaled agriculture-fisheries sensitivity vs. MSL ####
figS4a_SvsMSL_ByCountryWithStatsPlot <- ggplot(data = agFishSiteTBL) +
  geom_line(aes(x = MSL_rescale, y = fitSvsMSL)) +
  geom_ribbon(aes(x = MSL_rescale, ymin = fitSvsMSL-seBootSvsMSL, ymax = fitSvsMSL+seBootSvsMSL), alpha = 0.25) +
  geom_point(aes(x = MSL_rescale, y = S_agrfsh_rescale, fill = Country, shape = Country), size = 4) +
  geom_richtext(aes(x = 1,  y = 1),
                label = paste0("p = ", round(pvalueSvsMSL, digits = 3),
                               "<br>R<sup>2</sup> = ", round(r2_SvsMSL[1], digits = 3), " (m); ",
                               round(r2_SvsMSL[2], digits = 3), " (c)"),
                hjust = 1, vjust = 1, fill = NA, label.color = NA) +
  scale_shape_manual(values = countryShapes,
                     labels = countryLabels,
                     name = NULL) +
  scale_fill_manual(values = countryColors,
                    labels = countryLabels,
                    name = NULL) +
  scale_x_continuous(name = "Material style of life",
                     limits = c(0, 1),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_y_continuous(name = "Rescaled agriculture and fisheries sensitivity",
                     limits = c(0, 1),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_bw() +
  theme(legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

figS4_legend <- get_legend(figS4a_SvsMSL_ByCountryWithStatsPlot)

figS4a_SvsMSL_ByCountryWithStatsPlotNoLegend <- ggplot(data = agFishSiteTBL) +
  geom_line(aes(x = MSL_rescale, y = fitSvsMSL)) +
  geom_ribbon(aes(x = MSL_rescale, ymin = fitSvsMSL-seBootSvsMSL, ymax = fitSvsMSL+seBootSvsMSL), alpha = 0.25) +
  geom_point(aes(x = MSL_rescale, y = S_agrfsh_rescale, fill = Country, shape = Country), size = 4) +
  geom_richtext(aes(x = 1,  y = 1),
                label = paste0("p = ", round(pvalueSvsMSL, digits = 3),
                               "<br>R<sup>2</sup> = ", round(r2_SvsMSL[1], digits = 3), " (m); ",
                               round(r2_SvsMSL[2], digits = 3), " (c)"),
                hjust = 1, vjust = 1, fill = NA, label.color = NA) +
  scale_shape_manual(values = countryShapes,
                     labels = countryLabels,
                     name = NULL) +
  scale_fill_manual(values = countryColors,
                    labels = countryLabels,
                    name = NULL) +
  scale_x_continuous(name = "Material style of life",
                     limits = c(0, 1),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_y_continuous(name = "Rescaled agriculture and fisheries sensitivity",
                     limits = c(0, 1),
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

####      ii) Fig S4b - Agriculture-fisheries exposure vs. MSL ####
##Create labels for S4b plots
agFishSiteLongTBL <- agFishSiteLongTBL %>% 
  mutate(pvalLabel_ExpVsMSL = ifelse(test = Pvalue_ExpVsMSL < 0.001,
                                     yes = "p < 0.001",
                                     no = paste0("p = ", round(Pvalue_ExpVsMSL, digits = 3))),
         r2label_ExpVsMSL = paste0("R<sup>2</sup> = ", round(x = Rm2_ExpVsMSL, digits = 3), " (m); ", round(x = Rc2_ExpVsMSL, digits = 3), " (c)"),
         fullLabel_ExpVsMSL = paste0(pvalLabel_ExpVsMSL, "<br>", r2label_ExpVsMSL))

labelsExpVsMslTBL <- agFishSiteLongTBL %>% 
  dplyr::select(SSP, fullLabel_ExpVsMSL) %>% 
  distinct() 

figS4b_ExpvsMSL_ByCountrySspWithStatsNoLegendPlot <- ggplot(data = agFishSiteLongTBL) +
  facet_rep_wrap(~ SSP) +
  geom_line(aes(x = MSL_rescale, y = fit_ExpVsMSL)) +
  geom_ribbon(aes(x = MSL_rescale, ymin = fit_ExpVsMSL-seBoot_ExpVsMSL, ymax = fit_ExpVsMSL+seBoot_ExpVsMSL), alpha = 0.25) +
  geom_point(aes(x = MSL_rescale, y = E_agrfsh_rescale, fill = Country, shape = Country), size = 4) +
  geom_richtext(data = labelsExpVsMslTBL, aes(x = 1,  y = 0.85, label = fullLabel_ExpVsMSL), hjust = 1, vjust = 1, fill = NA, label.color = NA) +
  scale_shape_manual(values = countryShapes,
                     labels = countryLabels,
                     name = NULL) +
  scale_fill_manual(values = countryColors,
                    labels = countryLabels,
                    name = NULL) +
  scale_x_continuous(name = "Material style of life",
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  scale_y_continuous(name = "Agriculture and fisheries exposure",
                     breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  coord_cartesian(ylim = c(0,1), xlim = c(0,1)) +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text.x = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        aspect.ratio = 1)

######      iii) Assemble figure S4 ####
figS4_AgFishSensExpVsMSL_CompSsp126and585Plots <- ggarrange(figS4a_SvsMSL_ByCountryWithStatsPlotNoLegend,
                                                            figS4_legend,
                                                            figS4b_ExpvsMSL_ByCountrySspWithStatsNoLegendPlot,
                                                            nrow = 3, labels = c("A", "", "B"),
                                                            heights = c(1, 0.3, 1))
                                                                     

ggsave(plot = figS4_AgFishSensExpVsMSL_CompSsp126and585Plots, 
       filename = paste0(plotDir, "AgFish_FigS4_AgFishSensitivityAndExposureVsMSL_CompSsp126and585.tiff"),
       width = 8,
       height = 10)

####    e) Figure S5: Map of locations from which agriculture and fisheries data were extracted  ####
expTypeColors <- c("#62BD69", "#73A5C6")
names(expTypeColors) <- c("Agriculture", "Fisheries")

figS5_AgFishExposureLocMap <- ggplot() + 
  geom_polygon(data = mapWorld, aes(x = long, y = lat, group = group), fill = "lightgrey", color = "lightgrey") +
  coord_fixed(xlim = c(35, 155), ylim = c(-25, 20)) +
  geom_tile(data = agFishExpLocsTBL, aes(x = LonDD, y = LatDD, fill = Cell, width = Width, height = Height)) +
  guides(shape = guide_legend(override.aes = list(size = 3))) +
  scale_x_continuous("", expand = c(0,0)) +
  scale_y_continuous("", expand = c(0,0)) + 
  scale_fill_manual(values = expTypeColors, name = NULL) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.key = element_rect(fill = NA),
        aspect.ratio = 0.4)

ggsave(plot = figS5_AgFishExposureLocMap, 
       filename = paste0(plotDir, "AgFish_FigS5_AgFishExposureExtractionLocationMap.tiff"),
       width = 8,
       height = 6)

####    f) Figure S9: Finer scale maps showing locations of study sites in each country with average model agreement ####
#Get a fine scale coastline map to plot 
landZipLoc <- paste0(dataDir, "ne_10m_land.zip") 
download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_land.zip", landZipLoc)
utils::unzip(zipfile = landZipLoc, exdir = paste0(dataDir, "ne_10m_land"))
landSHP <- readOGR(dsn = paste0(dataDir, "ne_10m_land/ne_10m_land.shp"), stringsAsFactors = F)

smallIslandLoc <- paste0(dataDir, "ne_10m_minor_islands.zip") 
download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_minor_islands.zip", smallIslandLoc)
utils::unzip(zipfile = smallIslandLoc, exdir = paste0(dataDir, "ne_10m_minor_islands"))
smallIslandsSHP <- readOGR(dsn = paste0(dataDir, "ne_10m_minor_islands/ne_10m_minor_islands.shp"), stringsAsFactors = F)

coastlineLoc <- paste0(dataDir, "ne_10m_coastline.zip") 
download.file("https://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip", coastlineLoc)
utils::unzip(zipfile = coastlineLoc, exdir = paste0(dataDir, "ne_10m_coastline"))
coastlineSHP <- readOGR(dsn = paste0(dataDir, "ne_10m_coastline/ne_10m_coastline.shp"), stringsAsFactors = F)

####        i) Figure S9a - Tanzania and Madagascar ####
figS9a_TanzaniaMadagascarMap <- ggplot() + 
  geom_polygon(data = landSHP, aes(x = long, y = lat, group = group), fill = "lightgrey", color = "lightgrey") +
  geom_polygon(data = smallIslandsSHP, aes(x = long, y = lat, group = group), fill = "lightgrey", color = "lightgrey") +
  geom_path(data = coastlineSHP, aes(x = long, y = lat, group = group), color = "darkgrey") +
  coord_fixed(xlim = c(min(agFishSiteTBL$lonDD[agFishSiteTBL$Country %in% c("tanzania", "madagascar")]),
                       max(agFishSiteTBL$lonDD[agFishSiteTBL$Country %in% c("tanzania", "madagascar")])),
              ylim = c(min(agFishSiteTBL$latDD[agFishSiteTBL$Country %in% c("tanzania", "madagascar")]),
                       max(agFishSiteTBL$latDD[agFishSiteTBL$Country %in% c("tanzania", "madagascar")]))) +
  geom_point(data = agFishSiteTBL %>% filter(Country %in% c("tanzania", "madagascar")), 
             aes(x = lonDD, y = latDD, shape = Country, fill = E_agrfsh_Ssp585_AvgPercModelAgreement),
             size = 2) +
  scale_shape_manual(values = countryShapes,
                     labels = countryLabels,
                     name = NULL,
                     guide = NULL) +
  scale_x_continuous("") +
  scale_y_continuous("") + 
  scale_fill_distiller(palette = "YlGnBu", name = "Average % model agreement", direction = 1, limits = modelAgreeLims) +
  annotate("text",
           x = 40.5,
           y = -6.5,
           label = "Tanzania\n(n=6)",
           hjust = 0,
           vjust = 0.5) + 
  annotate("text",
           x = 47,
           y = -18,
           label = "Madagascar\n(n=6)",
           hjust = 0.5,
           vjust = 0.5) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

####        ii) Figure S9b - Indonesia, Philippines, and Papua New Guinea ####
figS9b_IndonesiaPhilippinesPngMap <- ggplot() + 
  geom_polygon(data = landSHP, aes(x = long, y = lat, group = group), fill = "lightgrey", color = "lightgrey") +
  geom_polygon(data = smallIslandsSHP, aes(x = long, y = lat, group = group), fill = "lightgrey", color = "lightgrey") +
  geom_path(data = coastlineSHP, aes(x = long, y = lat, group = group), color = "darkgrey") +
  coord_fixed(xlim = c(min(agFishSiteTBL$lonDD[agFishSiteTBL$Country %in% c("indonesia", "philippines", "papua new guinea")]),
                       max(agFishSiteTBL$lonDD[agFishSiteTBL$Country %in% c("indonesia", "philippines", "papua new guinea")])),
              ylim = c(min(agFishSiteTBL$latDD[agFishSiteTBL$Country %in% c("indonesia", "philippines", "papua new guinea")]),
                       max(agFishSiteTBL$latDD[agFishSiteTBL$Country %in% c("indonesia", "philippines", "papua new guinea")]))) +
  geom_point(data = agFishSiteTBL %>% filter(Country %in% c("indonesia", "philippines", "papua new guinea")), 
             aes(x = lonDD, y = latDD, shape = Country, fill = E_agrfsh_Ssp585_AvgPercModelAgreement),
             size = 2) +
  scale_shape_manual(values = countryShapes,
                     labels = countryLabels,
                     name = NULL,
                     guide = NULL) +
  scale_x_continuous("") +
  scale_y_continuous("") + 
  scale_fill_distiller(palette = "YlGnBu", name = "% model agreement", direction = 1, limits = modelAgreeLims) +
  annotate("text",
           x = 154,
           y = 2,
           label = "Papua New Guinea\n(n=10)",
           hjust = 1,
           vjust = 0.5) + 
  annotate("text",
           x = 127,
           y = 11,
           label = "Philippines\n(n=25)",
           hjust = 0,
           vjust = 0.5) + 
  annotate("text",
           x = 110,
           y = 0,
           label = "Indonesia\n(n=25)",
           hjust = 0.5,
           vjust = 0.5) + 
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom")

####        iii) Assemble figure S9 ####
figS9_FineScaleMapStudyLocsAvgAgreement <- ggpubr::ggarrange(
  #First column with Tanzania and Madagascar
  figS9a_TanzaniaMadagascarMap,
  #Second row with map
  figS9b_IndonesiaPhilippinesPngMap,
  ncol = 2, common.legend = T, 
  widths = c(0.5, 1), 
  legend = "bottom", labels = c("A", "B"), label.y = c(0.9, 0.8)
)

ggsave(plot = figS9_FineScaleMapStudyLocsAvgAgreement, 
       filename = paste0(plotDir, "AgFish_FigS9_FineScaleMapStudyLocsAvgModelAgreement.tiff"),
       width = 8,
       height = 6)

dev.off()
# end of the script













