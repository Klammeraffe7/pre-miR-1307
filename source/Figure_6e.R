"""
##_ Plotting Figure - 6e: Forest Plot of the COX-Regression Analysis_###############################
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS 13.6
Dependencies = c('R version 4.1.3 (2022-03-10)',
                 'RStudio Version 2023.06.1+524 – ©', 'dplyr_1.1.2', 'tidyr_1.3.0', 'ggplot2_3.3.6 ',
                 'survival_3.4-0','openxlsx_4.2.5','stringr_1.5.0','DT_0.25','survminer_0.4.9',
                 'data.table_1.14.2','forestmodel_0.6.2')
Description = 'This script is to perform COX-Regression Analysis on the TCGA data to produce Figure 6e 
####################################################################################################
"""

## Load Libraries
library(dplyr)
library(tidyr)
library(ggplot2)
library(survival)
library(openxlsx)
library(stringr)
library(DT)
library(survminer)
library(data.table)
library(forestmodel)
library(grid)

## Read Data

##  Read data from supplementary tables of the manuscript
# mir_rpms_median15 <- read.csv("~/miR1307/data/TCGA_BRCA_miRNA_rpm__cor_all_median15__collapsed__dup_rem.csv",
#                              sep = ",") 
# confounder <- read.csv("~/miR1307/data/confounder_miRNA_mRNA_dup_rem.csv",
#                       sep = ";")

## Reads miR1307 isomiR rpms and combines them into a dataframe

data_isomir <-cbind(mir_rpms_median15["ids"],
                    mir_rpms_median15["samp.ids"],
                    mir_rpms_median15["hsa.miR.1307.3p.0."],
                    mir_rpms_median15["hsa.miR.1307.3p.1."],
                    mir_rpms_median15["hsa.miR.1307.5p.0."])

colnames(data_isomir) <- c('ids_miRNA',
                           'sampleIDs',
                           'iso_0_3p',
                           'iso_1_3p',
                           'iso_5p')

## Integrates confounders to the miRNA expression dataset
data <- merge(data_isomir,confounder,by="ids_miRNA")

## Keeps only informative confounders and expressions
data <- select(data,c('ids_miRNA','sampleIDs','iso_0_3p','iso_1_3p',
                      'iso_5p','sample_type','subtype_mRNA','age_at_diagnosis',
                      'race','primary_diagnosis','OS.time','OS'))

colnames(data)[7] <- 'molecular_subtype'

### Removes normal samples from the data
data <- data[data$sample_type != "Normal",]
data <- data[data$molecular_subtype != "Normal_like",]

### Keep ductal and lobular types as primary diagnosis and change others to "Others"
data$primary_diagnosis <- ifelse(data$primary_diagnosis %in% c("Other, specify", 
                                                               "Mixed Histology (please specify)", 
                                                               "Mucinous Carcinoma",
                                                               "Metaplastic Carcinoma", 
                                                               "Infiltrating Carcinoma NOS",
                                                               "Medullary Carcinoma",""), 
                                 "other", data$primary_diagnosis)

## Convert age to years & remove patients without age information
data$age_at_diagnosis <- round(data$age_at_diagnosis/365.25)
data <- data[!is.na(data$age_at_diagnosis),]

# Convert overall survival information to numeric and remove missing survival time data 
colnames(data)[length(data)-1] <- 'OS_time'
data <- data[!is.na(data$OS),]
data <- data[!is.na(data$OS_time),]

data$OS <- as.numeric(data$OS)
data$OS_time <- as.numeric(data$'OS_time')

# Convert subtypes into factors
data$molecular_subtype <- as.factor(data$molecular_subtype)

# Convert primary_diagnosis information category into factors
data$primary_diagnosis <- as.factor(data$primary_diagnosis)

# Convert expression data into type numeric to make sure
data$iso_0_3p = as.numeric(data$iso_0_3p)
data$iso_1_3p = as.numeric(data$iso_1_3p)
data$iso_5p = as.numeric(data$iso_5p)

data$miR_1307_3p <- data$iso_0_3p + data$iso_1_3p
data$miR_1307_5p <- data$iso_5p

data$miR_1307_3p <- log2(data$miR_1307_3p)
data$miR_1307_5p<- log2(data$miR_1307_5p)

# Check Proportional Hazards
cox_res<- coxph(formula = Surv(OS_time, OS) ~ age_at_diagnosis  + miR_1307_3p + 
                miR_1307_5p +  molecular_subtype + primary_diagnosis, data = data)

plot(cox.zph(cox_res))
print(cox.zph(cox_res))

# Plot Figure 6f - the forest plot
cox_res2 <- coxph(formula = Surv(OS_time, OS) ~ age_at_diagnosis  + miR_1307_3p + 
                  miR_1307_5p + strata(molecular_subtype) + strata(primary_diagnosis), data = data)

cox_res_table2 <- cox_res2 %>% 
  broom::tidy(exp = FALSE) %>%
  as.data.frame(.) %>%
  dplyr::rename(coef = estimate, "se(coef)" = std.error, z = statistic) %>%
  mutate(exp_coef = exp(coef))

cox_res_table2


# pdf('Figure-COX-PH.pdf',paper='USr')
forest_model(cox_res2) + theme(axis.text.x = element_text(size=9, face="bold"),
                               plot.margin = unit(c(2.5,0.2,2.5,0.2), "cm"))
# dev.off()
sessionInfo()
