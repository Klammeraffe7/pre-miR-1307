## Load Libraries

library(tidyverse)
library(ggpubr)
library(rstatix)
library(openxlsx)

## Read Data
# Read data from supplementary tables of the manuscript

# mir_rpms_median15 <- read.csv("~/miR1307/data/TCGA_BRCA_miRNA_rpm__cor_all_median15__collapsed__dup_rem.csv",
#                              sep = ",") 

# confounder <- read.csv("~/miR1307/data/confounder_miRNA_mRNA_dup_rem.csv", 
#                        sep = ";")

# MCP_endo_data <- read.csv("~/miR1307/data/Immunedeconv.MCPCounter.UnLogTPM.BRCA.common.miRNAmRNA.TCGA.miRScoreProject.09112022.tsv",
#                            sep = "\t")%>%
#                    filter(cell_type == "Endothelial cell")%>%
#                    column_to_rownames("cell_type")%>%
#                    t()%>%
#                    as.data.frame()%>%
#                    rownames_to_column("samp.ids")%>%
#                    mutate(samp.ids = str_replace_all(samp.ids, "\\.", "-"),
#                           Endothelial_cell = `Endothelial cell`)

# endothelial_data <- read.xlsx("~/miR1307/data/TCGA_cell.type.composition.combined.xlsx"),
#                             startRow = 2)%>%
#                     dplyr::select(samp.ids = Sample_ID, endo_by_met_oct = "rel.Endothelial-cells")


## Figure 4D


#Calculate the correlation of all BRCA samples with MCP counter data

mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., hsa.miR.1307.3p.1.,  hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              dplyr::select(samp.ids, subtype_mRNA))%>%
  left_join(MCP_endo_data)%>%
  filter(subtype_mRNA != "Normal_like",
         subtype_mRNA != "Normal",
         !is.na(subtype_mRNA))%>%
  mutate(sum_1307 = (hsa.miR.1307.3p.0. + hsa.miR.1307.3p.1. + hsa.miR.1307.5p.0.))%>%
  cor_test(vars = sum_1307, vars2 = Endothelial_cell, method = "spearman")%>%
  mutate(panel = "4D BRCA all")

#Plotting

mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., hsa.miR.1307.3p.1.,  hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              dplyr::select(samp.ids, subtype_mRNA))%>%
  left_join(MCP_endo_data)%>%
  filter(subtype_mRNA != "Normal_like",
         subtype_mRNA != "Normal",
         !is.na(subtype_mRNA))%>%
  mutate(sum_1307 = (hsa.miR.1307.3p.0. + hsa.miR.1307.3p.1. + hsa.miR.1307.5p.0.))%>%
  ggplot(aes(x=sum_1307 , y = Endothelial_cell))+
  geom_point()+
  theme_bw()+
  theme(plot.title = element_text(hjust = .5),  
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 30, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.position = "none")+
  xlab("miR1307 expression")+
  ylab("relative endothelial cell content
(MCP Counter)")+
  scale_y_continuous(breaks =  c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))+
  scale_x_continuous(breaks =  c(0,2000,4000,6000,8000), limits = c(0,9000))+
  ggtitle("BRCA all")+
  geom_smooth(method = "lm", se = F)+
  annotate("text", x = 9000, y = 90, label = "R = -0.42  p =  0", size = 7, fontface = "bold", hjust = 1)




## Figure 4E



# Calculate the correlation with all BRCA samples with the methylation data

mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., hsa.miR.1307.3p.1.,  hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              dplyr::select(samp.ids, subtype_mRNA))%>%
  left_join(endothelial_data)%>%
  filter(subtype_mRNA != "Normal_like",
         subtype_mRNA != "Normal",
         !is.na(subtype_mRNA))%>%
  mutate(sum_1307 = (hsa.miR.1307.3p.0. + hsa.miR.1307.3p.1. + hsa.miR.1307.5p.0.))%>%
  cor_test(vars = sum_1307, vars2 = endo_by_met_oct, method = "spearman")%>%
  mutate(panel = "4E BRCA all")

#Plotting

mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., hsa.miR.1307.3p.1.,  hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              dplyr::select(samp.ids, subtype_mRNA))%>%
  left_join(endothelial_data)%>%
  filter(subtype_mRNA != "Normal_like",
         subtype_mRNA != "Normal",
         !is.na(subtype_mRNA))%>%
  mutate(sum_1307 = (hsa.miR.1307.3p.0. + hsa.miR.1307.3p.1. + hsa.miR.1307.5p.0.))%>%
  ggplot(aes(x = sum_1307 , y = endo_by_met_oct))+
  geom_point()+
  theme_bw()+
  theme(plot.title = element_text(hjust = .5),  
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 30, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.position = "none")+
  xlab("miR1307 expression")+
  ylab("% endothelial cells 
(methylation based)")+
  scale_y_continuous(breaks =  c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))+
  scale_x_continuous(breaks =  c(0,2000,4000,6000,8000), limits = c(0,9000))+
  ggtitle("BRCA all")+
  geom_smooth(method = "lm", se = F)+
  annotate("text", x = 9000, y = 90, label = "R = -0.41  p =  4.21e-29", size = 7, fontface = "bold", hjust = 1)

## Figure 4F


# Calculate the correlation of only basal samples with MCP counter data

mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., hsa.miR.1307.3p.1.,  hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              dplyr::select(samp.ids, subtype_mRNA))%>%
  left_join(MCP_endo_data)%>%
  filter(subtype_mRNA == "Basal")%>%
  mutate(sum_1307 = (hsa.miR.1307.3p.0. + hsa.miR.1307.3p.1. + hsa.miR.1307.5p.0.))%>%
  cor_test(vars = sum_1307, vars2 = Endothelial_cell, method = "spearman")%>%
  mutate(panel = "4F Basal only")

#Plotting

mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., hsa.miR.1307.3p.1.,  hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              dplyr::select(samp.ids, subtype_mRNA))%>%
  left_join(MCP_endo_data)%>%
  filter(subtype_mRNA == "Basal")%>%
  mutate(sum_1307 = (hsa.miR.1307.3p.0. + hsa.miR.1307.3p.1. + hsa.miR.1307.5p.0.))%>%
  ggplot(aes(x = sum_1307 , y = Endothelial_cell))+
  geom_point()+
  theme_bw()+
  theme(plot.title = element_text(hjust = .5),  
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 30, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.position = "none")+
  xlab("miR1307 expression")+
  ylab("relative endothelial cell content
(MCP Counter)")+
  scale_y_continuous(breaks =  c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))+
  scale_x_continuous(breaks =  c(0,2000,4000,6000,8000), limits = c(0,9000))+
  ggtitle("Basal only")+
  geom_smooth(method = "lm", se = F)+
  annotate("text", x = 9000, y = 90, label = "R = -0.34  p =  1.08e-5", size = 7, fontface = "bold", hjust = 1)


## Figure 4G



# Calculate the correlation with only basal samples with the methylation data

mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., hsa.miR.1307.3p.1.,  hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              dplyr::select(samp.ids, subtype_mRNA))%>%
  left_join(endothelial_data)%>%
  filter(subtype_mRNA == "Basal")%>%
  mutate(sum_1307 = (hsa.miR.1307.3p.0. + hsa.miR.1307.3p.1. + hsa.miR.1307.5p.0.))%>%
  cor_test(vars = sum_1307, vars2 = endo_by_met_oct, method = "spearman")%>%
  mutate(panel = "4G Basal only")

#Plotting

mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., hsa.miR.1307.3p.1.,  hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              dplyr::select(samp.ids, subtype_mRNA))%>%
  left_join(endothelial_data)%>%
  filter(subtype_mRNA == "Basal")%>%
  mutate(sum_1307 = (hsa.miR.1307.3p.0. + hsa.miR.1307.3p.1. + hsa.miR.1307.5p.0.))%>%
  ggplot(aes(x = sum_1307 , y = endo_by_met_oct))+
  geom_point()+
  theme_bw()+
  theme(plot.title = element_text(hjust = .5),  
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 30, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.position = "none")+
  xlab("miR1307 expression")+
  ylab("endothelial cells
(methylation-based)")+
  scale_y_continuous(breaks =  c(0,10,20,30,40,50,60,70,80,90,100), limits = c(0,100))+
  scale_x_continuous(breaks =  c(0,2000,4000,6000,8000), limits = c(0,9000))+
  ggtitle("Basal only") +
  geom_smooth(method = "lm", se = F)+
  annotate("text", x = 9000, y = 90, label = "R = -0.22  p =  0.0203", size = 7, fontface = "bold", hjust = 1)


# dev.off()
sessionInfo()
