## Load Libraries

library(tidyverse)
library(ggpubr)
library(rstatix)

## Read Data
# Read miRNA rpm and confounder data from supplementary tables of the manuscript
# mir_rpms_median15 <- read.csv("~/miR1307/data/TCGA_BRCA_miRNA_rpm__cor_all_median15__collapsed__dup_rem.csv",
                              sep = ",") 
# confounder <- read.csv("~/miR1307/data/confounder_miRNA_mRNA_dup_rem.csv",
                       sep = ";")

# Figure 1B


mir_rpms_median15%>%
  select(samp.ids,
         "miR1307 3p|0" = hsa.miR.1307.3p.0.,
         "miR1307 3p|1" = hsa.miR.1307.3p.1., 
         "miR1307 5p|0" = hsa.miR.1307.5p.0.)%>%
  mutate(status = ifelse(str_detect(samp.ids,"-01"),"tumor","normal"))%>%
  pivot_longer(-c("samp.ids","status"),
               names_to = "miR",
               values_to = "rpm")%>%
  mutate(rpm = log2(rpm))%>%
  ggplot(aes(x = miR,
             y = rpm,
             fill = status))+
  geom_violin()+
  geom_boxplot(width = .1, 
               position = position_dodge(.9))+
  scale_fill_manual(values = c("#F7F7F7","#969696"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = .5),  
        panel.grid  = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 30, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"),
        axis.text.x =  element_text(colour = c("#1f8ecf","#005485","#3ab027"))
  )+
  xlab("")+
  ylab("log2(rpm)")+
  scale_y_continuous(breaks = c(1,2,3,4,5,6,7,8,9,10,11,12))+
  stat_compare_means(method = "t.test", 
                     position = "identity", 
                     label = "p.signif",
                     size = 10)

## Figure 1C



mir_rpms_median15%>%
  filter(str_detect(samp.ids,"-01"))%>%
  select(samp.ids,
         "miR1307 3p|0" = hsa.miR.1307.3p.0.,
         "miR1307 3p|1" = hsa.miR.1307.3p.1., 
         "miR1307 5p|0" = hsa.miR.1307.5p.0.) %>%
  cor_test(method = "spearman", 
           vars = c("miR1307 3p|0","miR1307 3p|1", "miR1307 5p|0"))%>%
  select(var1,var2,cor)%>%
  filter(cor != 1)%>%
  distinct(cor, .keep_all = T)%>%
  mutate(miR = paste0(var1, " vs ", var2),
         miR = (str_remove_all(miR, "miR1307 ")))%>%
  select(miR, tumor = cor) %>%
  full_join(
    mir_rpms_median15%>%
      filter(str_detect(samp.ids,"-11"))%>%
      select(samp.ids,
             "miR1307 3p|0" = hsa.miR.1307.3p.0.,
             "miR1307 3p|1" = hsa.miR.1307.3p.1., 
             "miR1307 5p|0" = hsa.miR.1307.5p.0.) %>%
      cor_test(method = "spearman", 
               vars = c("miR1307 3p|0","miR1307 3p|1", "miR1307 5p|0"))%>%
      select(var1,var2,cor)%>%
      filter(cor != 1)%>%
      distinct(cor, .keep_all = T)%>%
      mutate(miR = paste0(var1, " vs ", var2),
             miR = (str_remove_all(miR, "miR1307 ")))%>%
      select(miR, normal = cor) 
  ) %>%
  pivot_longer(names_to = "status", values_to = "correlation", cols = -miR)%>%
  ggplot(aes(y = miR, x = correlation, fill = status))+
  geom_col(position = "dodge", color = "black")+
  theme_bw()+
  scale_fill_manual(values = c("#f7f7f7","#969696"))+
  theme(plot.title = element_text(hjust = .5),  
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 30, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"))+
  xlab("spearman correlation")+
  ylab("")


## Figure 1D


mir_rpms_median15%>%
  mutate(hsa.miR.1307.3p.0. = log2(hsa.miR.1307.3p.0.),
         hsa.miR.1307.3p.1. = log2(hsa.miR.1307.3p.1.),
         hsa.miR.1307.5p.0. = log2(hsa.miR.1307.5p.0.))%>%
  left_join(confounder%>%
              select(samp.ids, subtype_mRNA))%>%
  filter(!is.na(subtype_mRNA),
         subtype_mRNA != "Normal_like")%>%
  select(samp.ids,
         hsa.miR.1307.3p.0.,
         hsa.miR.1307.3p.1., 
         hsa.miR.1307.5p.0.,
         subtype_mRNA)%>%
  mutate(Subtype = ifelse(subtype_mRNA == "Normal", "Normal","Tumour"))%>%
  mutate(sum_1307_3p = hsa.miR.1307.3p.0.+ hsa.miR.1307.3p.1.)%>%
  arrange(desc(Subtype))%>%      #plot the normal subtype above the tumor subtype for better visibility
  ggplot(aes( x= sum_1307_3p, 
              y = hsa.miR.1307.5p.0.))+
  geom_point(size = 3,
             aes(fill = Subtype), 
             shape = 21, 
             alpha = .9)+
  scale_fill_manual(values = c("#f7f7f7","#969696"))+
  theme_bw()+
  xlab("log2(rpm) miR-1307-3p")+
  ylab("log2(rpm) miR-1307-5p")+
  theme(plot.title = element_text(hjust = .5),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 30, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"))+
  guides(fill = guide_legend(override.aes = list(size = 10)))


## Figure 1E


mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., 
         hsa.miR.1307.3p.1.,  
         hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              select(samp.ids,
                     status = sample_type, 
                     subtype_mRNA))%>%
  filter(!is.na(subtype_mRNA),
         subtype_mRNA != "Normal_like")%>%
  mutate(sum_1307 = hsa.miR.1307.3p.0. +  hsa.miR.1307.3p.1.)%>%
  dplyr::select(samp.ids, 
                sum_1307, 
                hsa.miR.1307.5p.0., 
                subtype_mRNA)%>%
  pivot_longer(cols = c("sum_1307","hsa.miR.1307.5p.0."),
               names_to = "miRNA", 
               values_to = "expression")%>%
  mutate(expression = log2(expression))%>%
  mutate(miRNA = case_when(
    miRNA == "sum_1307" ~ "miR-1307-3p",
    miRNA == "hsa.miR.1307.5p.0." ~ "miR-1307-5p"
  ))%>%
  ggplot(aes(x = factor(subtype_mRNA,
                        levels = c("Normal","LumA","LumB","Her2","Basal")), 
             y = expression, 
             fill = subtype_mRNA))+
  facet_wrap(~miRNA)+
  geom_violin()+
  geom_boxplot(width = .1)+
  theme_bw()+
  theme(plot.title = element_text(hjust = .5),  
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        strip.background = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 30, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")+
  xlab("subtype")+
  ylab("log2(rpm)")+
  scale_fill_manual(values = c("#525252","#737373","#D9D9D9",
                                        "#BDBDBD","white","white"))+
                                          stat_compare_means(method = "t.test", 
                                                             label = "p.signif",
                                                             ref.group = "LumA", 
                                                             label.y = c(14,14,16,18,20),
                                                             bracket.size = 1,
                                                             size = 6.5, 
                                                             comparisons = list(c("Normal","LumA"),
                                                                                c("LumA","LumB"),
                                                                                c("LumA","Her2"),
                                                                                c("LumA","Basal")))+
  ylim(0,21)


## Figure 1F


mir_rpms_median15%>%
  dplyr::select(samp.ids,
                hsa.miR.1307.3p.0., 
                hsa.miR.1307.3p.1.,  
                hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              select(samp.ids,
                     status = sample_type, 
                     subtype_mRNA))%>%
  mutate(percent_5p = ((hsa.miR.1307.5p.0.)/
                         (hsa.miR.1307.3p.0. + hsa.miR.1307.3p.1. + hsa.miR.1307.5p.0.)
  )*100)%>%
  filter(!is.na(subtype_mRNA),
         subtype_mRNA != "Normal_like")%>%
  ggplot(aes(x = factor(subtype_mRNA,
                        levels = c("Normal","LumA","LumB","Her2","Basal")), 
             y = percent_5p, fill = subtype_mRNA))+
  geom_violin()+
  geom_boxplot(width = .1)+
  theme_bw()+
  theme(plot.title = element_text(hjust = .5),  
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 30, color = "black"),
        axis.text = element_text(colour = "black", 
                                 face = "bold"),
        legend.position = "none")+
  xlab("subtype")+
  ylab("% 5p")+
  scale_fill_manual(values = c("#525252","#737373","#D9D9D9",
                                        "#BDBDBD","white","white"))+
                                          stat_compare_means(method = "t.test", 
                                                             label = "p.signif",
                                                             ref.group = "LumA", 
                                                             label.y = c(61,61,66,71,76),
                                                             bracket.size = 1,
                                                             size = 6.5, 
                                                             comparisons = list(c("Normal","LumA"),
                                                                                c("LumA","LumB"),
                                                                                c("LumA","Her2"),
                                                                                c("LumA","Basal")))+
  ylim(0,80)


# dev.off()
sessionInfo()
