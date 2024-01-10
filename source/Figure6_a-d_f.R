## Load Libraries

library(tidyverse)
library(ggpubr)
library(rstatix)
library(ComplexHeatmap)
library(openxlsx)

## Read Data

### Read data from supplementary tables of the manuscript

# mir_rpms_median15 <- read.csv("~/miR1307/data/TCGA_BRCA_miRNA_rpm__cor_all_median15__collapsed__dup_rem.csv",
#                              sep = ",") 

# confounder <- read.csv("~/miR1307/data/confounder_miRNA_mRNA_dup_rem.csv", 
#                       sep = ";")

# MCP_endo_data <- read.csv("~/miR1307/data/Immunedeconv.MCPCounter.UnLogTPM.BRCA.common.miRNAmRNA.TCGA.miRScoreProject.09112022.tsv",
#                          sep = "\t")%>%
#  filter(cell_type == "Endothelial cell")%>%
#  column_to_rownames("cell_type")%>%
#  t()%>%
#  as.data.frame()%>%
#  rownames_to_column("samp.ids")%>%
#  mutate(samp.ids = str_replace_all(samp.ids, "\\.", "-"),
#         Endothelial_cell = `Endothelial cell`)

# endothelial_data <- read.xlsx("~/miR1307/data/TCGA_cell.type.composition.combined.xlsx"),
# startRow = 2)%>%
 # dplyr::select(samp.ids = Sample_ID, endo_by_met_oct = "rel.Endothelial-cells")


# angiogenesis_scores <- read.csv("~/miR1307/data/ActivationScores_AllSamp.csv"),
#                                sep = ";", dec = ",")%>%
#                          select(samp.ids, ANGIOGENESIS.gene_activity_score)


# hypoxia_scores <- read.csv("~/miR1307/data/TCGA_hypoxia_scores.csv"),
#         sep = ",") 




## Figure 6A

mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., hsa.miR.1307.3p.1.,  hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              dplyr::select(samp.ids, subtype_mRNA))%>%
  left_join(MCP_endo_data)%>%
  filter(!is.na(subtype_mRNA),
         subtype_mRNA != "Normal_like")%>%
  mutate(sum_1307 = (hsa.miR.1307.3p.0. + hsa.miR.1307.3p.1. + hsa.miR.1307.5p.0.))%>%
  ggplot(aes(x = factor(subtype_mRNA,levels = c("Normal","LumA","LumB","Her2","Basal")),
             y = Endothelial_cell,
             fill = subtype_mRNA))+
  geom_violin()+
  geom_boxplot(width = .1)+
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5),  
        panel.grid.major.x = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 20, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"))+
  ylab("Endothelial cell MCP counter")+
  xlab("subtype")+
  scale_y_continuous(breaks = c(0,20,40,60,80,100))+
  scale_fill_manual(values = c("#525252","#737373","#D9D9D9","#BDBDBD","white","white"))+
  stat_compare_means(method = "t.test",label = "p.signif",paired = F,
                     label.y = c(120,110,114,118,125),
                     size = 7,vjust = 0.7,
                     bracket.size = .6,
                     comparisons = list(c("Normal","LumA"),
                                        c("LumA","LumB"),
                                        c("LumA","Her2"),
                                        c("LumA","Basal"),
                                        c("Normal","Basal")))


## Figure 6B

mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., hsa.miR.1307.3p.1.,  hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              dplyr::select(samp.ids, subtype_mRNA))%>%
  left_join(angiogenesis_scores)%>%
  filter(!is.na(subtype_mRNA),
         subtype_mRNA != "Normal_like")%>%
  ggplot(aes(x = factor(subtype_mRNA,levels = c("Normal", "LumA", "LumB","Her2" ,"Basal")),
             y = ANGIOGENESIS.gene_activity_score,
             fill = subtype_mRNA))+
  geom_violin()+
  geom_boxplot(width = .1)+
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5),  
        panel.grid.major.x = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 20, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"))+
  ylab("angiogenesis score")+
  xlab("subtype")+
  ylim(-1,1.2)+
  scale_fill_manual(values = c("#525252","#737373","#D9D9D9","#BDBDBD","white","white"))+
  stat_compare_means(method = "t.test",label = "p.signif",paired = F,
                     label.y = c(0.8,0.8,0.9,1,1.1),
                     size = 7,vjust = 0.7,
                     bracket.size = .6,
                     comparisons = list(c("Normal","LumA"),
                                        c("LumA","LumB"),
                                        c("LumA","Her2"),
                                        c("LumA","Basal"),
                                        c("Normal","Basal")))


## Figure 6C



mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., hsa.miR.1307.3p.1.,  hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              dplyr::select(samp.ids, subtype_mRNA))%>%
  left_join(hypoxia_scores)%>%
  filter(!is.na(subtype_mRNA),
         subtype_mRNA != "Normal_like")%>%
  mutate(mean_hypoxia_score = (Ragnum.Hypoxia.Score + Winter.Hypoxia.Score + Buffa.Hypoxia.Score) / 3)%>%
  ggplot(aes(x = factor(subtype_mRNA,levels = c("Normal","LumA","LumB","Her2","Basal")),
             y= mean_hypoxia_score,
             fill = subtype_mRNA))+
  geom_violin()+
  geom_boxplot(width = .1)+
  theme_bw()+
  theme(legend.position = "none",
        plot.title = element_text(hjust = .5),  
        panel.grid.major.x = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 20, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"))+
  ylab("mean hypoxia score")+
  xlab("subtype")+
  ylim(-50,60)+
  scale_fill_manual(values = c("#525252","#737373","#D9D9D9","#BDBDBD","white","white"))+
  stat_compare_means(method = "t.test",label = "p.signif",paired = F,
                     label.y = c(40,47,55),
                     size = 7,vjust = 0.7,
                     bracket.size = .6,
                     comparisons = list(c("LumA","LumB"),
                                        c("LumA","Her2"),
                                        c("LumA","Basal")))


## Figure 6D


# All BRCA samples:


mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., hsa.miR.1307.3p.1.,  hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              dplyr::select(samp.ids, subtype_mRNA))%>%
  left_join(MCP_endo_data%>%
              dplyr::select(samp.ids, Endothelial_cell))%>%
  left_join(endothelial_data%>%
              dplyr::select(samp.ids,endo_by_met_oct))%>%
  left_join(angiogenesis_scores)%>%
  distinct(samp.ids,.keep_all = T)%>%
  left_join(hypoxia_scores)%>%
  mutate(mean_hypoxia_score = (Ragnum.Hypoxia.Score + Winter.Hypoxia.Score + Buffa.Hypoxia.Score) / 3)%>%
  mutate(sum_1307_3p = hsa.miR.1307.3p.0. + hsa.miR.1307.3p.1.)%>%
  filter(!is.na(subtype_mRNA),
         subtype_mRNA != "Normal_like")%>%
  cor_test(vars = c("sum_1307_3p",
                    "hsa.miR.1307.5p.0."), 
           vars2 = c("Endothelial_cell",
                     "endo_by_met_oct",
                     "ANGIOGENESIS.gene_activity_score",
                     "mean_hypoxia_score"), 
           method = "spearman")%>%
  mutate(p.adjust = p.adjust(p, method = "BH"))%>%
  select(-c(statistic, method,p))%>%
  mutate(var1 = case_when(
    var1 == "sum_1307_3p" ~"miR-1307-3p",
    var1 == "hsa.miR.1307.5p.0." ~"miR-1307-5p"
  ),
  var2 = case_when(
    var2 == "ANGIOGENESIS.gene_activity_score" ~ "angiogenesis 
    score",
    var2 == "mean_hypoxia_score" ~ "hypoxia 
    score",
    var2 == "Endothelial_cell" ~ "[%] endothelial cell 
    MCP counter",
    var2 == "endo_by_met_oct" ~ "[%] endothelial cell
    by methylation"
  ))%>%
  ggplot(aes( x = var2, y = var1 ))+
  geom_point(aes(col = cor, size = -log10(p.adjust)))+
  theme_minimal()+
  xlab("")+
  ylab("")+
  labs(col = "spearman 
correlation")+
  theme(plot.title = element_text(hjust = .5),  
        panel.border = element_blank(),
        text = element_text(size = 20, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.box = "horizontal")+
  ggtitle("all BRCA patients")+
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab" , limits = c(-.6,.6))+
  scale_size_continuous(range = c(1,20),
                        breaks = c(10,20,40,60,80,100)) -> fig6d_all

#Same correlation but LumA patients only

mir_rpms_median15%>%
  select(samp.ids,hsa.miR.1307.3p.0., hsa.miR.1307.3p.1.,  hsa.miR.1307.5p.0.)%>%
  left_join(confounder%>%
              dplyr::select(samp.ids, subtype_mRNA))%>%
  left_join(MCP_endo_data%>%
              dplyr::select(samp.ids, Endothelial_cell))%>%
  left_join(endothelial_data%>%
              dplyr::select(samp.ids,endo_by_met_oct))%>%
  left_join(angiogenesis_scores)%>%
  distinct(samp.ids, .keep_all = T)%>%
  left_join(hypoxia_scores)%>%
  mutate(mean_hypoxia_score = (Ragnum.Hypoxia.Score + Winter.Hypoxia.Score + Buffa.Hypoxia.Score) / 3)%>%
  mutate(sum_1307_3p = hsa.miR.1307.3p.0. + hsa.miR.1307.3p.1.)%>%
  filter(subtype_mRNA == "LumA")%>%
  cor_test(vars = c("sum_1307_3p",
                    "hsa.miR.1307.5p.0."), 
           vars2 = c("Endothelial_cell",
                     "endo_by_met_oct",
                     "ANGIOGENESIS.gene_activity_score",
                     "mean_hypoxia_score"), 
           method = "spearman")%>%
  mutate(p.adjust = p.adjust(p, method = "BH"))%>%
  select(-c(statistic, method,p))%>%
  mutate(var1 = case_when(
    var1 == "sum_1307_3p" ~"miR-1307-3p",
    var1 == "hsa.miR.1307.5p.0." ~"miR-1307-5p"
  ),
  var2 = case_when(
    var2 == "ANGIOGENESIS.gene_activity_score" ~ "angiogenesis 
    score",
    var2 == "mean_hypoxia_score" ~ "hypoxia 
    score",
    var2 == "Endothelial_cell" ~ "[%] endothelial cell 
    MCP counter",
    var2 == "endo_by_met_oct" ~ "[%] endothelial cell
    by methylation"
  ))%>%
  ggplot(aes( x = var2, y = var1 ))+
  geom_point(aes(col = cor, size = -log10(p.adjust)))+
  theme_minimal()+
  xlab("")+
  ylab("")+
  labs(col = "spearman 
correlation")+
  theme(plot.title = element_text(hjust = .5),  
        panel.border = element_blank(),
        text = element_text(size = 20, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.box = "horizontal")+
  ggtitle("LumA patients only")+
  scale_color_gradient2(midpoint = 0, low = "blue", mid = "white",
                        high = "red", space = "Lab", limits = c(-.6,.6) )+
  scale_size_continuous(range = c(1,20), breaks = c(10,20,40,60,80,100)) -> fig6d_luma

# add the plots together:
ggpubr::ggarrange(fig6d_all, fig6d_luma, common.legend = T,ncol = 1, legend = "right")



## Figure6F

#scatter the correlation vs mean %

mir_rpms_median15%>%
  pivot_longer(cols = -c("samp.ids","ids"), names_to = "miRNA",values_to = "expression")%>%
  mutate(miR = sub("^(?:[^.]*\\.){2}(.*)", "\\1", miRNA))%>%
  separate(miR, into = c("miR","arm","subtype","rest"))%>%
  mutate(arm = ifelse(rest != "",subtype,arm))%>%
  left_join(confounder%>%dplyr::select(samp.ids, subtype_mRNA))%>%
  filter(subtype_mRNA != "Normal_like",
         !is.na(subtype_mRNA))%>%
  mutate(Subtype = ifelse(str_detect(subtype_mRNA,"Normal"),"Normal","Tumor"))%>%
  group_by(samp.ids,Subtype,miR, arm)%>%
  summarize(summed_expression = sum(expression, na.rm = T))%>%
  ungroup()%>%
  filter(!is.na(arm))%>%
  pivot_wider(names_from = arm, values_from = summed_expression)%>%
  filter(!is.na(`5p`),
         !is.na(`3p`))%>%
  mutate(percent_5p =    `5p`/(  `5p`+ `3p`)*100)%>%
  mutate(miR = paste0("miR_",miR))%>%
  dplyr::select(samp.ids,Subtype,miR,percent_5p)%>%
  pivot_wider(names_from = miR, values_from = percent_5p)%>%
  dplyr::select(-samp.ids)%>%
  group_by(Subtype)%>%
  cor_test(vars = miR_1307, method = "spearman") -> correlation_5p_percent

#Create a percentage matrix

mir_rpms_median15%>%
  pivot_longer(cols = -c("samp.ids","ids"), names_to = "miRNA",values_to = "expression")%>%
  mutate(miR = sub("^(?:[^.]*\\.){2}(.*)", "\\1", miRNA))%>%
  separate(miR, into = c("miR","arm","subtype","rest"))%>%
  mutate(arm = ifelse(rest != "",subtype,arm))%>%
  left_join(confounder%>%dplyr::select(samp.ids, subtype_mRNA))%>%
  filter(subtype_mRNA != "Normal_like",
         !is.na(subtype_mRNA))%>%
  mutate(Subtype = ifelse(str_detect(subtype_mRNA,"Normal"),"Normal","Tumor"))%>%
  group_by(samp.ids,miR, arm,)%>%
  summarize(summed_expression = sum(expression, na.rm = T))%>%
  ungroup()%>%
  filter(!is.na(arm))%>%
  pivot_wider(names_from = arm, values_from = summed_expression)%>%
  filter(!is.na(`5p`),
         !is.na(`3p`))%>%
  mutate(percent_5p =    `5p`/(  `5p`+ `3p`)*100)%>%
  dplyr::select(miR,percent_5p,samp.ids)%>%
  arrange(desc(samp.ids))%>%
  pivot_wider(names_from = samp.ids, values_from = percent_5p)%>%
  mutate(miR = paste0("miR-",miR))%>%
  arrange(desc(miR))%>%
  column_to_rownames("miR")%>%
  as.matrix()%>%
  t() -> percentage_matrix_5p

percentage_matrix_5p[1:5,1:5]

#Scaling of the matrix

percentage_matrix_5p%>%scale() -> percentage_matrix_5p_scaled

#Annotations for the heatmap

col_annotations <- ComplexHeatmap::HeatmapAnnotation(which = "col",
                                                     border = F,
                                                     annotation_name_side = "left",
                                                     show_legend = F,
                                                     col = list(Type = c("Normal" = "#f7f7f7",
                                                                         "Tumor" = "#969696")
                                                                ),
                                                     "Type" = confounder%>%
                                                       dplyr::select(samp.ids,subtype_mRNA)%>%
                                                       filter(samp.ids %in% rownames(percentage_matrix_5p))%>%
                                                       arrange(desc(samp.ids))%>%
                                                       mutate(Type = ifelse(subtype_mRNA == "Normal","Normal","Tumor"))%>%
                                                       pull(Type)
)

row_annotations <- ComplexHeatmap::HeatmapAnnotation(which = "row",
                                                     border = T,
                                                     "Correlation Tumor
with % miR-1307-5p" = ComplexHeatmap::anno_simple(
  correlation_5p_percent%>%
    filter(Subtype == "Tumor")%>%
    ungroup()%>%
    dplyr::select(var2, cor)%>%
    mutate(var2 = str_remove_all(var2,"miR_"))%>%
    mutate(var2 = paste0("miR-",var2))%>%
    filter(var2 %in% colnames(percentage_matrix_5p))%>%
    bind_rows(tibble(var2 = "miR-1307",cor = 1))%>%
    arrange(desc(var2))%>%
    pull("cor"),
  col = circlize::colorRamp2(c(-1, -.5,0,.5, 1), c("#8c510a","#d8b365", "#f5f5f5","#5ab4ac", "#01665e")),
  pch = correlation_5p_percent%>%
    mutate(p.adj = p.adjust(p, method = "BH"))%>%
    filter(Subtype =="Tumor")%>%
    mutate(var2 = str_replace_all(var2,"_","-"))%>%
    mutate(pch = case_when(p.adj <= 0.01 ~ "*", .default = NA))%>%
    dplyr::select(var2,pch)%>%
    bind_rows(tibble(var2 = "miR-1307",pch = NA) )%>%
    arrange(desc(var2))%>%
    pull(pch)))

#create the legend

lgd1<- ComplexHeatmap::Legend(title_gp = gpar(fontsize = 15, fontface = 'bold'),
                              title = "z score % 5p",
                              at = c(-4, -2,0,2, 4),
                              labels = c("-4","-2", "0","2", "4"),
                              col_fun  = circlize::colorRamp2(c(-4,0, 4), c('blue','white','red'))
)


lgd2<- ComplexHeatmap::Legend(title_gp = gpar(fontsize = 15, fontface = 'bold'),
                              title = "Type",
                              labels = c('Normal', 'Tumor'),
                              legend_gp = gpar(fill= c('#f7f7f7','#969696'), color = 'black')
)

lgd3<- ComplexHeatmap::Legend(title_gp = gpar(fontsize = 15, fontface = 'bold'),
  title = "Correlation",
  at = c(-1, -.5,0,.5, 1),
  labels = c("-1","-0.5", "0","0.5", "1"),
  col_fun  = circlize::colorRamp2(c(-1, -.5,0,.5, 1), c("#8c510a","#d8b365", "#f5f5f5","#5ab4ac", "#01665e"))
)

legend_param <- packLegend(lgd1, lgd2,lgd3)

#Plot the heatmap
ComplexHeatmap::draw(
  percentage_matrix_5p_scaled%>%t()%>%
    ComplexHeatmap::Heatmap(km = 2, column_km = 2, 
                            name = "z-score % 5p",
                            row_title = "miRNAs", show_row_names = T,row_title_side = "left",
                            row_title_gp = gpar(fontsize = 20),
                            column_title = "patients", show_column_names = F, column_title_side = "bottom",
                            column_title_gp = gpar(fontsize = 20),
                            right_annotation = row_annotations,
                            top_annotation = col_annotations,show_heatmap_legend = F
    ),
  annotation_legend_list = legend_param
)

