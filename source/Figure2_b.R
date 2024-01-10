## Load libraries

library(tidyverse)
library(openxlsx)
library(ggrepel)

## Load in the data
### Download DESeq2 result data from GEO repository (GSE227354) and use as input
## DESeq_data_1307_vs_ctrl <- read.xlsx("~/miR1307/data/20220506_miR1307_vs_wt_5iso_DESeq2.xlsx")

# Figure 2B

ggplot()+
  geom_point(data = DESeq_data_1307_vs_ctrl%>%                              #add all points
               filter(-log10(padj_miR1307_vs_ctrl) < -log10(0.05)),
             aes(x = log2FoldChange_miR1307_vs_ctrl , 
                 y = -log10(padj_miR1307_vs_ctrl)), 
             col = "grey", 
             size = 5,
             alpha = .75)+
  geom_point(data = DESeq_data_1307_vs_ctrl%>%                               #color significant points black
               filter(-log10(padj_miR1307_vs_ctrl) > -log10(0.05)),
             aes(x = log2FoldChange_miR1307_vs_ctrl , 
                 y = -log10(padj_miR1307_vs_ctrl)), 
             col = "black", 
             size = 5, 
             alpha = .75)+
  geom_point(data = DESeq_data_1307_vs_ctrl%>%                             #color isomiRs of 1307
               filter(str_detect(ids,"1307-3p\\|1")),
             aes(x = log2FoldChange_miR1307_vs_ctrl , 
                 y = -log10(padj_miR1307_vs_ctrl)),
             col = "#005485", 
             size = 6)+
  geom_point(data = DESeq_data_1307_vs_ctrl%>%
               filter(str_detect(ids,"1307-3p\\|0")),
             aes(x = log2FoldChange_miR1307_vs_ctrl , 
                 y = -log10(padj_miR1307_vs_ctrl)),
             col = "#1F8ECF", 
             size = 6)+
  geom_point(data = DESeq_data_1307_vs_ctrl%>%
               filter(str_detect(ids,"1307-5p\\|0")),
             aes(x = log2FoldChange_miR1307_vs_ctrl , 
                 y = -log10(padj_miR1307_vs_ctrl)),
             col = "#3AB027", 
             size = 6)+
  geom_text_repel(                                             #make the text annotations
    data = DESeq_data_1307_vs_ctrl%>%
      filter(str_detect(ids,"1307"),
             !str_detect(ids,"1307-3p\\|2")),
    aes(label = ids,
        x = log2FoldChange_miR1307_vs_ctrl , 
        y = -log10(padj_miR1307_vs_ctrl)),
    size = 8,
    box.padding = unit(2, "lines"),
    point.padding = unit(1.5, "lines")
  )+
  theme_bw()+
  theme(plot.title = element_text(hjust = .5),  
        panel.grid = element_blank(), 
        panel.border = element_blank(),
        axis.line = element_line(color = "black"),
        text = element_text(size = 30, color = "black"),
        axis.text = element_text(colour = "black", face = "bold"),
        legend.position = "none")+
  xlab("log2FoldChange")+
  ylab("-log10(p.adjust)")+
  ggtitle("miR1307 vs ctrl")+
  geom_hline(yintercept = -log10(0.05),linetype="dotted")+
  scale_y_continuous(breaks = c(0,2,4,6,8,10,12,14))
