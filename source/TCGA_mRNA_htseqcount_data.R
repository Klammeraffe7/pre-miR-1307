### Obtaining data frame with read counts of TCGA mRNA data

## libraries -----
library(plyr)
library(dplyr)
library(data.table)
library(biomaRt)
library(tidyr)

## directories -----
data_dir = "~/tcga_data/"
out_dir = "~/tcga_data/HTSeq_count_data/"


projects = c("TCGA-BRCA")

## iterate through projects -----
lapply(seq_along(projects), function(x){
 Project = projects[x]
 # the mapping file allows to find the patient ID that corresponds to the cryptic file name; to be downloaded from GDC portal
 mapping_df = read.csv(file=paste0(data_dir,Project,"/results_gene_expression.csv"), stringsAsFactors = F) 
 # define folders containing raw data for the respective TCGA cohort; 
 # files containing read counts for each genes determined by htseq-count to be downloaded from GDC portal
 folders = list.dirs(paste0(data_dir,Project,"/harmonized/Transcriptome_Profiling/Gene_Expression_Quantification"))
 folders = folders[-1]
 # folders <- folders[1:10] # for test purposes only
 # extract readcounts from htseq-count files
 reads<- lapply(seq_along(folders), function(z){
   if (length(list.files(folders[z],pattern = "*.htseq.counts.gz$"))>=1){
     path = list.files(folders[z],pattern = "*.htseq.counts.gz$", full.names=TRUE)
     file = list.files(folders[z],pattern = "*.htseq.counts.gz$")
     pid = mapping_df$cases[mapping_df$file_name == file]
     exprs_file = fread(path, data.table= FALSE) # read htseq-counts
     names(exprs_file) <- c("id", pid)
     return(exprs_file)
   }
 })
 
 reads <- reads[lengths(reads) > 0L] # remove empty entries from reads
 reads_df = do.call(cbind, reads) # make dataframe from reads list
 
 reads_df <- reads_df[-c((nrow(reads_df)-4):nrow(reads_df)),-seq(3, ncol(reads_df), by = 2)] # remove unneeded rows(htseq summary in the end of the file) and columns (duplications of the id column)
 print(dim(reads_df))
 
 ## Annotations -----
 # add HGNC symbols as annotations
 Separation = reads_df[,1:2]
 Separation$ids <- reads_df$id
 head(Separation)
 Separation <- separate(data = Separation, col = id, into = c("id", "notNeeded"), sep = "\\.")
 ForAnnotation<- Separation[,c(1,4)]
 head(ForAnnotation)
 ForAnnotation <- ForAnnotation[!duplicated(ForAnnotation$id), ]
 row.names(ForAnnotation)<- ForAnnotation$id
 ForAnnotation$id<- NULL
 head(ForAnnotation)
 
 # Human orthologues Annotations
 HumanAnnotations<- getBM(filters= "ensembl_gene_id",
                          attributes= c("ensembl_gene_id",
                                        "hgnc_symbol"),
                          values=row.names(ForAnnotation),
                          mart= useDataset("hsapiens_gene_ensembl", useMart("ensembl")))
 HumanAnnotations[1:10,]
 HumanAnnotations<-as.data.frame(HumanAnnotations)
 nrow(HumanAnnotations)
 
 idx1 <- match(row.names(ForAnnotation), HumanAnnotations$ensembl_gene_id)
 Annotations<-cbind(ForAnnotation, HumanAnnotations[idx1, ])
 head(Annotations)
 Annotations["ensembl_gene_id"] <- NULL
 names(Annotations)[1] = "id"
 reads_df <- join(Annotations, reads_df, type = "right")
 
 write.csv2(reads_df, paste0(out_dir,  "/htseq_count_", Project, ".csv"), row.names = FALSE)
})
