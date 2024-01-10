# Selective arm-usage of pre-miR-1307 dysregulates angiogenesis and contributes to breast cancer aggressiveness

Repository to reproduce analyses from "Selective arm-usage of pre-miR-1307 dysregulates angiogenesis and contributes to breast cancer aggressiveness". Code can be found in the /source folder, input data are in the supplementary tables of the paper or need to be downloaded from other repositories as mentioned in the respective script.

## Author affiliations
Oyku Ece Sumer<sup>1,2,#</sup>, Korbinian Schelzig^<sup>1,2,#</sup>, Janine Jung<sup>1,2</sup>, Xiaoya Li<sup>1,3</sup>, Janina Moros<sup>1,4</sup>, Luisa Schwarzmüller<sup>1,2</sup>, Ezgi Sen<sup>1,2</sup>, Sabine Karolus<sup>1</sup>, Angelika Wörner<sup>1</sup>, Nishanth Belugali Nataraj<sup>5</sup>, Efstathios-Iason Vlachavas<sup>1</sup>, Clarissa Gerhäuser<sup>6</sup>, Karin Müller-Decker<sup>7</sup>, Dominic Helm<sup>8</sup>, Yosef Yarden<sup>5</sup>, Birgitta E. Michels<sup>1</sup> and Cindy Körner<sup>1,*</sup>

<sup>1</sup>	Division of Molecular Genome Analysis, German Cancer Research Center (DKFZ), Im Neuenheimer Feld 580, 69120 Heidelberg, Germany

<sup>2</sup>	Faculty of Biosciences, University of Heidelberg, Im Neuenheimer Feld 234, 69120 Heidelberg, Germany

<sup>3</sup>	Medical Faculty Heidelberg, University of Heidelberg, Im Neuenheimer Feld 672, 69120 Heidelberg, Germany

<sup>4</sup>	MCBI program, Department of Biology, Faculty of Science, University of Tübingen, 72074 Tübingen

<sup>5</sup>	Department of Biological Regulation, Weizmann Institute of Science, 76100 Rehovot, Israel

<sup>6</sup>	Division of Cancer Epigenomics, German Cancer Research Center (DKFZ), Im Neuenheimer Feld 280, 69120 Heidelberg, Germany

<sup>7</sup>	Tumor Models Core Facility, German Cancer Research Center (DKFZ), Im Neuenheimer Feld 280, 69120 Heidelberg, Germany

<sup>8</sup>	Proteomics Core Facility, German Cancer Research Center (DKFZ), Im Neuenheimer Feld 280, 69120 Heidelberg, Germany

<sup>#</sup>	OES and KS contributed equally to the study.

<sup>*</sup>	Correspondence: c.koerner@dkfz-heidelberg.de

## Requirements
Scripts were run in R using a variety of additional packages as mentioned in the supplementary methods of the manuscript. 

## Scripts (/source)

### TCGA-BRCA data formating
* `TCGA_BRCA_isomiR_collapse.Rmd`: Script to recapitulate the collapsing of TCGA-BRCA isomiR data to either 5p isomiR or miRNA-arm level. Requires input from GEO repository GSE164767.
* `TCGA_mRNA_htseqcount_data.R`: Script to combine TCGA-BRCA mRNA htseq-count data from the GDC portal in one table.
* `TCGA_BRCA_mRNA_removeMTgenes_TPMcalculations.Rmd`: Script to filter and format mRNA read count data and to convert them to transcripts per million (tpm).
* `TCGA_activity_scores.Rmd`: Script to calculate activity scores from the TCGA-BRCA mRNA expression data. Additionally requires download of data from the MSigDB database.

### Formating of sequencing data from overexpression of pre-miR-1307 in MDA-MB-231 cells
* `MDA_MB_231_Formating_miRNA_seq_counts_rpms.Rmd`: Script to reproduce formatting of read counts derived from small RNA sequencing of the cells. More for informational purposes. The results of this file can be downloaded from GSE227354 and used for further analyses.
* `MDA_MB_231_DESeq2.Rmd`: Script to reproduce differential expression analysis to compare cells overexpressing pre-miR-1307 and controls. Input file can be downloaded from GSE227354.

### Figure scripts
All scripts starting with "Figure" are to reproduce the indicated Figures in the manuscript. 

  


