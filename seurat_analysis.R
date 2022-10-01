######################################  
## install packages
######### ---------------------------------------------------
##
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("glmGamPoi") 
######### ---------------------------------------------------
##
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages("remotes")
# }
# remotes::install_github("mojaveazure/seurat-disk")
######### ---------------------------------------------------
##
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("SingleR", force = TRUE)
######### ---------------------------------------------------
##
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("scRNAseq")
######### ---------------------------------------------------
##
# install.packages("celldex")
# install.packages("remotes")
# remotes::install_github("LTLA/celldex")
######### ---------------------------------------------------
# BiocManager::install("scater")
# BiocManager::install("harmony")
###################################### 




######################################
## libraries
######### ---------------------------------------------------
library(Seurat)
library(magrittr)
library(harmony)
library(dplyr)
library(fs)
#library("glmGamPoi")
## my libraries 
library(SeuratDisk)
library(tidyverse)
library(readr)
# reticulate::py_install(packages = 'umap-learn')
library(Matrix)
library(sjmisc)
library(ggplot2)
##
library(SingleR)
#library(scRNAseq)
library(celldex)
#library(scater)
###################################### 



############## SEURAT
## setting work dir

work_dir = file.path("scRNA_data")
file_dir = file.path(work_dir)
file_dir
plot_dir = file.path(work_dir, "output")
plot_dir

setwd(file_dir)

## download date from
## https://data.humancellatlas.org/explore/projects/a62dae2e-cd69-4d5c-b5f8-4f7e8abdbafa
## https://data.humancellatlas.org/explore/projects/a62dae2e-cd69-4d5c-b5f8-4f7e8abdbafa/project-matrices

## reading in the 10X cell ranger .H5 output
## works with library(SeuratDisk)
h5_object = Read10X_h5(
  # filename = "output_filtered_hw1.h5", 
  filename = "GSM4186975_HTAPP-394-SMP-1561_CST_channel1_raw_gene_bc_matrices_h5.h5",
  use.names = TRUE, 
  unique.features = TRUE)

str(h5_object)
head(h5_object, n=3)


seurat_h5 = CreateSeuratObject(counts = h5_object)
str(seurat_h5)
seurat_h5

# An object of class Seurat 
# 33694 features across 737280 samples within 1 assay 
# Active assay: RNA (33694 features, 0 variable features)


## counts data 
## it should be a sparse matrix
count = seurat_h5$nCount_RNA
head(count)
## barcode
barcode = seurat_h5$orig.ident
head(barcode)
## feature
feature = seurat_h5$nFeature_RNA
head(feature)


###### 1. QC --------------------------------------------------------------
# percent mitochondrial reads


seurat_h5[["percent.mt"]] = PercentageFeatureSet(seurat_h5, pattern = "^MT-")
# View(seurat_h5@meta.data)

print(head(seurat_h5[["percent.mt"]], n=5))
head(seurat_h5@meta.data, 5)

# Visualize QC metrics as a violin plot
violin_plot = 
  VlnPlot(seurat_h5, 
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
          ncol = 3)


violin_plot



ggsave(
  file.path(plot_dir, violin_plot),
  plot = last_plot(),
  device = NULL,
  path = NULL,
  scale = 1,
  width = NA,
  height = NA,
  units = c("in", "cm", "mm", "px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL
)




# ggsave(
#   violin_plot,
#   plot = last_plot(),
#   device = NULL,
#   path = NULL,
#   scale = 1,
#   width = NA,
#   height = NA,
#   units = c("in", "cm", "mm", "px"),
#   dpi = 300,
#   limitsize = TRUE,
#   bg = NULL,
#   ...
# )




# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
plotA = 
  FeatureScatter(seurat_h5, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotB = 
  FeatureScatter(
    seurat_h5, 
    feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

plotA + plotB


####### 2. Filtering --------------------------------------------------------
seurat_h5_subset = 
  subset(seurat_h5, 
         subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
seurat_h5_subset

# An object of class Seurat 
# 33694 features across 6367 samples within 1 assay 
# Active assay: RNA (33694 features, 0 variable features)


# An object of class Seurat 
# 38436 features across 538 samples within 1 assay 
# Active assay: RNA (38436 features, 0 variable features)


plotA1 = 
  FeatureScatter(seurat_h5_subset, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotB1 = 
  FeatureScatter(
    seurat_h5_subset, 
    feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

plotA1 + plotB1



####### 3. Normalize Data ---------------------------------------------------
## ## normalize
seurat_h5_normalised = NormalizeData(seurat_h5_subset)
seurat_h5_normalised

# An object of class Seurat 
# 33694 features across 6367 samples within 1 assay 
# Active assay: RNA (33694 features, 0 variable features)


# An object of class Seurat 
# 38436 features across 538 samples within 1 assay 
# Active assay: RNA (38436 features, 0 variable features)

str(seurat_h5_normalised)


plotA2 = 
  FeatureScatter(seurat_h5_normalised, feature1 = "nCount_RNA", feature2 = "percent.mt")
plotB2 = 
  FeatureScatter(
    seurat_h5_normalised, 
    feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')

plotA2 + plotB2




####### 4. Identification of highly variable features (feature selection) ------
### find variable features

seurat_h5_hvf = 
  FindVariableFeatures(
    seurat_h5_normalised, 
    selection.method = "vst", 
    nfeatures = 2000)

seurat_h5_hvf

# An object of class Seurat 
# 33694 features across 6367 samples within 1 assay 
# Active assay: RNA (33694 features, 2000 variable features)


# An object of class Seurat 
# 38436 features across 538 samples within 1 assay 
# Active assay: RNA (38436 features, 2000 variable features)

class(seurat_h5_hvf)
# [1] "Seurat"
# attr(,"package")
# [1] "SeuratObject"

top10_hvf = head(VariableFeatures(seurat_h5_hvf), 10)

top10_hvf

# [1] "SERPINE1"     "PCDH9"        "PPP2R2B"      "MUCL1"        "RP11-281J9.2" "LRP1B"       
# [7] "COL1A2"       "SCGB2A2"      "FN1"          "CXCL8"     


# [1] "ADAMTS9"   "MCTP1"     "PTPRC"     "AKAP12"    "MEDAG"     "LINC00278" "DCN"       "CFD"      
# [9] "FLT1"      "PHLDB2"  

# plot variable features with and without labels
plot1 = VariableFeaturePlot(seurat_h5_hvf)
plot2 = LabelPoints(plot = plot1, points = top10_hvf, repel = TRUE)
plot1 + plot2


####### 5. Scaling ---------------------------------------------------------
# Scaling the data

all.genes = rownames(seurat_h5_hvf)
all.genes %>% head()
# [1] "MIR1302-2HG" "FAM138A"     "OR4F5"       "AL627309.1"  "AL627309.3"  "AL627309.2" 


seurat_h5_scaled = ScaleData(seurat_h5_hvf)
seurat_h5_scaled

# An object of class Seurat 
# 33694 features across 6367 samples within 1 assay 
# Active assay: RNA (33694 features, 2000 variable features)

# An object of class Seurat 
# 38436 features across 538 samples within 1 assay 
# Active assay: RNA (38436 features, 2000 variable features)

seurat_h5_scaled_centered = 
  ScaleData(
    seurat_h5_scaled, 
    vars.to.regress = "percent.mt"
  )
seurat_h5_scaled_centered
# An object of class Seurat 
# 33694 features across 6367 samples within 1 assay 
# Active assay: RNA (33694 features, 2000 variable features)



# An object of class Seurat 
# 38436 features across 538 samples within 1 assay 
# Active assay: RNA (38436 features, 2000 variable features)


####### 6. Perform linear dimensional reduction -----------------------------
# Perform linear dimensional reduction
seurat_h5_reduced_PCA = 
  RunPCA(
    seurat_h5_scaled_centered, 
    features = VariableFeatures(object = seurat_h5_scaled_centered)
  )
seurat_h5_reduced_PCA
# An object of class Seurat 
# 33694 features across 6367 samples within 1 assay 
# Active assay: RNA (33694 features, 2000 variable features)
# 1 dimensional reduction calculated: pca






# An object of class Seurat 
# 38436 features across 538 samples within 1 assay 
# Active assay: RNA (38436 features, 2000 variable features)
# 1 dimensional reduction calculated: pca

# # A Seurat-tibble abstraction: 538 × 10
# # Features=38436 | Cells=538 | Active assay=RNA | Assays=RNA
# .cell              orig.ident    nCount_RNA nFeature_RNA percent.mt    PC_1   PC_2   PC_3   PC_4     PC_5
# <chr>              <fct>              <dbl>        <int>      <dbl>   <dbl>  <dbl>  <dbl>  <dbl>    <dbl>
#   1 AGACAAACAGTTCCAA-1 SeuratProject       2923         1729      0.958  -5.92  -27.2    7.55   7.99 -19.0   
# 2 AGTTAGCTCCTACACC-1 SeuratProject       4080         2219      0.392 -26.3     7.76  23.8   13.9    5.08  
# 3 GTGTTCCTCTTGTTAC-1 SeuratProject       4240         2443      2.41  -35.7    -7.72  -9.34   7.69  -0.625 
# 4 GTCTGTCAGGGCAAGG-1 SeuratProject       3245         2256      4.68  -24.6   -12.8   -5.52   2.43  -0.0312
# 5 GCATGATGTTGGCCGT-1 SeuratProject       3508         2329      0.884 -26.7     6.32  18.4    7.21   2.10  
# 6 ACTTATCCAGTGTGCC-1 SeuratProject       3734         2442      0.107   1.71  -22.7    1.63  -6.45   7.32  
# 7 TCGTCCAAGCGCAATG-1 SeuratProject       3752         2409      3.78  -32.9   -12.1   -8.54   3.85  -0.0670
# 8 TTTCACAGTTTCGCTC-1 SeuratProject       3365         2271      0.267 -17.2    15.6   17.9  -17.8   -5.59  
# 9 TATCAGGTCTGGAAGG-1 SeuratProject       3792         2326      2.85  -36.9    -9.39 -12.6    2.27   0.239 
# 10 AGGATAAAGAGCAAGA-1 SeuratProject       3734         1911      0.348   0.928 -32.2    9.26   2.02 -16.0   
# # … with 528 more rows
# # ℹ Use `print(n = ...)` to see more rows

# PC_ 1 
# Positive:  MECOM, FAM155A, MYRIP, ANO2, ST6GALNAC3, MCTP1, EMCN, CCSER1, PLCB4, EVA1C 
# TLL1, VWF, PTPRB, NRG3, ADGRL4, RYR3, CDH13, AL589693.1, LINC02147, PLEKHG1 
# CNTNAP3B, DISC1FP1, STOX2, CPXM2, ERG, EGFL7, RALGAPA2, FLT1, SLCO2A1, RAPGEF4 
# Negative:  CD44, SVEP1, LAMA2, EGFR, NOVA1, COL12A1, DLC1, COL5A2, REV3L, ANK2 
# ABCA9, DCLK1, MAN1A1, RTN4, EBF1, TNFAIP6, GFPT2, CELF2, PRRX1, COL6A2 
# ABCA10, PCDH9, PLAGL1, SLC19A2, MEDAG, CCDC80, NEGR1, XIST, ANXA1, FGFR1 
# PC_ 2 
# Positive:  MT-RNR2, PDE3B, TMEM132C, SLC19A3, PCDH9, TENM3, ACSS3, MLXIPL, AC004160.1, ADH1B 
# LPL, GRK3, BCL2, PDZRN3, LGR4, TRHDE-AS1, ITSN1, ACSL1, AQP7, LRP1B 
# SOX5, GHR, GYG2, GABRE, ANGPTL4, GPAM, ACVR1C, CPM, AL136119.1, LINC01239 
# Negative:  FLT1, MCTP1, ADAMTS9, MIR222HG, KLF7, ERG, EVA1C, ADGRL4, PKP4, PLEKHG1 
# PODXL, ARL15, PDLIM1, ICAM1, EGFL7, LDB2, ADAMTS1, ZNF385D, DOCK4, TCF4 
# RHOJ, CDK17, VEGFC, SYNE2, ABLIM1, MSN, BMPR2, FLI1, RALGAPA2, MECOM 
# PC_ 3 
# Positive:  PDE3B, SORBS1, DIRC3, TMEM131L, LPL, CD36, LIPE-AS1, ATP1B3, MAST4, SLC19A3 
# BCL2, ACACB, PHLDB2, SAT1, PPARG, PPP2R1B, FABP4, WDPCP, SOS1, TNS1 
# GPAM, MLXIPL, CRIM1, PDZRN3, RTN3, LIPE, RANBP2, EEPD1, OGA, LENG8 
# Negative:  DCLK1, KAZN, NEGR1, LINC02511, BICC1, MEDAG, GFPT2, COL1A2, ROBO2, CFH 
# UST, LAMA2, CHL1, DCN, LPAR1, TNFAIP6, CACNA2D3, NOVA1, ABI3BP, C3 
# GPRC5A, MAP1B, MEG8, AL110292.1, PRRX1, TNXB, VCAN, CCDC80, KCND2, SLIT2 
# PC_ 4 
# Positive:  ATP13A3, ACSL4, MCL1, NAMPT, ZNF331, ATP1B3, GPRC5A, COQ10B, TSC22D2, ELL2 
# SLC4A7, CLIC4, ELOVL5, IL6R, NR4A3, CREM, ANXA1, LDLR, TMED5, THBS1 
# MARCH3, EMP1, OSMR, NABP1, TUBB6, RTN4, NUP98, VDAC1, MEDAG, SLC38A2 
# Negative:  IGF1, PTEN, PIK3R1, GAB1, HSPA1A, ADH1B, FGF10, RHOBTB3, LINC00278, TRHDE-AS1 
# PTGER3, FAM13A, PECR, IRS2, CACHD1, IGFBP5, MME, ORC2, LINC01239, MAPK14 
# PLSCR4, BCL6, ACVR1C, ANTXR2, ABCA8, SOS1, LRP1B, GRK3, CDON, FBLN5 
# PC_ 5 
# Positive:  IL7R, PTPRC, LINC02147, CYTIP, SH2B3, BCL11B, ST6GAL1, ID1, ATP2B1, NRG3 
# SAMSN1, ZNF831, PARP8, NR2F2-AS1, DOCK2, PTGIS, HSPA14.1, SKAP1, MEOX2, PCED1B-AS1 
# NOSTRIN, ARHGAP15, AL589693.1, DGKA, HEG1, PGM5, MYRIP, ESYT2, CDC42SE2, SPNS2 
# Negative:  ADAMTS4, CSRP2, PNP, IL6, PKIG, RND1, AKAP12, RASGRF2, UGCG, TUSC3 
# GCH1, MAST4, TM4SF1, SERPINB1, KDR, NFKBIA, RNF122, LINC01695, COL4A1, PARP14 
# AL583785.1, ATP11A, ASPH, RELL1, CXCL2, ADAMTS9, ETS2, BACE2, ADAMTS1, UBASH3B 



seurat_h5_reduced_CCA = 
  RunCCA(object1 = seurat_h5_scaled_centered, 
         object2 = seurat_h5_scaled_centered,
  )

seurat_h5_reduced_CCA
# An object of class Seurat 
# 33694 features across 12734 samples within 1 assay 
# Active assay: RNA (33694 features, 0 variable features)
# 1 dimensional reduction calculated: cca







# # A Seurat-tibble abstraction: 1,076 × 10
# # Features=38436 | Cells=1076 | Active assay=RNA | Assays=RNA
# .cell                orig.ident    nCount_RNA nFeature_RNA percent.mt    CC_1     CC_2    CC_3     CC_4     CC_5
# <chr>                <chr>              <dbl>        <int>      <dbl>   <dbl>    <dbl>   <dbl>    <dbl>    <dbl>
#   1 AGACAAACAGTTCCAA-1_1 SeuratProject       2923         1729      0.958  0.0398  0.0884   0.0572  0.0746   0.144  
# 2 AGTTAGCTCCTACACC-1_1 SeuratProject       4080         2219      0.392 -0.0945 -0.00541  0.140   0.110   -0.105  
# 3 GTGTTCCTCTTGTTAC-1_1 SeuratProject       4240         2443      2.41  -0.101   0.109   -0.0392  0.0525  -0.00464
# 4 GTCTGTCAGGGCAAGG-1_1 SeuratProject       3245         2256      4.68  -0.0554  0.110   -0.0201  0.00924 -0.0124 
# 5 GCATGATGTTGGCCGT-1_1 SeuratProject       3508         2329      0.884 -0.104   0.00614  0.128   0.0654  -0.0730 
# 6 ACTTATCCAGTGTGCC-1_1 SeuratProject       3734         2442      0.107  0.0778  0.0863   0.0199 -0.0652  -0.0964 
# 7 TCGTCCAAGCGCAATG-1_1 SeuratProject       3752         2409      3.78  -0.0801  0.116   -0.0353  0.0226  -0.0109 
# 8 TTTCACAGTTTCGCTC-1_1 SeuratProject       3365         2271      0.267 -0.0868 -0.0625   0.107  -0.134    0.0479 
# 9 TATCAGGTCTGGAAGG-1_1 SeuratProject       3792         2326      2.85  -0.0903  0.104   -0.0555  0.0113   0.00732
# 10 AGGATAAAGAGCAAGA-1_1 SeuratProject       3734         1911      0.348  0.0853  0.103    0.0807  0.0410   0.158  
# # … with 1,066 more rows
# # ℹ Use `print(n = ...)` to see more rows



# Examine and visualize PCA results a few different ways
print(seurat_h5_reduced_PCA[["pca"]], dims = 1:5, nfeatures = 5)
# PC_ 1 
# Positive:  SLC4A8, DACH1, MAGI3, C8orf4, ESR1 
# Negative:  PLXDC2, MAML2, SLC1A3, FMNL2, ZEB2 
# PC_ 2 
# Positive:  PCDH9, PPP2R2B, PCDH9-AS2, EDIL3, SLC24A2 
# Negative:  PLXDC2, MSR1, RGS1, CELF2, APBB1IP 
# PC_ 3 
# Positive:  TOP2A, DIAPH3, ARHGAP11B, ASPM, MELK 
# Negative:  C8orf4, RP11-624L4.1, CYP4Z1, SNX9, DLG5 
# PC_ 4 
# Positive:  TPRG1, KCNE4, PTPRT, NOVA1, CLSTN2 
# Negative:  COL4A1, COL4A2, IGFBP7, SPARC, ITGA1 
# PC_ 5 
# Positive:  SLC12A2, FGF13, MACC1, HS6ST3, SEMA3C 
# Negative:  KCNE4, TPRG1, PTPRT, NOVA1, RP11-624L4.1 








# PC_ 1 
# Positive:  MECOM, FAM155A, MYRIP, ANO2, ST6GALNAC3 
# Negative:  CD44, SVEP1, LAMA2, EGFR, NOVA1 
# PC_ 2 
# Positive:  MT-RNR2, PDE3B, TMEM132C, SLC19A3, PCDH9 
# Negative:  FLT1, MCTP1, ADAMTS9, MIR222HG, KLF7 
# PC_ 3 
# Positive:  PDE3B, SORBS1, DIRC3, TMEM131L, LPL 
# Negative:  DCLK1, KAZN, NEGR1, LINC02511, BICC1 
# PC_ 4 
# Positive:  ATP13A3, ACSL4, MCL1, NAMPT, ZNF331 
# Negative:  IGF1, PTEN, PIK3R1, GAB1, HSPA1A 
# PC_ 5 
# Positive:  IL7R, PTPRC, LINC02147, CYTIP, SH2B3 
# Negative:  ADAMTS4, CSRP2, PNP, IL6, PKIG 


# Examine and visualize PCA results a few different ways
print(seurat_h5_reduced_CCA[["cca"]], dims = 1:5, nfeatures = 5)

# CC_ 1 
# Positive:  MCTP1, MECOM, FAM155A, EVA1C, ADGRL4 
# Negative:  SVEP1, PCDH9, ANK2, NOVA1, SOX5 
# CC_ 2 
# Positive:  CBLB, SAMD4A, MIR222HG, FBXL7, CDK17 
# Negative:  MT-RNR2, PDE3B, TENM3, SLC19A3, ACSS3 
# CC_ 3 
# Positive:  SORBS1, MAST4, ATP1B3, PDE3B, LPL 
# Negative:  KAZN, DCLK1, ROBO2, NEGR1, BICC1 
# CC_ 4 
# Positive:  ATP13A3, NAMPT, ACSL4, MCL1, ATP1B3 
# Negative:  PTEN, PIK3R1, IGF1, ADH1B, GAB1 
# CC_ 5 
# Positive:  ADAMTS4, PNP, IL6, RASGRF2, UGCG 
# Negative:  LINC02147, ID1, RANBP2, ESYT2, NOSTRIN 

## plots PCA
VizDimLoadings(seurat_h5_reduced_PCA, dims = 1:2, reduction = "pca")
DimPlot(seurat_h5_reduced_PCA, reduction = "pca")
DimHeatmap(seurat_h5_reduced_PCA, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(seurat_h5_reduced_PCA, dims = 1:15, cells = 500, balanced = TRUE)



## plots CCA
VizDimLoadings(seurat_h5_reduced_CCA, dims = 1:2, reduction = "cca")
DimPlot(seurat_h5_reduced_CCA, reduction = "cca")
# DimHeatmap(seurat_h5_reduced_CCA, dims = 1, cells = 500, balanced = TRUE)
# DimHeatmap(seurat_h5_reduced_CCA, dims = 1:15, cells = 500, balanced = TRUE)




# Determine the 'dimensionality' of the dataset
seurat_h5_reduced_dim_PCA = 
  JackStraw(seurat_h5_reduced_PCA, num.replicate = 100)
seurat_h5_reduced_dim_PCA
# An object of class Seurat 
# 33694 features across 6367 samples within 1 assay 
# Active assay: RNA (33694 features, 2000 variable features)
# 1 dimensional reduction calculated: pca




# An object of class Seurat 
# 38436 features across 538 samples within 1 assay 
# Active assay: RNA (38436 features, 2000 variable features)
# 1 dimensional reduction calculated: pca


ElbowPlot(seurat_h5_reduced_dim_PCA)






####### 7. Clustering -----------------------------------------------------
# Cluster the cells
seurat_h5_neighbors_PCA = 
  FindNeighbors(seurat_h5_reduced_dim_PCA, dims = 1:10)

seurat_h5_neighbors_PCA
# An object of class Seurat 
# 33694 features across 6367 samples within 1 assay 
# Active assay: RNA (33694 features, 2000 variable features)
# 1 dimensional reduction calculated: pca






# An object of class Seurat 
# 38436 features across 538 samples within 1 assay 
# Active assay: RNA (38436 features, 2000 variable features)
# 1 dimensional reduction calculated: pca

# # A Seurat-tibble abstraction: 538 × 10
# # Features=38436 | Cells=538 | Active assay=RNA | Assays=RNA
# .cell              orig.ident    nCount_RNA nFeature_RNA percent.mt    PC_1   PC_2   PC_3   PC_4     PC_5
# <chr>              <fct>              <dbl>        <int>      <dbl>   <dbl>  <dbl>  <dbl>  <dbl>    <dbl>
#   1 AGACAAACAGTTCCAA-1 SeuratProject       2923         1729      0.958  -5.92  -27.2    7.55   7.99 -19.0   
# 2 AGTTAGCTCCTACACC-1 SeuratProject       4080         2219      0.392 -26.3     7.76  23.8   13.9    5.08  
# 3 GTGTTCCTCTTGTTAC-1 SeuratProject       4240         2443      2.41  -35.7    -7.72  -9.34   7.69  -0.625 
# 4 GTCTGTCAGGGCAAGG-1 SeuratProject       3245         2256      4.68  -24.6   -12.8   -5.52   2.43  -0.0312
# 5 GCATGATGTTGGCCGT-1 SeuratProject       3508         2329      0.884 -26.7     6.32  18.4    7.21   2.10  
# 6 ACTTATCCAGTGTGCC-1 SeuratProject       3734         2442      0.107   1.71  -22.7    1.63  -6.45   7.32  
# 7 TCGTCCAAGCGCAATG-1 SeuratProject       3752         2409      3.78  -32.9   -12.1   -8.54   3.85  -0.0670
# 8 TTTCACAGTTTCGCTC-1 SeuratProject       3365         2271      0.267 -17.2    15.6   17.9  -17.8   -5.59  
# 9 TATCAGGTCTGGAAGG-1 SeuratProject       3792         2326      2.85  -36.9    -9.39 -12.6    2.27   0.239 
# 10 AGGATAAAGAGCAAGA-1 SeuratProject       3734         1911      0.348   0.928 -32.2    9.26   2.02 -16.0   
# # … with 528 more rows
# # ℹ Use `print(n = ...)` to see more rows



########   ## change resolution here 
seurat_h5_clusters = 
  FindClusters(seurat_h5_neighbors_PCA, 
               # resolution = 0.3
               ## change resolution here 
               resolution = 0.7
  )

# View(seurat_h5_clusters@meta.data)
colnames(seurat_h5_clusters@meta.data)

# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.7"
# [6] "seurat_clusters"



# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.7"
# [6] "seurat_clusters"


DimPlot(seurat_h5_clusters, 
        group.by = "RNA_snn_res.0.7",
        # group.by = "RNA_snn_res.0.3",
        label = TRUE)




##########################


seurat_h5_clusters_0.1 = 
  FindClusters(seurat_h5_neighbors_PCA, 
               resolution = 0.1)
colnames(seurat_h5_clusters_0.1@meta.data)
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.3" "seurat_clusters"
DimPlot(seurat_h5_clusters_0.1, 
        group.by = "RNA_snn_res.0.1",
        label = TRUE)



seurat_h5_clusters_0.2 = 
  FindClusters(seurat_h5_neighbors_PCA, 
               resolution = 0.2)
colnames(seurat_h5_clusters_0.2@meta.data)
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.3" "seurat_clusters"
DimPlot(seurat_h5_clusters_0.2, 
        group.by = "RNA_snn_res.0.2",
        label = TRUE)


seurat_h5_clusters_0.3 = 
  FindClusters(seurat_h5_neighbors_PCA, 
               resolution = 0.3)
colnames(seurat_h5_clusters_0.3@meta.data)
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.3" "seurat_clusters"
DimPlot(seurat_h5_clusters_0.3, 
        group.by = "RNA_snn_res.0.3",
        label = TRUE)



seurat_h5_clusters_0.4 = 
  FindClusters(seurat_h5_neighbors_PCA, 
               resolution = 0.4)
colnames(seurat_h5_clusters_0.4@meta.data)
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.4" "seurat_clusters"
DimPlot(seurat_h5_clusters_0.4, 
        group.by = "RNA_snn_res.0.4",
        label = TRUE)



seurat_h5_clusters_0.5 = 
  FindClusters(seurat_h5_neighbors_PCA, 
               resolution = 0.5)
colnames(seurat_h5_clusters_0.5@meta.data)
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.5" "seurat_clusters"
DimPlot(seurat_h5_clusters_0.5, 
        group.by = "RNA_snn_res.0.5",
        label = TRUE)



seurat_h5_clusters_0.6 = 
  FindClusters(seurat_h5_neighbors_PCA, 
               resolution = 0.6)
colnames(seurat_h5_clusters_0.6@meta.data)
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.6" "seurat_clusters"
DimPlot(seurat_h5_clusters_0.6, 
        group.by = "RNA_snn_res.0.6",
        label = TRUE)


###########  WORKS THE BEST
seurat_h5_clusters_0.7 = 
  FindClusters(seurat_h5_neighbors_PCA, 
               resolution = 0.7)
colnames(seurat_h5_clusters_0.7@meta.data)
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.7" "seurat_clusters"
DimPlot(seurat_h5_clusters_0.7, 
        group.by = "RNA_snn_res.0.7",
        label = TRUE)


seurat_h5_clusters_0.8 = 
  FindClusters(seurat_h5_neighbors_PCA, 
               resolution = 0.8)
colnames(seurat_h5_clusters_0.8@meta.data)
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.8" "seurat_clusters"
DimPlot(seurat_h5_clusters_0.8, 
        group.by = "RNA_snn_res.0.8",
        label = TRUE)

seurat_h5_clusters_0.9 = 
  FindClusters(seurat_h5_neighbors_PCA, 
               resolution = 0.9)
colnames(seurat_h5_clusters_0.9@meta.data)
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.8" "seurat_clusters"
DimPlot(seurat_h5_clusters_0.9, 
        group.by = "RNA_snn_res.0.9",
        label = TRUE)





##########################





# Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
# 
# Number of nodes: 538
# Number of edges: 14887
# 
# Running Louvain algorithm...
# 0%   10   20   30   40   50   60   70   80   90   100%
#   [----|----|----|----|----|----|----|----|----|----|
#      **************************************************|
#      Maximum modularity in 10 random starts: 0.8160
#    Number of communities: 5
#    Elapsed time: 0 seconds

seurat_h5_clusters
# An object of class Seurat 
# 33694 features across 6367 samples within 1 assay 
# Active assay: RNA (33694 features, 2000 variable features)
# 1 dimensional reduction calculated: pca



# An object of class Seurat 
# 38436 features across 538 samples within 1 assay 
# Active assay: RNA (38436 features, 2000 variable features)
# 1 dimensional reduction calculated: pca


# Look at cluster IDs of the first 5 cells
head(Idents(seurat_h5_clusters), 5)
# AAACCTGAGAACTGTA-1 AAACCTGAGAGTAAGG-1 AAACCTGAGCGCCTCA-1 AAACCTGCAAACGTGG-1 AAACCTGCAATGGAGC-1 
# 5                  5                  0                 11                  1 
# Levels: 0 1 2 3 4 5 6 7 8 9 10 11



# AGACAAACAGTTCCAA-1 AGTTAGCTCCTACACC-1 GTGTTCCTCTTGTTAC-1 GTCTGTCAGGGCAAGG-1 GCATGATGTTGGCCGT-1 
# 4                  1                  2                  2                  1 
# Levels: 0 1 2 3 4 5
Idents(seurat_h5_clusters)





# Look at cluster IDs of the first 5 cells
head(Idents(seurat_h5_clusters_0.7), 5)
# AAACCTGAGAACTGTA-1 AAACCTGAGAGTAAGG-1 AAACCTGAGCGCCTCA-1 AAACCTGCAAACGTGG-1 AAACCTGCAATGGAGC-1 
# 5                  5                  0                 11                  1 
# Levels: 0 1 2 3 4 5 6 7 8 9 10 11



# AGACAAACAGTTCCAA-1 AGTTAGCTCCTACACC-1 GTGTTCCTCTTGTTAC-1 GTCTGTCAGGGCAAGG-1 GCATGATGTTGGCCGT-1 
# 4                  1                  2                  2                  1 
# Levels: 0 1 2 3 4 5
Idents(seurat_h5_clusters_0.7)

# Levels: 0 1 2 3 4 
Idents(seurat_h5_clusters_0.5)
Idents(seurat_h5_clusters_0.6)
# Levels: 0 1 2 3 4 5 
Idents(seurat_h5_clusters_0.7)
Idents(seurat_h5_clusters_0.8)
# Levels: 0 1 2 3 4 5 6
Idents(seurat_h5_clusters_0.9)


seurat_h5_clusters_0.1
Idents(seurat_h5_clusters_0.1)
# Levels: 0 1
seurat_h5_clusters_0.2
Idents(seurat_h5_clusters_0.2)
# Levels: 0 1 2 3


# Run non-linear dimensional reduction (UMAP/tSNE)
# If you haven't installed UMAP, 
# you can do so via 
# reticulate::py_install(packages = 'umap-learn')
# seurat_h5_umap = 
#   RunUMAP(seurat_h5_clusters, dims = 1:10)
# seurat_h5_umap
# An object of class Seurat 
# 38436 features across 538 samples within 1 assay 
# Active assay: RNA (38436 features, 2000 variable features)
# 2 dimensional reductions calculated: pca, umap

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
## use "dev.off()" in case you have the following error 
## while running the DimPlot function

## Error in .Call.graphics(C_palette2, .Call(C_palette2, NULL)) : 
##   invalid graphics state
## REF: https://stackoverflow.com/questions/20155581/persistent-invalid-graphics-state-error-when-using-ggplot2/20627536#20627536?newreg=54712cb43d4746c795b1ed41a8f6bc16

# DimPlot(seurat_h5_umap, reduction = "umap", label = TRUE)


seurat_h5_umap_0.1 = 
  RunUMAP(seurat_h5_clusters_0.1, dims = 1:10)
seurat_h5_umap_0.1
umap_0.1 = 
  DimPlot(seurat_h5_umap_0.1, reduction = "umap", label = TRUE) +
  ggtitle("umap_0.1")
umap_0.1


seurat_h5_umap_0.2 = 
  RunUMAP(seurat_h5_clusters_0.2, dims = 1:10)
seurat_h5_umap_0.2
umap_0.2 = 
  DimPlot(seurat_h5_umap_0.2, reduction = "umap", label = TRUE) +
  ggtitle("umap_0.2")
umap_0.2







seurat_h5_umap_0.7 = 
  RunUMAP(seurat_h5_clusters_0.7, dims = 1:10)
seurat_h5_umap_0.7
umap_0.7 = 
  DimPlot(seurat_h5_umap_0.7, reduction = "umap", label = TRUE) +
  ggtitle("umap_0.7")
umap_0.7


seurat_h5_umap_0.3 = 
  RunUMAP(seurat_h5_clusters_0.3, dims = 1:10)
seurat_h5_umap_0.3
umap_0.3 = 
  DimPlot(seurat_h5_umap_0.3, reduction = "umap", label = TRUE) +
  ggtitle("umap_0.3")
umap_0.3


umap_0.7 + umap_0.3

DimPlot(seurat_h5_clusters_0.6, 
        group.by = "RNA_snn_res.0.6",
        label = TRUE)


seurat_h5_umap_0.6 = 
  RunUMAP(seurat_h5_clusters_0.6, dims = 1:10)


umap_0.6 = 
  DimPlot(seurat_h5_umap_0.6, 
          reduction = "umap", label = TRUE) +
  ggtitle("umap_0.6")

umap_0.6 + umap_0.7



####### 8. Finding DEGs -----------------------------------------------------
# Finding differentially expressed features (cluster biomarkers)

# find all markers of cluster 1
## potential    
cluster1.markers =  
  FindMarkers(seurat_h5_umap, 
              ident.1 = 1, 
              min.pct = 0.25)
head(cluster1.markers, n = 5)

DimPlot(seurat_h5_umap, 
        reduction = "umap", 
        group.by = "labs_new",
        # label = row.names(cluster1.markers)
        label = TRUE
)



# find all markers of cluster 2
## potential pre   
cluster2.markers =  
  FindMarkers(seurat_h5_umap, 
              ident.1 = 2, 
              min.pct = 0.25)
head(cluster2.markers, n = 5)



# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers = 
  FindMarkers(
    seurat_h5_umap, 
    ident.1 = 5, 
    ident.2 = c(0, 3), 
    min.pct = 0.25)
head(cluster5.markers, n = 5)

# # find markers for every cluster compared to all remaining cells, report only the positive
# # ones
# seurat_h5.markers = 
#   FindAllMarkers(
#     seurat_h5_umap, 
#     only.pos = TRUE, 
#     min.pct = 0.25, 
#     logfc.threshold = 0.25
#   )
# 
# seurat_h5.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 2, order_by = avg_log2FC)



# find markers for every cluster compared to all remaining cells, report only the positive
# ones
seurat_h5.markers_0.7 = 
  FindAllMarkers(
    seurat_h5_umap_0.7, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    logfc.threshold = 0.25
  )


seurat_h5.markers_0.7 %>%
  group_by(cluster)


seurat_h5.markers_0.7 %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)






seurat_h5.markers_0.3 = 
  FindAllMarkers(
    seurat_h5_umap_0.3, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    logfc.threshold = 0.25
  )

seurat_h5.markers_0.3 %>%
  group_by(cluster)


seurat_h5.markers_0.3 %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)




# cluster0.markers = 
#   FindMarkers(
#     seurat_h5_umap, 
#     ident.1 = 0, 
#     logfc.threshold = 0.25, 
#     test.use = "roc", 
#     only.pos = TRUE
#   )
# 
# VlnPlot(seurat_h5_umap, 
#         features = c("TTN", "NPPA", 
#                      "AP000251.3", "EGFL7", 
#                      "IFI27", "VWF"))
# 
# FeaturePlot(seurat_h5_umap, 
#             features = c("TTN", "NPPA"
#                          # , 
#                          # "AP000251.3", "EGFL7", 
#                          # "IFI27", "VWF"
#             ))
# 
# seurat_h5.markers %>%
#   group_by(cluster) %>%
#   top_n(n = 10, wt = avg_log2FC) -> top10
# DoHeatmap(seurat_h5_umap, features = top10$gene) + NoLegend()


# cluster_marker_info = 
#   seurat_h5.markers %>%
#   group_by(cluster) %>%
#   slice_max(n = 2, order_by = avg_log2FC) %>% 
#   select(cluster, gene)
# 
# cluster_marker_info
# 
# new.cluster.ids = c("cluster_0", "cluster_1", "cluster_2", 
#                     "cluster_3", "cluster_4", "cluster_5"
#                     # ,
#                     # "cluster_6", "cluster_7", 
#                     # "cluster_8", "cluster_9"
# )
# names(new.cluster.ids) = levels(seurat_h5_umap)
# new.cluster.ids
# 
# seurat_h5_cluster_id = RenameIdents(seurat_h5_umap, new.cluster.ids)
# 
# DimPlot(seurat_h5_cluster_id, 
#         reduction = "umap", label = TRUE,
#         pt.size = 0.5
# ) 
# # + NoLegend()


cluster_marker_info_0.7 = 
  seurat_h5.markers_0.7 %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) %>% 
  select(cluster, gene)

cluster_marker_info_0.7

new.cluster.ids_0.7 = c("cluster_0", "cluster_1", "cluster_2", 
                        "cluster_3", "cluster_4", "cluster_5")
names(new.cluster.ids_0.7) = levels(seurat_h5_umap_0.7)
new.cluster.ids_0.7

seurat_h5_cluster_id_0.7 = 
  RenameIdents(seurat_h5_umap_0.7, new.cluster.ids_0.7)

umap_cluster_0.7 = 
  DimPlot(seurat_h5_cluster_id_0.7, 
          reduction = "umap", label = TRUE) + 
  ggtitle("umap_cluster_0.7")

umap_cluster_0.7




cluster_marker_info_0.7 = 
  seurat_h5.markers_0.3 %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) %>% 
  select(cluster, gene)

cluster_marker_info_0.7

new.cluster.ids_0.7 = c("cluster_0", "cluster_1", "cluster_2", 
                        "cluster_3", "cluster_4", "cluster_5")
names(new.cluster.ids_0.7) = levels(seurat_h5_umap_0.7)
new.cluster.ids_0.7

seurat_h5_cluster_id_0.7 = 
  RenameIdents(seurat_h5_umap_0.7, new.cluster.ids_0.7)

umap_cluster_0.7 = 
  DimPlot(seurat_h5_cluster_id_0.7, 
          reduction = "umap", label = TRUE) + 
  ggtitle("umap_cluster_0.7")

umap_cluster_0.7





seurat_h5.markers_0.6 = 
  FindAllMarkers(
    seurat_h5_umap_0.6, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    logfc.threshold = 0.25
  )


cluster_marker_info_0.6 = 
  seurat_h5.markers_0.6 %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) %>% 
  select(cluster, gene)

cluster_marker_info_0.6

new.cluster.ids_0.6 = c("cluster_0", "cluster_1", "cluster_2", 
                        "cluster_3", "cluster_4")
names(new.cluster.ids_0.6) = levels(seurat_h5_umap_0.6)
new.cluster.ids_0.6

seurat_h5_cluster_id_0.6 = 
  RenameIdents(seurat_h5_umap_0.6, new.cluster.ids_0.6)

umap_cluster_0.6 = 
  DimPlot(seurat_h5_cluster_id_0.6, 
          reduction = "umap", label = TRUE) + 
  ggtitle("umap_cluster_0.6")

umap_cluster_0.6



umap_cluster_0.7 + umap_cluster_0.3
umap_cluster_0.7 + umap_cluster_0.6


new.cluster





#########################


seurat_h5.markers_0.2 = 
  FindAllMarkers(
    seurat_h5_umap_0.2, 
    only.pos = TRUE, 
    min.pct = 0.25, 
    logfc.threshold = 0.25
  )


cluster_marker_info_0.2 = 
  seurat_h5.markers_0.2 %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC) %>% 
  select(cluster, gene)

cluster_marker_info_0.2


# # A tibble: 8 × 2
# # Groups:   cluster [4]
# cluster gene   
# <fct>   <chr>  
#   1 0       CNTNAP2
# 2 0       RBFOX1 
# 3 1       PHLDB2 
# 4 1       LPL    
# 5 2       DCLK1  
# 6 2       NOVA1  
# 7 3       ADAMTS9
# 8 3       MCTP1  



new.cluster.ids_0.2 = c("cluster_0", "cluster_1", "cluster_2", 
                        "cluster_3")
names(new.cluster.ids_0.2) = levels(seurat_h5_umap_0.2)
new.cluster.ids_0.2
# 0           1           2           3 
# "cluster_0" "cluster_1" "cluster_2" "cluster_3"




seurat_h5_cluster_id_0.2 = 
  RenameIdents(seurat_h5_umap_0.2, new.cluster.ids_0.2)

umap_cluster_0.2 = 
  DimPlot(seurat_h5_cluster_id_0.2, 
          reduction = "umap", 
          # label = TRUE
          label = FALSE
  ) + 
  ggtitle("umap_cluster_0.2")

umap_cluster_0.2

















### SingleR
# SingleR() expects reference datasets to be normalized and log-transformed.
# pred.hesc <- SingleR(test = hESCs, ref = hpca.se, assay.type.test=1,
#                      labels = hpca.se$label.main)
# 
# pred.grun <- SingleR(test=sceG, ref=sceM, labels=sceM$label, de.method="wilcox")
# table(pred.grun$labels)

library(SingleR)
library(celldex)
# View(seurat_h5_reduced@meta.data)
ref = celldex::HumanPrimaryCellAtlasData()
ref
colData(ref)
# write_rds(seurat_h5_reduced, "../../output/seurat_h5_reduced.rds")
results = 
  SingleR(
    test = as.SingleCellExperiment(seurat_h5_cluster_id), 
    ref = ref,
    labels = ref$label.main
    # ref$label.fine
  )
results
# write_rds(results, "../../output/results.rds")
seurat_h5_cluster_id$singleR_lable = results$labels

###############
### adding the manual lables here 
## annotated according to 
## Supplemental Table 8 
## Single nucleus RNAseq of human    white    tissue RNAseq 
###############
## cells of interest
## Pre    (e.g., gene PNISR  maps to cluster 2 and 4 in seurat_h5.markers)
##     (e.g., gene RTN3  maps to cluster 1 in seurat_h5.markers)
meta_new = 
  seurat_h5_cluster_id@meta.data %>% 
  as_tibble() %>% 
  mutate(labs_t8 = 
           case_when(seurat_clusters == 1 ~ "   ", 
                     seurat_clusters == 2 ~ "Pre   ", 
                     # seurat_clusters == 4 ~ "Pre   ", # Endothelial_cells     
                     TRUE ~ as.character(seurat_clusters)
           )) 
# %>% 
# filter(labs_t8 == "Pre   ")



meta_new %>% 
  as_tibble() %>% 
  # filter(labs_t8 == "   ") %>% 
  filter(labs_t8 == "Pre   ") %>% 
  dplyr::count(singleR_lable, sort = TRUE)

labs_new = 
  meta_new %>% 
  mutate(labs_new = 
           case_when(
             singleR_lable == "Smooth_muscle_cells" ~ "SMC", 
             singleR_lable == "Endothelial_cells" ~ "EC", 
             # singleR_lable == "Fibroblasts" ~ "Fibroblasts", 
             singleR_lable == "Tissue_stem_cells" ~ "TSC",
             TRUE ~ "x"
             # TRUE ~ as.character(seurat_clusters)
           ))

# Smooth_muscle_cells
# Endothelial_cells
# Fibroblasts   
# Tissue_stem_cells


# meta_new = CreateSeuratObject(meta_new)

seurat_h5_cluster_id$labs_t8 = meta_new$labs_t8
seurat_h5_cluster_id$labs_new = labs_new$labs_new


seurat_h5_cluster_id[[]] %>% head(n=3) %>% colnames()
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.7" "seurat_clusters"
# [7] "singleR_lable"   "labels"          "labs_t8" 

# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.3" "seurat_clusters" "singleR_lable"  
# [8] "labs_t8"         "labs_new"



## works without group by
umap1 = DimPlot(
  seurat_h5_cluster_id,
  reduction = "umap", 
  ## group.by adds the lables
  # group.by = "singleR_lable",
  group.by = "labs_t8",
  # group.by = "seurat_clusters",
  label = TRUE)
umap1


umap2 = DimPlot(
  seurat_h5_cluster_id,
  reduction = "umap", 
  ## group.by adds the lables
  group.by = "singleR_lable",
  # group.by = "labs_t8",
  # group.by = "seurat_clusters",
  label = TRUE)
umap2

umap1 + umap2

umap3 = DimPlot(
  seurat_h5_cluster_id,
  reduction = "umap", 
  ## group.by adds the lables
  group.by = "labs_new",
  label = TRUE)
umap3

umap1 + umap3






DimPlot(seurat_h5_umap, 
        group.by = "seurat_clusters", 
        reduction = "umap", label = TRUE)



Pre    = 
  seurat_h5_cluster_id@meta.data %>% 
  filter(labs_t8 == "Pre   ") %>% 
  CreateSeuratObject()
DimPlot(seurat_h5_umap, group.by = "labs_t8", 
        reduction = "umap", label = TRUE)






## check 
FeaturePlot(seurat_h5_cluster_id, 
            features = c("Ptprc", "Cd3e"))
seurat_h5_cluster_id[[]] %>% 
  # as_tibble() %>% 
  # group_by(singleR_lable) %>% 
  dplyr::count(singleR_lable) %>% 
  arrange(desc(n)) 


table(seurat_h5_cluster_id@meta.data$singleR_lable) 
# %>% View()
# View(seurat_h5.markers)

###########################################################################


## marker genes of all clusters 
seurat_h5.markers

## annotated according to 
## Supplemental Table 8 
## Single nucleus RNAseq of human    white    tissue RNAseq 
###############
## cells of interest
## Pre    (e.g., gene PNISR  maps to cluster 2 and 4 in seurat_h5.markers)
##     (e.g., gene RTN3  maps to cluster 1 in seurat_h5.markers)

seurat_h5.markers_new 


seurat_h5.markers %>% 
  filter(cluster == 4) %>% 
  count(cluster)

# cluster   n
# 1       4 872

colnames(seurat_h5.markers)
# [1] "p_val"      "avg_log2FC" "pct.1"      "pct.2"      "p_val_adj"  "cluster"    "gene"   

cluster4.markers_gene = 
  cluster4.markers %>% 
  tibble::rownames_to_column(., var = "gene") 

colnames(cluster4.markers_gene)
# [1] "gene"       "myAUC"      "avg_diff"   "power"      "avg_log2FC" "pct.1"      "pct.2" 

cluster4.markers_gene %>% 
  distinct(gene, .keep_all = TRUE) %>% 
  dim()
# [1] 1628    7


# find markers for every cluster compared to all remaining cells, report only the positive
# ones


cluster1.markers = 
  FindMarkers(
    seurat_h5_umap, 
    ## point of change
    ident.1 = 1, 
    logfc.threshold = 0.25, 
    test.use = "roc", 
    only.pos = TRUE
  )

cluster1.markers

cluster2.markers = 
  FindMarkers(
    seurat_h5_umap,
    ## point of change
    ident.1 = 2, 
    logfc.threshold = 0.25, 
    test.use = "roc", 
    only.pos = TRUE
  )


cluster4.markers = 
  FindMarkers(
    seurat_h5_umap,
    ## point of change
    ident.1 = 4, 
    logfc.threshold = 0.25, 
    test.use = "roc", 
    only.pos = TRUE
  )


###########################################################################

###########################################################################
## REF: https://bioconductor.org/packages/devel/bioc/vignettes/SingleR/inst/doc/SingleR.html
# table(results$labels) %>% View()
# plotScoreHeatmap(results)
# plotDeltaDistribution(results, ncol = 3)
# summary(is.na(results$pruned.labels))
# # 
# # Mode   FALSE 
# # logical     538
# all.markers <- metadata(results)$de.genes
# 
# SingleCellExperiment(seurat_h5_cluster_id)$labels <- results$labels
# 
# # Beta cell-related markers
# library(scater)
# 
# #The function plotHeatmap comes from scater
# 
# # plotHeatmap(sceG, order_columns_by="labels",
# #             features=unique(unlist(all.markers$beta))) 
# 
# plotHeatmap(SingleCellExperiment(seurat_h5_cluster_id), 
#             order_columns_by="labels",
#             features=unique(unlist(all.markers$beta))) 






































