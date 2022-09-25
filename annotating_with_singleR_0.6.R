### SingleR
# SingleR() expects reference datasets to be normalized and log-transformed.


library(SingleR)
library(celldex)



ref = celldex::HumanPrimaryCellAtlasData()
ref
colData(ref)


results_0.6 = 
  SingleR(
    test = as.SingleCellExperiment(seurat_h5_cluster_id_0.6), 
    ref = ref,
    labels = ref$label.main
    # ref$label.fine
  )

results_0.6
# write_rds(results, "../../output/results.rds")
seurat_h5_cluster_id_0.6$singleR_lable = results_0.6$labels

###############
### adding the manual lables here 
## annotated according to 
## Supplemental Table 8 
## Single nucleus RNAseq of human subcutaneous white adipose tissue RNAseq 
###############
## cells of interest
## Preadipocytes (e.g., gene PNISR  maps to cluster 2 and 4 in seurat_h5.markers)
## Adipocytes (e.g., gene RTN3  maps to cluster 1 in seurat_h5.markers)
meta_new_0.6 = 
  seurat_h5_cluster_id_0.6@meta.data %>% 
  as_tibble() %>% 
  mutate(labs_t8 = 
           case_when(seurat_clusters == 1 ~ "Adipocytes", 
                     seurat_clusters == 2 ~ "Preadipocytes", 
                     # seurat_clusters == 4 ~ "Preadipocytes", # Endothelial_cells     
                     TRUE ~ as.character(seurat_clusters)
           )) 


labs_new_0.6 = 
  meta_new_0.6 %>% 
  mutate(labs_new = 
           case_when(
             singleR_lable == "Smooth_muscle_cells" ~ "SMC", 
             singleR_lable == "Endothelial_cells" ~ "EC", 
             # singleR_lable == "Fibroblasts" ~ "Fibroblasts", 
             singleR_lable == "Tissue_stem_cells" ~ "TSC",
             TRUE ~ "x"
             # TRUE ~ as.character(seurat_clusters)
           ))

labs_new_0.6

# Smooth_muscle_cells
# Endothelial_cells
# Fibroblasts   
# Tissue_stem_cells



labs_t8_new = 
  meta_new_0.6 %>% 
  mutate(labs_t8_new = 
           case_when(
             seurat_clusters == 0 ~ "Fibroblast cluster 1",
             seurat_clusters == 1 ~ "Adipocytes", 
             seurat_clusters == 2 ~ "Preadipocytes", 
             seurat_clusters == 3 ~ "Fibroblast cluster 2",
             seurat_clusters == 4 ~ "Endothelial cells"
           ))

labs_t8_new



labs_t8_new2 = 
  meta_new_0.6 %>% 
  mutate(labs_t8_new2 = 
           case_when(
             seurat_clusters == 0 ~ "Fibroblast",
             seurat_clusters == 1 ~ "Adipocytes", 
             seurat_clusters == 2 ~ "Preadipocytes", 
             seurat_clusters == 3 ~ "Fibroblast",
             seurat_clusters == 4 ~ "Endothelial cells"
           ))

labs_t8_new2








# meta_new = CreateSeuratObject(meta_new)

seurat_h5_cluster_id_0.6$labs_t8 = meta_new_0.6$labs_t8
seurat_h5_cluster_id_0.6$labs_new = labs_new_0.6$labs_new
seurat_h5_cluster_id_0.6$labs_t8_new = labs_t8_new$labs_t8_new
seurat_h5_cluster_id_0.6$labs_t8_new2 = labs_t8_new2$labs_t8_new2



seurat_h5_cluster_id_0.6[[]] %>% head(n=3) %>% colnames()
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"     
# [5] "RNA_snn_res.0.6" "seurat_clusters" "singleR_lable"   "labs_t8"        
# [9] "labs_new"        "labs_t8_new"     "labs_t8_new2"   



umap_cluster_0.6 = 
  DimPlot(seurat_h5_cluster_id_0.6, 
          reduction = "umap", 
          # label = TRUE
          label = FALSE
          ) + 
  ggtitle("umap_cluster_0.6")

umap_cluster_0.6




## works without group by
umap_labs_t8_0.6 = 
  DimPlot(
    seurat_h5_cluster_id_0.6,
    reduction = "umap", 
    group.by = "labs_t8",
    label = TRUE) +
  ggtitle("umap1_0.6 + labs_t8")
umap_labs_t8_0.6


umap_singleR_0.6 = 
  DimPlot(
    seurat_h5_cluster_id_0.6,
    reduction = "umap", 
    group.by = "singleR_lable",
    label = TRUE) +
  ggtitle("umap1_0.6 + singleR")
umap_singleR_0.6

umap_labs_t8_0.6 + umap_singleR_0.6

umap_labs_new_0.6 = 
  DimPlot(
    seurat_h5_cluster_id_0.6,
    reduction = "umap", 
    group.by = "labs_new",
    label = TRUE) +
  ggtitle("umap1_0.6 + labs_new")
umap_labs_new_0.6


umap_labs_t8_new_0.6 = 
  DimPlot(
    seurat_h5_cluster_id_0.6,
    reduction = "umap", 
    group.by = "labs_t8_new",
    # label = TRUE
    label = FALSE
    ) +
  ggtitle("umap1_0.6 + labs_t8_new")
umap_labs_t8_new_0.6



umap_labs_t8_new2_0.6 = 
  DimPlot(
    seurat_h5_cluster_id_0.6,
    reduction = "umap", 
    group.by = "labs_t8_new2",
    # label = TRUE
    label = FALSE
  ) +
  ggtitle("umap1_0.6 + labs_t8_new2")
umap_labs_t8_new2_0.6








umap_labs_t8_0.6 + umap_labs_new_0.6
umap_cluster_0.6 + umap_cluster_0.3
umap_cluster_0.6 + umap_labs_t8_new_0.6
umap_cluster_0.6 + umap_labs_t8_new2_0.6
umap_cluster_0.2 + umap_labs_t8_new2_0.6




# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster1.markers_0.6_adipo = 
  FindMarkers(seurat_h5_umap_0.6, ident.1 = 1, 
              ident.2 = c(0,2,3,4), 
              min.pct = 0.25,
              only.pos = T)

cluster1.markers_0.6_adipo




cluster2.markers_0.6_preadipo = 
  FindMarkers(seurat_h5_umap_0.6, ident.1 = 2, 
              ident.2 = c(0,1,3,4), 
              min.pct = 0.25,
              only.pos = T)

cluster2.markers_0.6_preadipo


# write.csv(cluster1.markers_0.6_adipo, "../../output/cluster1.markers_0.6_adipo.csv")
# write.csv(cluster2.markers_0.6_preadipo, "../../output/cluster2.markers_0.6_preadipo.csv")

cluster_marker_info_0.6
# write.csv(cluster_marker_info_0.6, "../../output/cluster_marker_info_0.6.csv")




#######################  12 Sep 2022 ####################

cluster0.markers_0.6 = 
  FindMarkers(seurat_h5_umap_0.6, ident.1 = 0, 
              ident.2 = c(1,2,3,4), 
              min.pct = 0.25,
              only.pos = T)

cluster0.markers_0.6

cluster3.markers_0.6 = 
  FindMarkers(seurat_h5_umap_0.6, ident.1 = 3, 
              ident.2 = c(0,1,2,4), 
              min.pct = 0.25,
              only.pos = T)

cluster3.markers_0.6


cluster4.markers_0.6 = 
  FindMarkers(seurat_h5_umap_0.6, ident.1 = 4, 
              ident.2 = c(0,1,2,3), 
              min.pct = 0.25,
              only.pos = T)
cluster4.markers_0.6



# write.csv(cluster0.markers_0.6, "../../output/cluster0.markers_0.6.csv")
# write.csv(cluster3.markers_0.6, "../../output/cluster3.markers_0.6.csv")
# write.csv(cluster4.markers_0.6, "../../output/cluster4.markers_0.6.csv")




adipocyte_0.6 = 
  FeaturePlot(
    seurat_h5_umap_0.6, 
    features = c("PHLDB2", "LPL")) 
adipocyte_0.6


preadipocyte_0.6 = 
  FeaturePlot(
    seurat_h5_umap_0.6, 
    features = c("DCLK1", "NOVA1")) 
preadipocyte_0.6



umap_cluster_0.6


adipocyte_0.6 + umap_labs_t8_0.6 + umap_cluster_0.6
preadipocyte_0.6 + umap_labs_t8_0.6 + umap_cluster_0.6


hist_adi = 
  hist(
    cluster1.markers_0.6_adipo$p_val, 
    xlab = "Seurat analysis (Adipocytes): P-value",
    labels = TRUE,
    col = "Blue", border = "black") #, breaks = seq(0,1,by=0.2)

hist(cluster2.markers_preadipo$p_val, 
     xlab = "Seurat analysis (Preadipocytes): P-value",
     col = "Blue", border = "black", )


hist(cluster1.markers_adipo$p_val, xlab = "Seurat analysis (Adipocytes): P-value",
     col = "Blue", border = "black", labels = TRUE)

hist(cluster1.markers_adipo$p_val, 
     xlab = "Seurat analysis (Adipocytes): P-value",
     col = "Blue", border = "black",
     # xlim=c(0,1), ylim=c(0,370), 
     # labels = TRUE
)


umap_cluster_0.3 = 
  DimPlot(seurat_h5_cluster_id_0.3, 
          reduction = "umap", label = TRUE) + 
  ggtitle("umap_cluster_0.3")




seurat_h5_cluster_id_0.4
hist(cluster1.markers_0.6_adipo$p_val)







