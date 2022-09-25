
















### SingleR
# SingleR() expects reference datasets to be normalized and log-transformed.


library(SingleR)
library(celldex)



ref = celldex::HumanPrimaryCellAtlasData()
ref
colData(ref)


results_0.2 = 
  SingleR(
    test = as.SingleCellExperiment(seurat_h5_cluster_id_0.2), 
    ref = ref,
    labels = ref$label.main
    # ref$label.fine
  )

results_0.2
# write_rds(results, "../../output/results.rds")
seurat_h5_cluster_id_0.2$singleR_lable = results_0.2$labels

###############
### adding the manual lables here 
## annotated according to 
## Supplemental Table 8 
## Single nucleus RNAseq of human subcutaneous white adipose tissue RNAseq 
###############
## cells of interest
## Preadipocytes (e.g., gene PNISR  maps to cluster 2 and 4 in seurat_h5.markers)
## Adipocytes (e.g., gene RTN3  maps to cluster 1 in seurat_h5.markers)
meta_new_0.2 = 
  seurat_h5_cluster_id_0.2@meta.data %>% 
  as_tibble() %>% 
  mutate(labs_t8 = 
           case_when(seurat_clusters == 1 ~ "Adipocytes", 
                     seurat_clusters == 2 ~ "Preadipocytes", 
                     # seurat_clusters == 4 ~ "Preadipocytes", # Endothelial_cells     
                     TRUE ~ as.character(seurat_clusters)
           )) 

labs_t8_new2 = 
  meta_new_0.2 %>% 
  mutate(labs_t8_new2 = 
           case_when(
             seurat_clusters == 0 ~ "Fibroblast",
             seurat_clusters == 1 ~ "Adipocytes", 
             seurat_clusters == 2 ~ "Preadipocytes", 
             seurat_clusters == 3 ~ "Endothelial cells"
           ))

labs_t8_new2

####################
seurat_h5_cluster_id_0.2$labs_t8_new2 = labs_t8_new2$labs_t8_new2
seurat_h5_cluster_id_0.2[[]] %>% head(n=3) %>% colnames()
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"     
# [5] "RNA_snn_res.0.2" "seurat_clusters" "singleR_lable"   "labs_t8"        
# [9] "labs_new"        "labs_t8_new"     "labs_t8_new2"   



umap_cluster_0.2 = 
  DimPlot(seurat_h5_cluster_id_0.2, 
          reduction = "umap", 
          # label = TRUE
          label = FALSE
          ) + 
  ggtitle("umap_cluster_0.2")

umap_cluster_0.2




## works without group by

umap_labs_t8_new2_0.2 = 
  DimPlot(
    seurat_h5_cluster_id_0.2,
    reduction = "umap", 
    group.by = "labs_t8_new2",
    # label = TRUE
    label = FALSE
  ) +
  ggtitle("umap1_0.2 + labs_t8_new2")
umap_labs_t8_new2_0.2


umap_cluster_0.2 + umap_labs_t8_new2_0.2




# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster1.markers_0.2_adipo = 
  FindMarkers(seurat_h5_umap_0.2, ident.1 = 1, 
              ident.2 = c(0,2,3), 
              min.pct = 0.25,
              only.pos = T)

cluster1.markers_0.2_adipo




cluster2.markers_0.2_preadipo = 
  FindMarkers(seurat_h5_umap_0.2, ident.1 = 2, 
              ident.2 = c(0,1,3), 
              min.pct = 0.25,
              only.pos = T)

cluster2.markers_0.2_preadipo


# write.csv(cluster1.markers_0.2_adipo, "../../output/cluster1.markers_0.2_adipo.csv")
# write.csv(cluster2.markers_0.2_preadipo, "../../output/cluster2.markers_0.2_preadipo.csv")

cluster_marker_info_0.2
# write.csv(cluster_marker_info_0.2, "../../output/cluster_marker_info_0.2.csv")




#######################   ####################

adipocyte_0.2 = 
  FeaturePlot(
    seurat_h5_umap_0.2, 
    features = c("PHLDB2", "LPL")) 
adipocyte_0.2


preadipocyte_0.2 = 
  FeaturePlot(
    seurat_h5_umap_0.2, 
    features = c("DCLK1", "NOVA1")) 
preadipocyte_0.2



umap_cluster_0.2


adipocyte_0.2 + umap_labs_t8_new2_0.2 + umap_cluster_0.2
preadipocyte_0.2 + umap_labs_t8_new2_0.2 + umap_cluster_0.2









