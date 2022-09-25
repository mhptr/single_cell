### SingleR
# SingleR() expects reference datasets to be normalized and log-transformed.


library(SingleR)
library(celldex)



ref = celldex::HumanPrimaryCellAtlasData()
ref
colData(ref)


results_0.7 = 
  SingleR(
    test = as.SingleCellExperiment(seurat_h5_cluster_id_0.7), 
    ref = ref,
    labels = ref$label.main
    # ref$label.fine
  )

results_0.7
# write_rds(results, "../../output/results.rds")
seurat_h5_cluster_id_0.7$singleR_lable = results_0.7$labels

###############
### adding the manual lables here 
## annotated according to 
## Supplemental Table 8 
## Single nucleus RNAseq of human subcutaneous white adipose tissue RNAseq 
###############
## cells of interest
## Preadipocytes (e.g., gene PNISR  maps to cluster 2 and 4 in seurat_h5.markers)
## Adipocytes (e.g., gene RTN3  maps to cluster 1 in seurat_h5.markers)
meta_new_0.7 = 
  seurat_h5_cluster_id_0.7@meta.data %>% 
  as_tibble() %>% 
  mutate(labs_t8 = 
           case_when(seurat_clusters == 1 ~ "Adipocytes", 
                     seurat_clusters == 2 ~ "Preadipocytes", 
                     # seurat_clusters == 4 ~ "Preadipocytes", # Endothelial_cells     
                     TRUE ~ as.character(seurat_clusters)
           )) 
# %>% 
# filter(labs_t8 == "Preadipocytes")



meta_new %>% 
  as_tibble() %>% 
  # filter(labs_t8 == "Adipocytes") %>% 
  filter(labs_t8 == "Preadipocytes") %>% 
  dplyr::count(singleR_lable, sort = TRUE)

labs_new_0.7 = 
  meta_new_0.7 %>% 
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

seurat_h5_cluster_id_0.7$labs_t8 = meta_new_0.7$labs_t8
seurat_h5_cluster_id_0.7$labs_new = labs_new_0.7$labs_new


seurat_h5_cluster_id_0.7[[]] %>% head(n=3) %>% colnames()
# [1] "orig.ident"      "nCount_RNA"      "nFeature_RNA"    "percent.mt"      "RNA_snn_res.0.7" "seurat_clusters"
# [7] "singleR_lable"   "labs_t8"         "labs_new" 



## works without group by
umap_labs_t8_0.7 = 
  DimPlot(
    seurat_h5_cluster_id_0.7,
    reduction = "umap", 
    group.by = "labs_t8",
    label = TRUE) +
  ggtitle("umap1_0.7 + labs_t8")
umap_labs_t8_0.7


umap_singleR_0.7 = 
  DimPlot(
    seurat_h5_cluster_id_0.7,
    reduction = "umap", 
    group.by = "singleR_lable",
    label = TRUE) +
  ggtitle("umap1_0.7 + singleR")
umap_singleR_0.7

umap_labs_t8_0.7 + umap_singleR_0.7

umap_labs_new_0.7 = 
  DimPlot(
    seurat_h5_cluster_id_0.7,
    reduction = "umap", 
    group.by = "labs_new",
    label = TRUE) +
  ggtitle("umap1_0.7 + labs_new")
umap_labs_new_0.7



umap_labs_t8_0.7 + umap_labs_new_0.7
umap_cluster_0.7 + umap_cluster_0.3







Preadipocytes = 
  seurat_h5_cluster_id@meta.data %>% 
  filter(labs_t8 == "Preadipocytes") %>% 
  CreateSeuratObject()
DimPlot(seurat_h5_umap, group.by = "labs_t8", 
        reduction = "umap", label = TRUE)



# find markers for every cluster compared to all remaining cells, report only the positive
# ones

## cluster 1 = tentative adipocyte
cluster1.markers_adipo =
  FindMarkers(
    seurat_h5_umap_0.7,
    ## point of change
    ident.1 = 1,
    logfc.threshold = 0.25,
    test.use = "roc",
    only.pos = TRUE
  )

# cluster1.markers_adipo
# 
# cluster1.markers_preadipo = 
#   FindMarkers(
#     seurat_h5_umap_0.7,
#     ## point of change
#     ident.1 = 2, 
#     logfc.threshold = 0.25, 
#     test.use = "roc", 
#     only.pos = TRUE
#   )
# cluster1.markers_preadipo


# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster1.markers_adipo = 
  FindMarkers(seurat_h5_umap_0.7, ident.1 = 1, 
              ident.2 = c(0,2,3,4,5), 
              min.pct = 0.25,
              only.pos = T)

cluster2.markers_preadipo = 
  FindMarkers(seurat_h5_umap_0.7, ident.1 = 2, 
              ident.2 = c(0,1,3,4,5), 
              min.pct = 0.25,
              only.pos = T)


cluster2.markers_preadipo


# write.csv(cluster1.markers_adipo, "../../output/cluster1.markers_adipo.csv")
# write.csv(cluster2.markers_preadipo, "../../output/cluster2.markers_preadipo.csv")

cluster_marker_info_0.7
# write.csv(cluster_marker_info_0.7, "../../output/cluster_marker_info_0.7.csv")



# new.cluster.ids_0.7 = 
#   c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono",
#                      "NK", "DC", "Platelet")
# names(new.cluster.ids) <- levels(pbmc)
# pbmc <- RenameIdents(pbmc, new.cluster.ids)
# DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
# 
# cluster_marker_info_0.7


#######################  12 Sep 2022 ####################

cluster0.markers = 
  FindMarkers(seurat_h5_umap_0.7, ident.1 = 0, 
              ident.2 = c(1,2,3,4,5), 
              min.pct = 0.25,
              only.pos = T)

cluster0.markers

cluster3.markers = 
  FindMarkers(seurat_h5_umap_0.7, ident.1 = 3, 
              ident.2 = c(0,1,2,4,5), 
              min.pct = 0.25,
              only.pos = T)

cluster3.markers


cluster4.markers = 
  FindMarkers(seurat_h5_umap_0.7, ident.1 = 4, 
              ident.2 = c(0,1,2,3,5), 
              min.pct = 0.25,
              only.pos = T)
cluster4.markers



cluster5.markers = 
  FindMarkers(seurat_h5_umap_0.7, ident.1 = 5, 
              ident.2 = c(0,1,2,3,4), 
              min.pct = 0.25,
              only.pos = T)

cluster5.markers



# write.csv(cluster0.markers, "../../output/cluster0.markers.csv")
# write.csv(cluster3.markers, "../../output/cluster3.markers.csv")
# write.csv(cluster4.markers, "../../output/cluster4.markers.csv")
# write.csv(cluster5.markers, "../../output/cluster5.markers.csv")



cluster0.markers_all =
  FindMarkers(
    seurat_h5_umap_0.7,
    ## point of change
    ident.1 = 0,
    logfc.threshold = 0.25,
    test.use = "roc",
    only.pos = TRUE
  )

cluster0.markers_all
# write.csv(cluster0.markers_all, "../../output/cluster0.markers_all.csv")



adipocyte = 
  FeaturePlot(
    seurat_h5_umap_0.7, 
    features = c("PHLDB2", "LPL")) 
adipocyte


preadipocyte = 
  FeaturePlot(
    seurat_h5_umap_0.7, 
    features = c("DCLK1", "NOVA1")) 
preadipocyte



umap_cluster_0.7


adipocyte + umap_labs_t8_0.7 + umap_cluster_0.7
preadipocyte + umap_labs_t8_0.7 + umap_cluster_0.7


hist_adi = 
  hist(
    cluster1.markers_adipo$p_val, 
    xlab = "Seurat analysis (Adipocytes): P-value",
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








