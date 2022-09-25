## https://stemangiola.github.io/tidyseurat/

# From CRAN
# install.packages("tidyseurat")

# From Github
# devtools::install_github("stemangiola/tidyseurat")


library(dplyr)
library(tidyr)
library(purrr)
library(magrittr)
library(ggplot2)
library(Seurat)
library(tidyseurat)


pbmc_small = seurat_h5 ## seurat object as the starting inputs
pbmc_small
pbmc_small@assays
# Use colourblind-friendly colours
friendly_cols <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC")

# Set theme
my_theme <-
  list(
    scale_fill_manual(values = friendly_cols),
    scale_color_manual(values = friendly_cols),
    theme_bw() +
      theme(
        panel.border = element_blank(),
        axis.line = element_line(),
        panel.grid.major = element_line(size = 0.2),
        panel.grid.minor = element_line(size = 0.1),
        text = element_text(size = 12),
        legend.position = "bottom",
        aspect.ratio = 1,
        strip.background = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10)),
        axis.title.y = element_text(margin = margin(t = 10, r = 10, b = 10, l = 10))
      )
  )


pbmc_small %>%
  tidyseurat::ggplot(aes(nFeature_RNA , fill = orig.ident
                         # , fill = groups
                         )) +
  geom_histogram() +
  my_theme


pbmc_small %>%
  tidyseurat::ggplot(aes(orig.ident, nCount_RNA , fill = orig.ident
                         # fill = groups
                         )) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  my_theme



pbmc_small %>%
  join_features(features = c("HLA-DRA", "LYZ")) %>%
  ggplot(aes(orig.ident, .abundance_RNA + 1, fill = orig.ident)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(size = nCount_RNA), alpha = 0.5, width = 0.2) +
  scale_y_log10() +
  my_theme


pbmc_small_pca <-
  pbmc_small %>%
  SCTransform(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  RunPCA(verbose = FALSE)

pbmc_small_pca



pbmc_small_cluster <-
  pbmc_small_pca %>%
  FindNeighbors(verbose = FALSE) %>%
  FindClusters(method = "igraph", verbose = FALSE)

pbmc_small_cluster

pbmc_small_cluster %>%
  tidyseurat::count(groups, seurat_clusters)




# Identify top 10 markers per cluster
markers <-
  pbmc_small_cluster %>%
  FindAllMarkers(only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25) %>%
  group_by(cluster) %>%
  top_n(10, avg_log2FC)

# Plot heatmap
pbmc_small_cluster %>%
  DoHeatmap(
    features = markers$gene,
    group.colors = friendly_cols
  )


pbmc_small_UMAP <-
  pbmc_small_cluster %>%
  RunUMAP(reduction = "pca", dims = 1:15, n.components = 3L)


pbmc_small_UMAP %>%
  plot_ly(
    x = ~`UMAP_1`,
    y = ~`UMAP_2`,
    z = ~`UMAP_3`,
    color = ~seurat_clusters,
    colors = friendly_cols[1:4]
  )



# Get cell type reference data
blueprint <- celldex::BlueprintEncodeData()

# Infer cell identities
cell_type_df <-
  GetAssayData(pbmc_small_UMAP, slot = 'counts', assay = "SCT") %>%
  log1p() %>%
  Matrix::Matrix(sparse = TRUE) %>%
  SingleR::SingleR(
    ref = blueprint,
    labels = blueprint$label.main,
    method = "single"
  ) %>%
  as.data.frame() %>%
  as_tibble(rownames = "cell") %>%
  select(cell, first.labels)

cell_type_df

# Join UMAP and cell type info
pbmc_small_cell_type <-
  pbmc_small_UMAP %>%
  left_join(cell_type_df, by = "cell")

# Reorder columns
pbmc_small_cell_type %>%
  tidyseurat::select(cell, first.labels, everything())


pbmc_small_cell_type %>%
  count(seurat_clusters, first.labels)

pbmc_small_cell_type %>%
  
  # Reshape and add classifier column
  pivot_longer(
    cols = c(seurat_clusters, first.labels),
    names_to = "classifier", values_to = "label"
  ) %>%
  
  # UMAP plots for cell type and cluster
  ggplot(aes(UMAP_1, UMAP_2, color = label)) +
  geom_point() +
  facet_wrap(~classifier) +
  my_theme


pbmc_small_cell_type %>%
  
  # Add some mitochondrial abundance values
  mutate(mitochondrial = rnorm(n())) %>%
  
  # Plot correlation
  join_features(features = c("CST3", "LYZ"), shape = "wide") %>%
  ggplot(aes(CST3 + 1, LYZ + 1, color = groups, size = mitochondrial)) +
  geom_point() +
  facet_wrap(~first.labels, scales = "free") +
  scale_x_log10() +
  scale_y_log10() +
  my_theme


pbmc_small_nested <-
  pbmc_small_cell_type %>%
  filter(first.labels != "Erythrocytes") %>%
  mutate(cell_class = if_else(`first.labels` %in% c("Macrophages", "Monocytes"), "myeloid", "lymphoid")) %>%
  nest(data = -cell_class)

pbmc_small_nested


pbmc_small_nested_reanalysed <-
  pbmc_small_nested %>%
  mutate(data = map(
    data, ~ .x %>%
      FindVariableFeatures(verbose = FALSE) %>%
      RunPCA(npcs = 10, verbose = FALSE) %>%
      FindNeighbors(verbose = FALSE) %>%
      FindClusters(method = "igraph", verbose = FALSE) %>%
      RunUMAP(reduction = "pca", dims = 1:10, n.components = 3L, verbose = FALSE)
  ))

pbmc_small_nested_reanalysed

pbmc_small_nested_reanalysed %>%
  
  # Convert to tibble otherwise Seurat drops reduced dimensions when unifying data sets.
  mutate(data = map(data, ~ .x %>% as_tibble())) %>%
  unnest(data) %>%
  
  # Define unique clusters
  unite("cluster", c(cell_class, seurat_clusters), remove = FALSE) %>%
  
  # Plotting
  ggplot(aes(UMAP_1, UMAP_2, color = cluster)) +
  geom_point() +
  facet_wrap(~cell_class) +
  my_theme








































