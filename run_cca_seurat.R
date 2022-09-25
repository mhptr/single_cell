###########################################################################
###########################################################################
###########################################################################

## reference : https://satijalab.org/seurat/reference/runcca

data("pbmc_small")
pbmc_small = seurat_h5
#> An object of class Seurat 
#> 230 features across 80 samples within 1 assay 
#> Active assay: RNA (230 features, 20 variable features)
#>  2 dimensional reductions calculated: pca, tsne# As CCA requires two datasets, we will split our test object into two just for this example
# pbmc1 <- subset(pbmc_small, cells = colnames(pbmc_small)[1:40])
pbmc1 <- subset(pbmc_small, cells = colnames(pbmc_small)[1:100])
# pbmc2 <- subset(pbmc_small, cells = colnames(x = pbmc_small)[41:80])
pbmc2 <- subset(pbmc_small, cells = colnames(x = pbmc_small)[101:200])
pbmc1[["group"]] <- "group1"
pbmc2[["group"]] <- "group2"
pbmc_cca <- RunCCA(object1 = pbmc1, object2 = pbmc2)
#> Warning: Fewer than 50 features used as input for CCA.#> Running CCA#> Warning: You're computing too large a percentage of total singular values, use a standard svd instead.#> Merging objects# Print results
print(x = pbmc_cca[["cca"]])
#> Warning: Requested number is larger than the number of available items (16). Setting to 16.#> CC_ 1 
#> Positive:  GNLY, PF4, TUBB1, SDPR, VDAC3, PPBP, PGRMC1, TREML1 
#> Negative:  HLA-DPB1, S100A9, S100A8, HLA-DQA1, RP11-290F20.3, CD1C, PARVB, RUFY1 #> Warning: Requested number is larger than the number of available items (16). Setting to 16.#> CC_ 2 
#> Positive:  SDPR, PF4, TUBB1, RUFY1, PPBP, TREML1, HLA-DPB1, PGRMC1 
#> Negative:  GNLY, S100A9, S100A8, RP11-290F20.3, VDAC3, HLA-DQA1, CD1C, PARVB #> Warning: Requested number is larger than the number of available items (16). Setting to 16.#> CC_ 3 
#> Positive:  S100A9, S100A8, PPBP, TUBB1, SDPR, TREML1, PF4, PGRMC1 
#> Negative:  GNLY, HLA-DQA1, HLA-DPB1, CD1C, VDAC3, RUFY1, RP11-290F20.3, PARVB #> Warning: Requested number is larger than the number of available items (16). Setting to 16.#> CC_ 4 
#> Positive:  RP11-290F20.3, HLA-DQA1, VDAC3, CD1C, HLA-DPB1, PARVB, RUFY1, S100A9 
#> Negative:  PGRMC1, SDPR, TUBB1, TREML1, PF4, GNLY, PPBP, S100A8 #> Warning: Requested number is larger than the number of available items (16). Setting to 16.#> CC_ 5 
#> Positive:  HLA-DQA1, VDAC3, CD1C, PF4, RUFY1, TUBB1, SDPR, PPBP 
#> Negative:  RP11-290F20.3, GNLY, HLA-DPB1, S100A8, PGRMC1, S100A9, TREML1, PARVB 
--