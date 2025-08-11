process pseudobulk_singlecell {
    tag "${single_cell_file}"
    label "process_high_memory"
    publishDir "${params.outdir}/expression_matrices/", mode: "copy"

    input:
    path single_cell_file
    path source_R

    output:
    path "ct_names.txt", emit: ct_names
    path "*pseudobulk.csv", emit: pseudobulk_counts
    path "gene_locations.csv", emit: gene_locations


    script:
    """
    #!/usr/bin/env Rscript
    library(Seurat)
    library(BPCells)
    library(dplyr)
    library(Matrix)
    library(data.table)

    source("$source_R")

    # Load object
    seuratobj <- readRDS("$single_cell_file")

    # (A) Per–cell type pseudobulk
    celltypelist <- Seurat::SplitObject(seuratobj, split.by = "${params.celltype_column}")

    aggregated_by_ct <- pseudobulk_counts(
      seuratlist = celltypelist,
      min.cells  = as.numeric(${params.min_cells}),
      indiv_col  = "${params.individual_column}",
      assay      = "${params.counts_assay}",
      slot       = "${params.counts_slot}"
    )

    for (i in seq_along(aggregated_by_ct)) {
      df <- aggregated_by_ct[[i]] %>% mutate(geneid = rownames(.))
      ct <- names(aggregated_by_ct)[i]
      data.table::fwrite(df, paste0(ct, "_pseudobulk.csv"))
    }

    # (B) Whole “Bulk” pseudobulk (all cells, still per individual)
    bulk_list <- list(Bulk = seuratobj)
    aggregated_bulk <- pseudobulk_counts(
      seuratlist = bulk_list,
      min.cells  = as.numeric(${params.min_cells}),
      indiv_col  = "${params.individual_column}",
      assay      = "${params.counts_assay}",
      slot       = "${params.counts_slot}"
    )
    bulk_df <- aggregated_bulk[[1]] %>% mutate(geneid = rownames(.))
    data.table::fwrite(bulk_df, "Bulk_pseudobulk.csv")

    # Gene locations (robust, based on assay rownames)
    counts_mat <- Seurat::GetAssayData(
      object = seuratobj,
      assay  = "${params.counts_assay}",
      slot   = "${params.counts_slot}"
    )
    gene_locations <- get_gene_locations(counts_mat)
    data.table::fwrite(gene_locations, "gene_locations.csv")

    # List of “cell types” emitted (now includes Bulk)
    ct_names <- c(names(aggregated_by_ct), "Bulk")
    writeLines(ct_names, "ct_names.txt")
    """
}