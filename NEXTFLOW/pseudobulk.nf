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
    suppressPackageStartupMessages({
      library(Seurat)
      library(dplyr)
      library(Matrix)
      library(data.table)
    })

    source("$source_R")  # defines pseudobulk_counts(), get_gene_locations(), etc.

    # Load object
    seuratobj <- readRDS("$single_cell_file")

    # Per–cell type pseudobulk
    celltypelist <- Seurat::SplitObject(seuratobj, split.by = "${params.celltype_column}")

    aggregated_by_ct <- pseudobulk_counts(
      seuratlist = celltypelist,
      min.cells  = as.numeric(${params.min_cells}),
      indiv_col  = "${params.individual_column}",
      assay      = "${params.counts_assay}",
      slot       = "${params.counts_slot}"
    )

    # write only non-empty CT matrices
    valid_cts <- names(aggregated_by_ct)[vapply(aggregated_by_ct, ncol, integer(1)) > 0]
    for (ct in valid_cts) {
      df <- aggregated_by_ct[[ct]] %>% mutate(geneid = rownames(.))
      data.table::fwrite(df, paste0(ct, "_pseudobulk.csv"))
    }

    # Whole “Bulk” pseudobulk (all cells, still per individual)
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

    # Gene locations: use first non-empty CT, else Bulk (avoids undefined counts_mat)
    ref_mat <- if (length(valid_cts) > 0) aggregated_by_ct[[ valid_cts[1] ]] else aggregated_bulk[[1]]
    gene_locations <- get_gene_locations(ref_mat)
    data.table::fwrite(gene_locations, "gene_locations.csv")

    # ct_names = only non-empty CTs + Bulk
    ct_names <- c(valid_cts, "Bulk")
    writeLines(ct_names, "ct_names.txt")
    """
}