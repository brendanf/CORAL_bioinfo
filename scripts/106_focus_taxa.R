if (length(target_taxa) > 0) {
  target_taxa_plan <- tar_map(
    values = post_cluster_meta,
    names = .conf_level,
    tar_fst_tbl(
      target_otus,
      find_target_taxa(
        target_taxa,
        asv_all_tax_prob,
        taxon_table_ingroup,
        otu_taxonomy
      )
    ),
    tar_file_fast(
      write_target_otus,
      write_and_return_file(
        target_otus,
        sprintf("output/target_taxon_otus_%s.rds", .conf_level),
        type = "rds"
      )
    )
  )
  optimotu_plan <- c(optimotu_plan, target_taxa_plan)
}
