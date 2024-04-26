library(tarchetypes)

protax_plan <- list(

  #### protax_dir ####
  # character : directory name
  #
  # the main Protax directory (often a symlink). Here to be sure that it is
  # present and has not changed
  tar_file_fast(
    protax_dir,
    protax_root,
    deployment = "main"
  ),

  if (protax_aligned) {
    #### aligned protax ####
    ##### all_tax_prob #####
    # tibble:
    #  `seq_idx` integer : index of sequence in seq_all_trim
    #  `rank` ordered factor : rank of taxonomic assignment (phylum ... species)
    #  `parent_taxonomy` character : comma-separated taxonomy of parent to this taxon
    #  `taxon` character : name of the taxon
    #  `prob` numeric : probability that the asv in `seq_id` belongs to `taxon`
    #
    # In contrast to the unaligned case, each ASV may or may not have at least
    # one row at each rank; if no assignment at all was made at that rank,
    # then it will be missing.  If `taxon` is `NA`, this indicates an actual
    # prediction of "unknown taxon at this rank", and has an associated
    # `parent_taxon` and `prob`.
    # When alternative assignments are each above the probability threshold (10%)
    # then all are included on different rows.
    tar_target(
      all_tax_prob,
      fastx_split(
        asv_model_align,
        n = local_cpus(),
        outroot = tempfile(tmpdir = withr::local_tempdir())
      ) |>
        run_protax_animal(modeldir = protax_dir, id_is_int = TRUE, min_p = 0.02, info = TRUE, options = c("-m", "300")) |>
        dplyr::transmute(
          seq_idx,
          rank = int2rankfactor(rank),
          parent_taxonomy = paste(paste(KNOWN_TAXA, collapse = ","), taxonomy, sep = ",") |>
            sub(",[^,]+$", "", x = _),
          taxon = sub(".*,", "", taxonomy),
          prob,
          best_id,
          best_dist,
          second_id,
          second_dist
        )
        ,
      pattern = map(asv_model_align) # per seqbatch
    )
  } else {
    #### unaligned protax ####
    list(
      ##### protax_script #####
      # character: path and file name (executable)
      #
      # main protax script.  Slightly modified to accept various directories as
      # command line arguments
      tar_file_fast(
        protax_script,
        "scripts/runprotax"
      ),
      ##### protax #####
      # character of length 24 : path and filename for all protax output files
      tar_file_fast(
        protax,
        {
          protax_dir # dependency
          protax_script # dependency
          run_protax(
            seqs = fastx_gz_extract(
              infile = !!seq_all_trim,
              index = seq_index,
              i = seqbatch$seq_idx,
              outfile = withr::local_tempfile(fileext=".fasta"),
              hash = seqbatch_hash
            ),
            outdir = file.path(protax_path, tar_name()),
            modeldir = protax_model
          )
        },
        pattern = map(seqbatch, seqbatch_hash), # per seqbatch
        iteration = "list"
      ),

      ##### all_tax_prob #####
      # tibble:
      #  `seq_idx` integer : index of sequence in seq_all_trim
      #  `rank` ordered factor : rank of taxonomic assignment (phylum ... species)
      #  `parent_taxonomy` character : comma-separated taxonomy of parent to this taxon
      #  `taxon` character : name of the taxon
      #  `prob` numeric : probability that the asv in `seq_id` belongs to `taxon`
      #
      # Each ASV should have at least one row at each rank; if no assignment was
      # made at that rank, then `taxon` will be `NA`, `parent_taxon` may be `NA`,
      # and `prob` will be 0.
      # When alternative assignments are each above the probability threshold (10%)
      # then all are included on different rows.
      tar_fst_tbl(
        all_tax_prob,
        lapply(protax, grep, pattern = "query\\d.nameprob", value = TRUE) |>
          purrr::map_dfr(parse_protax_nameprob, id_is_int = TRUE)
      )
    )
  },

  #### asv_all_tax_prob ####
  # tibble:
  #  `seq_id` character : unique asv id
  #  `rank` ordered factor : rank of taxonomic assignment (phylum ... species)
  #  `parent_taxonomy` character : comma-separated taxonomy of parent to this taxon
  #  `taxon` character : name of the taxon
  #  `prob` numeric : probability that the asv in `seq_id` belongs to `taxon`
  tar_fst_tbl(
    asv_all_tax_prob,
    all_tax_prob |>
      dplyr::inner_join(asv_names, by = "seq_idx") |>
      dplyr::select(seq_id, everything() & !seq_idx),
    pattern = map(all_tax_prob),
    deployment = "main"
  ),

  #### asv_tax ####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {ROOT_RANK} character : taxonomic assignment at ROOT_RANK (e.g. kingdom)
  #  ... character: taxonomic assignments at intermediate ranks
  #  {TIP_RANK} character : taxonomic assignment at TIP_RANK (e.g. species)
  #
  # The most probable assignment for each ASV at each rank.  NA if there was no
  # assignment above Protax's reporting threshold
  tar_fst_tbl(
    asv_tax,
    asv_all_tax_prob %>%
      dplyr::group_by(rank, seq_id) %>%
      dplyr::summarize(taxon = dplyr::first(taxon), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = rank, values_from = taxon) %>%
      dplyr::bind_cols(as.list(magrittr::set_names(KNOWN_TAXA, KNOWN_RANKS))) %>%
      dplyr::select("seq_id", all_of(TAX_RANKS)),
    pattern = map(asv_all_tax_prob),
    deployment = "main"
  ),

  #### asv_tax_prob ####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {ROOT_RANK} numeric : probability for taxonomic assignment at ROOT_RANK
  #    (e.g., kingdom)
  #  ... numeric : probability for taxonomic assignments at intermediate ranks
  #  {TIP_RANK} numeric : probability for taxonomic assignment at TIP_RANK (e.g.,
  #    species)
  #
  # Associated probaility for the most probable assignment for each ASV at each
  # rank.  0 if there was no assignment above Protax's reporting threshold
  tar_fst_tbl(
    asv_tax_prob,
    asv_all_tax_prob %>%
      dplyr::group_by(rank, seq_id) %>%
      dplyr::summarize(prob = dplyr::first(prob), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = rank, values_from = prob) %>%
      dplyr::bind_cols(as.list(set_names(rep_len(1, length(KNOWN_RANKS)), KNOWN_RANKS))) %>%
      dplyr::select("seq_id", all_of(TAX_RANKS)),
    pattern = map(asv_all_tax_prob),
    deployment = "main"
  ),

  #### asv_unknown_prob ####
  tar_fst_tbl(
    asv_unknown_prob,
    asv_all_tax_prob |>
      dplyr::summarize(
        novel_prob = sum(prob[taxon == "unk"]),
        known_prob = max(prob[taxon != "unk"], 0),
        .by = c(seq_id, rank)
      ) |>
      tidyr::complete(seq_id, rank, fill = list(novel_prob = 0, known_prob = 0)) |>
      dplyr::filter(!rank %in% KNOWN_RANKS) |>
      dplyr::arrange(seq_id, desc(rank)) |>
      dplyr::mutate(
        novel_prob = cumsum(novel_prob),
        .by = seq_id
      ),
      pattern = map(asv_all_tax_prob),
      deployment = "main"
  ),

  #### asv_tax_prob_reads ####
  # tibble:
  #  `seq_id` character : unique asv ID
  #  `rank` character : taxonomic rank (e.g., kingdom...species)
  #  `taxon` character : name of taxon assigned at rank
  #  `prob` numeric : probability that taxon assignment is correct
  #  `nread` integer : number of reads for the ASV
  tar_fst_tbl(
    asv_tax_prob_reads,
    dplyr::full_join(
      tidyr::pivot_longer(asv_tax, all_of(TAX_RANKS), names_to = "rank", values_to = "taxon"),
      tidyr::pivot_longer(asv_tax_prob, all_of(TAX_RANKS), names_to = "rank", values_to = "prob"),
      by = c("seq_id", "rank")
    ) %>%
      dplyr::inner_join(asv_reads, by = "seq_id"),
    deployment = "main"
  )
)

optimotu_plan <- c(optimotu_plan, protax_plan)
