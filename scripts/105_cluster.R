# values for several variables to be used inside the tar_map
rank_meta <- tibble::tibble(
  .rank = tail(TAX_RANKS, -1), # e.g. phylum:species
  .rank_sym = rlang::syms(.rank),
  .parent_rank = head(TAX_RANKS, -1), # e.g. kingdom:genus
  .parent_rank_sym = rlang::syms(.parent_rank),
  .super_ranks = purrr::accumulate(.parent_rank, c),
  .parent_taxa = rlang::syms(paste0("taxon_table_", .parent_rank)), # for recursion
  .parent_pseudotaxa = rlang::syms(paste0("pseudotaxon_table_", .parent_rank)) # for recursion
)

taxon_table_TIP_RANK <- rlang::sym(sprintf("taxon_table_%s", TIP_RANK))
pseudotaxon_table_TIP_RANK <- rlang::sym(sprintf("pseudotaxon_table_%s", TIP_RANK))

#### rank_plan ####
# this ends up inside the reliablility_plan
rank_plan <- tar_map(
  values = rank_meta,
  names = .rank,

  ##### known_taxon_table_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV ID
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #
  # taxonomy as known before we start clustering at this rank
  tar_fst_tbl(
    known_taxon_table,
    asv_tax_prob_reads %>%
      dplyr::filter(
        rank == .rank,
        prob >= .prob_threshold,
        taxon != "unk"
      ) %>%
      dplyr::select(seq_id, .rank_sym := taxon) %>%
      dplyr::left_join(.parent_taxa, ., by = "seq_id") # results from previous rank
  ),

  ##### preclosed_taxon_table_{.rank}_{.conf_level} #####
  # grouped tibble:
  #  `seq_id` character : unique ASV ID
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #  `seq_idx` integer: index of the sequence in asv_taxsort_seq
  #  `tar_group` integer : grouping variable for dispatch to multiple dynamic
  #    jobs (matches values in .parent_rank)
  #
  # find only the groups which need to be closed-ref clustered
  # i.e. they have some known and some unknown
  tar_fst_tbl(
    preclosed_taxon_table,
    {
      out <- known_taxon_table |>
        dplyr::mutate(seq_idx_in = readr::parse_number(seq_id)) |>
        dplyr::left_join(asv_taxsort, by = "seq_idx_in") |>
        dplyr::arrange(dplyr::pick(all_of(.super_ranks)), seq_idx) |>
        dplyr::select(-seq_idx_in) |>
        dplyr::group_by(.parent_rank_sym) |>
        dplyr::filter(any(is.na(.rank_sym)) & !all(is.na(.rank_sym))) |>
        tar_group()
      # we can't dynamically map over an empty data frame
      # so give a single row.
      if (nrow(out) == 0) {
        known_taxon_table[1,] |>
          dplyr::mutate(seq_idx_in = readr::parse_number(seq_id)) |>
          dplyr::left_join(asv_taxsort, by = "seq_idx_in") |>
          dplyr::arrange(dplyr::pick(all_of(.super_ranks)), seq_idx) |>
          dplyr::select(-seq_idx_in) |>
          dplyr::group_by(.parent_rank_sym) |>
          tar_group()
      } else {
        out
      }
    },
    iteration = "group"
  ),

  ##### thresholds_{.rank}_{.conf_level} #####
  # named numeric : optimal clustering threshold at .rank for different taxa
  #   at the rank of .parent_rank (taxa are given by names)
  tar_target(
    thresholds,
    calc_taxon_thresholds(
      rank = .parent_rank,
      conf_level = "plausible",
      taxon_table = known_taxon_table,
      fmeasure_optima = fmeasure_optima
    )
  ),

  ##### clusters_closed_ref_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id of an "unknown" sequence
  #  `cluster` character : unique ASV id of a "known" sequence
  #
  # find matches within the chosen threshold between "unknown" and "known"
  # sequences
  tar_fst_tbl(
    clusters_closed_ref,
    {
      unknowns <- is.na(preclosed_taxon_table[[.rank]])
      taxon <- preclosed_taxon_table[[.parent_rank]][1]
      if (any(unknowns) && !all(unknowns)) {
        !!(if (protax_aligned) {
          quote(
            protax_besthit_closedref(
              infile = aligned_taxsort_seq,
              index = aligned_taxsort_seq_index,
              i = preclosed_taxon_table$seq_idx,
              unknowns = unknowns,
              thresh = thresholds[taxon],
              seq_width = amplicon_model_length
            )
          )
        } else {
          quote(
            vsearch_usearch_global_closed_ref(
              query =
                fastx_gz_random_access_extract(
                  asv_taxsort_seq,
                  asv_taxsort_seq_index,
                  preclosed_taxon_table$seq_idx[unknowns],
                  outfile = withr::local_tempfile(fileext=".fasta")
                ),
              ref =
                fastx_gz_random_access_extract(
                  asv_taxsort_seq,
                  asv_taxsort_seq_index,
                  preclosed_taxon_table$seq_idx[!unknowns],
                  outfile = withr::local_tempfile(fileext=".fasta")
                ),
              threshold = thresholds[taxon]/100
            )
          )
        })
      } else {
        tibble::tibble(seq_id = character(), cluster = character(), dist = double())
      }
    },
    pattern = map(preclosed_taxon_table) # per taxon at rank .parent_rank
  ),

  ##### closedref_taxon_table_{.rank}_{.conf_level} #####
  # grouped tibble:
  #  `seq_id` character : unique ASV ID
  #  {ROOT_RANK} character : taxonomic assignment at ROOT_RANK (e.g. kingdom)
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #
  # incorporates information from the closed-ref clustering into the
  # taxon table
  tar_fst_tbl(
    closedref_taxon_table,
    dplyr::left_join(
      known_taxon_table,
      clusters_closed_ref,
      by = "seq_id"
    ) %>%
      dplyr::left_join(
        dplyr::select(known_taxon_table, cluster = seq_id, cluster_taxon = .rank_sym),
        by = "cluster"
      ) %>%
      dplyr::mutate(
        .rank_sym := dplyr::coalesce(.rank_sym, cluster_taxon)
      ) %>%
      dplyr::select(-cluster, -cluster_taxon, -dist)
  ),


  ##### predenovo_taxon_table_{.rank}_{.conf_level} #####
  # grouped tibble:
  #  `seq_id` character : unique ASV ID
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #  `tar_group` integer : grouping variable for dispatch to multiple dynamic
  #    jobs (matches values in .parent_rank)
  #
  # find only the ASVs which need to be de novo clustered
  # i.e. they are still unknown
  tar_fst_tbl(
    predenovo_taxon_table,
    {
      out <- closedref_taxon_table |>
        dplyr::filter(is.na(.rank_sym)) |>
        dplyr::mutate(seq_idx_in = readr::parse_number(seq_id)) |>
        dplyr::left_join(asv_taxsort, by = "seq_idx_in") |>
        dplyr::arrange(dplyr::pick(all_of(.super_ranks)), seq_idx) |>
        dplyr::select(-seq_idx_in) |>
        dplyr::group_by(.parent_rank_sym) |>
        dplyr::filter(dplyr::n() > 1) |>
        tar_group()
      # we can't dynamically map over an empty data frame
      # so give a single row.
      if (nrow(out) == 0) {
        closedref_taxon_table[1,] |>
          dplyr::mutate(seq_idx_in = readr::parse_number(seq_id)) |>
          dplyr::left_join(asv_taxsort, by = "seq_idx_in") |>
          dplyr::arrange(dplyr::pick(all_of(.super_ranks)), seq_idx) |>
          dplyr::select(-seq_idx_in) |>
          dplyr::group_by(.parent_rank_sym) |>
          tar_group()
      } else {
        out
      }
    },
    iteration = "group"
  ),

  ##### denovo_thresholds_{.rank}_{.conf_level} #####
  # named list of named numeric
  #  list names: taxa at .parent_rank
  #  numeric values: optimal clustering thresholds
  #  numeric names: ranks from .rank to {TIP_RANK} (usually species)
  #
  # once we are doing de novo clustering at one rank, we will have to do it at
  # all subranks as well, and there is no new information that can change our
  # thresholds.  So it is easiest to just calculate the distance matrix once and
  # do all the clustering now.
  tar_target(
    denovo_thresholds,
    calc_subtaxon_thresholds(
      rank = .parent_rank,
      conf_level = "plausible",
      taxon_table = predenovo_taxon_table,
      fmeasure_optima = fmeasure_optima,
    ),
  ),

  ##### clusters_denovo_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : additional taxonomic assignments down to .parent_rank
  #  ... integer : unique cluster index, for ranks from .rank to TIP_RANK (usually species)
  tar_target(
    clusters_denovo,
    if (nrow(predenovo_taxon_table) > 1) {
      !!(if (protax_aligned) {
        quote(
          seq_cluster_protax(
            aln_seq = aligned_taxsort_seq,
            aln_index = aligned_taxsort_seq_index,
            which = predenovo_taxon_table$seq_idx,
            aln_len = amplicon_model_length,
            thresh = tryCatch(
              denovo_thresholds[[unique(predenovo_taxon_table[[.parent_rank]])]],
              error = function(e) denovo_thresholds[["_NA_"]]
            )
          ) |>
            t() |>
            tibble::as_tibble() %>%
            dplyr::bind_cols(
              dplyr::select(predenovo_taxon_table, -.rank_sym, -tar_group, -seq_idx),
              .
            )
        )
      } else {
        quote(
          optimotu::seq_cluster_usearch(
            seq = fastx_gz_extract(
              infile = asv_taxsort_seq,
              index = asv_taxsort_seq_index,
              i = predenovo_taxon_table$seq_idx,
              outfile = withr::local_tempfile(fileext = ".fasta")
            ),
            threshold_config = optimotu::threshold_set(
              tryCatch(
                denovo_thresholds[[unique(predenovo_taxon_table[[.parent_rank]])]],
                error = function(e) denovo_thresholds[["_NA_"]]
              )
            ),
            clust_config = optimotu::clust_tree(),
            parallel_config = optimotu::parallel_concurrent(2),
            usearch = "bin/usearch",
            usearch_ncpu = local_cpus()
          ) %>%
            t() %>%
            dplyr::as_tibble() %>%
            dplyr::bind_cols(
              dplyr::select(predenovo_taxon_table, -.rank_sym, -tar_group, -seq_idx),
              .
            )
        )
      })
    } else {
      c(
        c("seq_id", superranks(.rank)) %>%
          magrittr::set_names(., .) %>%
          purrr::map(~character(0)),
        c(.rank, subranks(.rank)) %>%
          magrittr::set_names(., .) %>%
          purrr::map(~integer(0))
      ) %>%
        tibble::as_tibble()
    },
    pattern = map(predenovo_taxon_table) # per taxon at .parent_rank
        ),

  ##### taxon_table_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV ID
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : taxon assigned at additional ranks (down to .rank)
  #
  # taxonomy for all named taxa down to .rank
  tar_fst_tbl(
    taxon_table,
    dplyr::filter(closedref_taxon_table, !is.na(.rank_sym))
  ),

  ##### pseudotaxon_table_{.rank}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #  ... character : additional taxonomic assignments down to .rank
  #  ... integer : unique cluster index, for ranks below .rank
  #
  # assign "pseudotaxon" names to de novo clusters at this .rank
  tar_fst_tbl(
    pseudotaxon_table,
    dplyr::bind_rows(
      clusters_denovo,
      .parent_pseudotaxa # results from previous rank
    ) %>%
      dplyr::arrange(seq_id) %>% # pseudotaxon numbers are ordered by ASV numbers
      dplyr::mutate(
        .rank_sym := paste(.parent_rank_sym, .rank_sym) %>%
          forcats::fct_inorder() %>%
          forcats::fct_relabel(
            ~names(name_seqs(., paste0("pseudo", .rank, "_")))
          ) %>%
          as.character()
      )
  )
)

#### reliability_plan ####
# repeats entire clustering process for different assignment probability
# thresholds

reliability_meta <- c(
  plausible = 0.5,
  reliable = 0.9
) %>%
  tibble::enframe(name = ".conf_level", value = ".prob_threshold")

reliability_plan <- tar_map(
  values = reliability_meta,
  names = .conf_level,

  ##### taxon_table_{ROOT_RANK}_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV ID
  #  {ROOT_RANK} character : taxon assigned at ROOT_RANK (e.g. kingdom)
  #
  # values for other ranks are calculated recursively
  # this should be everything, because PROTAX-fungi assigns all sequences
  # 100% probability of being fungi
  tar_target_raw(
    sprintf("taxon_table_%s", ROOT_RANK),
    quote(
      asv_tax_prob_reads %>%
      dplyr::filter(rank == ROOT_RANK) %>%
      dplyr::mutate(
        taxon = ifelse(prob < .prob_threshold, NA_character_, taxon)
      ) %>%
      dplyr::select(seq_id, {{ ROOT_RANK }} := taxon)
    ),
    format = "fst_tbl"
  ),

  ##### pseudotaxon_table_{ROOT_RANK}_{.conf_level} #####
  # NULL
  #
  # this is required because pseudotaxon_table_{SECOND_RANK} will try to access its
  # parent, but there are no pseudotaxa at the root level, so it is empty.
  # they are combined with dplyr::bind_row(), which will be fine if we give it
  # NULL instead of a 0-row tibble with the correct columns.
  tar_target_raw(
    sprintf("pseudotaxon_table_%s", ROOT_RANK),
    NULL
  ),

  rank_plan,

  ##### taxon_table_ingroup_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {ROOT_RANK} character : taxonomic assignment at ROOT_RANK (e.g. kingdom)
  #  ... character : taxonomic assignments at all intermediate levels
  #  {TIP_RANK} character : taxonomic assignment at TIP_RANK (e.g. species)
  #
  # combine the taxon table and pseudotaxon table, and select only phyla
  # which a) are named ingroup phyla; or b) contain more known ingroup (based on
  # reference) than known outgroup or ASVs with no known ROOT_RANK.
  tar_fst_tbl(
    taxon_table_ingroup,
    dplyr::bind_rows(
      !!(taxon_table_TIP_RANK),
      !!(pseudotaxon_table_TIP_RANK)
    ) %>%
      dplyr::mutate(
        known_outgroup = seq_id %in% asv_known_outgroup$seq_id,
        known_ingroup = seq_id %in% asv_known_ingroup$seq_id,
        unknown_outin = seq_id %in% asv_unknown_outin$seq_id
      ) %>%
      dplyr::group_by(dplyr::across({{SECOND_RANK_VAR}})) %>%
      dplyr::filter(
        !startsWith({{SECOND_RANK_VAR}}, "pseudo") |
          sum(known_ingroup) > sum(known_outgroup) + sum(unknown_outin)
      ) %>%
      dplyr::select(!where(is.logical)) %>%
      dplyr::arrange(seq_id)
  ),

  ##### asv_otu_map_{.conf_level} #####
  tar_fst_tbl(
    asv_otu_map,
    dplyr::semi_join(
      taxon_table_ingroup,
      asv_table,
      by = "seq_id"
    ) |>
      dplyr::left_join(
        dplyr::select(otu_taxonomy, OTU = seq_id, {{ TIP_RANK }}),
        by = TIP_RANK
      ) |>
      dplyr::select(ASV = seq_id, OTU)
  ),

  ##### otu_taxonomy_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique OTU ID
  #  `ref_seq_id` character : unique ASV ID of the most prevalent/abundant ASV
  #    in this OTU
  #  `nsample` integer : the number of samples in which this OTU occurs
  #  `nread` integer : the number of reads of this OTU across all samples
  #  {ROOT_RANK} character : taxonomic assignment at ROOT_RANK (e.g. kingdom)
  #  ... character : taxonomic assignments at all intermediate levels
  #  {TIP_RANK} character : taxonomic assignment at TIP_RANK (e.g. species)
  tar_fst_tbl(
    otu_taxonomy,
    asv_table %>%
      dplyr::group_by(seq_id) %>%
      dplyr::mutate(asv_nsample = dplyr::n(), asv_nread = sum(nread)) %>%
      dplyr::inner_join(taxon_table_ingroup, by = "seq_id") %>%
      dplyr::group_by(dplyr::across(all_of(TAX_RANKS))) %>%
      dplyr::arrange(dplyr::desc(asv_nsample), dplyr::desc(asv_nread)) %>%
      dplyr::summarize(
        nsample = as.integer(dplyr::n_distinct(sample)),
        nread = sum(nread),
        ref_seq_id = dplyr::first(seq_id)
      ) %>%
      dplyr::arrange(dplyr::desc(nsample), dplyr::desc(nread)) %>%
      name_seqs("OTU", "seq_id") %>%
      dplyr::select(seq_id, ref_seq_id, nsample, nread, everything())
  ),

  ##### otu_table_sparse_{.conf_level} #####
  # tibble:
  #  `seq_id` character : unique OTU id
  #  `sample` character : sample name
  #  `seqrun` character : name of sequencing run
  #  `nread` integer : number of reads
  #
  # OTU sample/abundance matrix, in sparse format (0's are not included)
  tar_fst_tbl(
    otu_table_sparse,
    asv_table %>%
      dplyr::inner_join(taxon_table_ingroup, by = "seq_id") |>
      dplyr::inner_join(asv_otu_map, by = c("seq_id" = "ASV")) |>
      dplyr::group_by(OTU, sample, seqrun) |>
      dplyr::summarise(nread = sum(nread), .groups = "drop") |>
      dplyr::select(seq_id = OTU, sample, seqrun, nread)
  )
)

#### clust_plan ####
# the outer plan which contains the reliability_plan

clust_plan <- list(

  ##### asv_known_outgroup #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {KNOWN_RANKS} character : taxonomic assignment at KNOWN_RANKS (e.g. kingdom)
  #
  # ASVs whose best match is to a species of known outgroup
  tar_fst_tbl(
    asv_known_outgroup,
    {
      out <- asv_best_hit_taxon
      outgroup_cols <- character(length(KNOWN_RANKS))
      for (i in seq_along(KNOWN_RANKS)) {
        outgroup_cols[i] <- paste0(KNOWN_RANKS[i], "_outgroup")
        out[[outgroup_cols[i]]] <-
          !is.na(out[[KNOWN_RANKS[i]]]) &
          !out[[KNOWN_RANKS[i]]] %in% c(KNOWN_TAXA[i], "unspecified", "Eukaryota_kgd_Incertae_sedis", "None")
      }
      dplyr::filter(out, dplyr::if_any(all_of(outgroup_cols))) |>
        dplyr::select(seq_id, dplyr::all_of(KNOWN_RANKS))
    }
  ),

  ##### asv_known_ingroup #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {KNOWN_RANKS} character : taxonomic assignment at ROOT_RANK (e.g. kingdom)
  #
  # ASVs whose best match is to an ingroup
  tar_fst_tbl(
    asv_known_ingroup,
    dplyr::filter(asv_best_hit_taxon, {{INGROUP_RANK_VAR}} == INGROUP_TAXON) |>
      dplyr::select(seq_id, dplyr::all_of(KNOWN_RANKS))
  ),

  ##### asv_unknown_outin #####
  # tibble:
  #  `seq_id` character : unique ASV id
  #  {ROOT_RANK} character : taxonomic assignment at ROOT_RANK (e.g. kingdom)
  #
  # ASVs whose best match is to a species whose identity at ROOT_RANK is unknown
  tar_target(
    asv_unknown_outin,
    asv_best_hit_taxon |>
      dplyr::anti_join(asv_known_outgroup, by = "seq_id") |>
      dplyr::anti_join(asv_known_ingroup, by = "seq_id") |>
      dplyr::select(seq_id, dplyr::all_of(KNOWN_RANKS))
  ),

  reliability_plan
)

optimotu_plan <- c(optimotu_plan, clust_plan)

post_cluster_meta <-
  purrr::map_dfc(reliability_plan, tar_select_names, everything()) |>
  dplyr::mutate_all(rlang::syms) |>
  dplyr::bind_cols(reliability_meta)
