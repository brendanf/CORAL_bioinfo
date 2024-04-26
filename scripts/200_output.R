output_plan <- list(

  #### spike_summary ####
  # tibble:
  #  `seq_idx` integer: index of sequence in seq_all_trim
  #  `spike_id` character : name of matchign spike sequence
  #  `nsample` integer : number of samples the sequence was found in
  #  `nseqrun` integer : number of seqruns the sequence was found in
  #  `nread` integer: total number of reads for the sequence
  tar_fst_tbl(
    spike_summary,
    dplyr::left_join(spikes, seqtable_merged, by = "seq_idx") |>
      dplyr::rename(sample_key = sample, spike_id = cluster) |>
      dplyr::left_join(sample_table_key, by = "sample_key") |>
      dplyr::select(-sample_key) |>
      dplyr::summarize(
        nsample = dplyr::n_distinct(sample),
        nseqrun = dplyr::n_distinct(seqrun),
        nread = sum(nread),
        .by = c(seq_idx, spike_id)
      ),
    deployment = "main"
  ),

  #### write_spike_seqs ####
  # character: path + file name (fasta)
  #
  # ASV sequences which were identified as spikes
  tar_file_fast(
    write_spike_seqs,
    fastx_rename(
      infile = fastx_gz_extract(
        !!seq_all_trim,
        seq_index,
        spikes$seq_idx,
        withr::local_tempfile(fileext = ".fasta")
      ),
      names = glue::glue_data(
        spike_summary,
        "{seq_idx};{spike_id};nsample={nsample};nseqrun={nseqrun};nread={nread}"
      ),
      outfile = "output/spike_asvs.fasta"
    )
  ),

  #### write_asvtable ####
  # character: path + file name
  #
  # write the sparse ASV table to the output directory
  tar_file_fast(
    write_asvtable,
    write_and_return_file(asv_table, file.path("output", "asv_tab.rds"), type = "rds"),
    deployment = "main"
  ),

  tar_map(
    values = post_cluster_meta,
    names = .conf_level,

    ##### write_taxonomy_{.conf_level} #####
    # character : path and file name (.rds)
    #
    # write the ASV taxonomy to a file in the output directory
    tar_file_fast(
      write_taxonomy,
      tibble::column_to_rownames(taxon_table_ingroup, "seq_id") %>%
        write_and_return_file(sprintf("output/asv2tax_%s.rds", .conf_level), type = "rds")
    ),

    ##### write_duplicate_species_{.conf_level} #####
    # character : path and file name
    #
    # for testing purposes, write any species which exist in multiple places in
    # the taxonomy.  This file should be empty if everything has gone correctly.
    tar_file_fast(
      write_duplicate_species,
      dplyr::group_by(taxon_table_ingroup, !!TIP_RANK_VAR) %>%
        dplyr::filter(
          # !!TIP_RANK_VAR != "unk",
          dplyr::n_distinct(!!!rlang::syms(superranks(TIP_RANK))) > 1
        ) %>%
        dplyr::mutate(
          seq_idx = readr::parse_number(seq_id),
          classification = paste(!!!rlang::syms(superranks(TIP_RANK)), sep = ";") %>%
            ifelse(
              length(.) > 0L,
              sub(Biobase::lcPrefix(.), "", .),
              .
            ),
          name = sprintf("%s (%s) %s", !!TIP_RANK_VAR, classification, seq_id)
        ) %>%
        (
          \(x) {
            outfile <- sprintf("output/duplicates_%s.fasta", .conf_level)
            if (nrow(x) == 0) {
              if (file.exists(outfile)) unlink(outfile)
              character()
            } else {
              fastx_rename(
                infile = fastx_gz_extract(
                  infile = asv_seq,
                  index = asv_seq_index,
                  i = x$seq_idx,
                  outfile = withr::local_tempfile(fileext = ".fasta")
                ),
                names = write_and_return_file(
                  x$name,
                  withr::local_tempfile(fileext = ".txt")
                ),
                outfile = outfile
              )
            }
          }
        )()
    ),

    ##### write_otu_taxonomy_{.conf_level} #####
    # character : path and file name
    #
    # write the otu taxonomy to a file in the output directory
    tar_file_fast(
      write_otu_taxonomy,
      c(
        tibble::column_to_rownames(otu_taxonomy, "seq_id") %>%
          write_and_return_file(sprintf("output/otu_taxonomy_%s.rds", .conf_level), type = "rds"),
        dplyr::rename(otu_taxonomy, OTU = seq_id) %>%
          write_and_return_file(sprintf("output/otu_taxonomy_%s.tsv", .conf_level), type = "tsv")
      )
    ),
    if (do_dense_otu_table) {
      ##### write_otu_table_dense_{.conf_level} #####
      # character (length 2) : path and file name (.rds and .tsv)
      #
      # output the otu table in "dense" format, as required by most community
      # ecology analysis software
      tar_file_fast(
        write_otu_table_dense,
        otu_table_sparse %>%
          dplyr::mutate(sample = factor(sample, levels = unique(sample_table$sample))) %>%
          dplyr::summarize(nread = sum(nread), .by = c(sample, seq_id)) %>%
          tidyr::pivot_wider(names_from = seq_id, values_from = nread, values_fill = list(nread = 0L)) %>%
          tidyr::complete(sample) %>%
          dplyr::mutate(dplyr::across(where(is.integer), \(x) tidyr::replace_na(x, 0L))) %>%
          tibble::column_to_rownames("sample") %>%
          t() %>% {
            c(
              write_and_return_file(., sprintf("output/otu_table_%s.rds", .conf_level)),
              write_and_return_file(tibble::as_tibble(., rownames = "OTU"),
                                    sprintf("output/otu_table_%s.tsv", .conf_level),
                                    "tsv")
            )
          }
      )
    },

    ##### write_otu_refseq_{.conf_level} #####
    # character : path and file name (.fasta.gz)
    #
    # reference sequence for each OTU
    tar_file_fast(
      write_otu_refseq,
      fastx_rename(
        fastx_gz_random_access_extract(
          infile = asv_seq,
          index = asv_seq_index,
          i = readr::parse_number(otu_taxonomy$ref_seq_id),
          outfile = withr::local_tempfile(fileext = ".fasta")
        ),
        otu_taxonomy$seq_id,
        sprintf("output/otu_%s.fasta.gz", .conf_level)
      )
    ),

    ##### read_counts_{.conf_level} #####
    # tibble:
    #  `sample` character : sample name
    #  `seqrun` character : sequencing run name
    #  `raw_nread` integer : number of read pairs in input files
    #  `trim_nread` integer : number of read pairs remaining after adapter trimming
    #  `filt_nread` integer : number of read pairs remaining after quality filtering
    #  `denoise_nread` integer : number of merged reads remaining after denoising
    #  `uncross_nread` integer : number of merged reads remaining after removing tag-jumps
    #  `nochim1_nread` integer : number of merged reads remaining after de novo
    #    chimera removal
    #  `nochim2_nread` integer : number of merged reads remaining after reference
    #    based chimera removal
    #  `nospike_nread` integer : number of merged reads remaining after spike
    #    removal
    #  `full_length` integer : number of merged reads remaining after model scan for
    #    full-length amplicons
    #  `ingroup_nread` integer : number of merged reads remaining after outgroup
    #    removal
    if (nrow(orient_meta) > 1L) {
      tar_fst_tbl(
        read_counts,
        dplyr::bind_rows(
          (!!tar_map_bind_rows(seqrun_plan$dada2_meta_fwd)) |>
            dplyr::mutate(fastq_file = file.path(raw_path, fastq_R1)) |>
            dplyr::left_join(
              !!tar_map_bind_rows(seqrun_plan$raw_read_counts_fwd),
              by = "fastq_file"
            ) |>
            dplyr::left_join(
              !!tar_map_bind_rows(seqrun_plan$trim_read_counts_fwd),
              by = "trim_R1"
              ) |>
            dplyr::left_join(
              !!tar_map_bind_rows(seqrun_plan$filt_read_counts_fwd),
              by = "filt_R1"
            ) |>
            dplyr::mutate(sample_key = file_to_sample_key(filt_R1)),
          (!!tar_map_bind_rows(seqrun_plan$dada2_meta_rev)) |>
            dplyr::mutate(fastq_file = file.path(raw_path, fastq_R1)) |>
            dplyr::left_join(
              !!tar_map_bind_rows(seqrun_plan$raw_read_counts_rev),
              by = "fastq_file"
            ) |>
            dplyr::left_join(
              !!tar_map_bind_rows(seqrun_plan$trim_read_counts_rev),
              by = "trim_R1"
            ) |>
            dplyr::left_join(
              !!tar_map_bind_rows(seqrun_plan$filt_read_counts_rev),
              by = "filt_R1"
            ) |>
            dplyr::mutate(sample_key = file_to_sample_key(filt_R1))
        ) |>
          dplyr::summarize(
            raw_nread = max(raw_nread),
            dplyr::across(ends_with("nread") & !raw_nread, \(x) sum(x, na.rm = TRUE)),
            .by = c(sample, seqrun, sample_key)
          ) |>
          dplyr::left_join(
            !!tar_map_bind_rows(seqrun_plan$denoise_read_counts),
            by = "sample_key"
          ) |>
          dplyr::left_join(
            !!(if (isTRUE(do_uncross)) {
              tar_map_bind_rows(seqrun_plan$uncross_read_counts)
            } else {
              quote(tibble::tibble(sample_key = character()))
            }),
            by = "sample_key"
          ) |>
          dplyr::left_join(nochim1_read_counts, by = "sample_key") |>
          dplyr::left_join(
            nochim2_read_counts |>
              dplyr::summarize(dplyr::across(everything(), sum), .by = sample_key),
            by = "sample_key"
          ) |>
          dplyr::left_join(
            nospike_read_counts |>
              dplyr::summarize(dplyr::across(everything(), sum), .by = sample_key),
            by = "sample_key"
          ) |>
          dplyr::left_join(
            !!(if (isTRUE(do_model_filter)) {
              quote(
                full_length_read_counts |>
                  dplyr::summarize(dplyr::across(everything(), sum), .by = sample_key)
              )
            } else {
              quote(tibble::tibble(sample_key = "character"))
            }),
            by = "sample_key"
          ) |>
          dplyr::left_join(
            dplyr::group_by(otu_table_sparse, sample, seqrun) |>
              dplyr::summarize(ingroup_nread = sum(nread)),
            by = c("sample", "seqrun")
          ) |>
          dplyr::mutate(
            dplyr::across(
              dplyr::where(is.numeric),
              \(x) as.integer(tidyr::replace_na(x, 0L))
            )
          ) |>
          dplyr::select(sample, seqrun, raw_nread, trim_nread, filt_nread,
                        denoise_nread, any_of("uncross_nread"),
                        nochim1_nread, nochim2_nread, nospike_nread,
                        any_of("full_length_nread"), ingroup_nread)
      )
    } else {
      tar_fst_tbl(
        read_counts,
        (!!tar_map_bind_rows(seqrun_plan$dada2_meta)) |>
          dplyr::mutate(fastq_file = file.path(raw_path, fastq_R1)) |>
          dplyr::left_join(
            !!tar_map_bind_rows(seqrun_plan$raw_read_counts),
            by = "fastq_file"
          ) |>
          dplyr::left_join(
            !!tar_map_bind_rows(seqrun_plan$trim_read_counts),
            by = "trim_R1"
          ) |>
          dplyr::left_join(
            !!tar_map_bind_rows(seqrun_plan$filt_read_counts),
            by = "filt_R1"
          ) |>
          dplyr::mutate(sample_key = file_to_sample_key(filt_R1)) |>
          dplyr::left_join(
            !!tar_map_bind_rows(seqrun_plan$denoise_read_counts),
            by = "sample_key"
          ) |>
          dplyr::left_join(
            !!(if (isTRUE(do_uncross)) {
              tar_map_bind_rows(seqrun_plan$uncross_read_counts)
            } else {
              quote(tibble::tibble(sample_key = character()))
            }),
            by = "sample_key"
          ) |>
          dplyr::left_join(nochim1_read_counts, by = "sample_key") |>
          dplyr::left_join(
            nochim2_read_counts %>%
              dplyr::summarize(dplyr::across(everything(), sum), .by = sample_key),
            by = "sample_key"
          ) %>%
          dplyr::left_join(
            nospike_read_counts %>%
              dplyr::summarize(dplyr::across(everything(), sum), .by = sample_key),
            by = "sample_key"
          ) %>%
          dplyr::left_join(
            !!(if (isTRUE(do_model_filter)) {
              quote(
                full_length_read_counts |>
                  dplyr::summarize(dplyr::across(everything(), sum), .by = sample_key)
              )
            } else {
              quote(tibble::tibble(sample_key = "character"))
            }),
            by = "sample_key"
          ) %>%
          dplyr::left_join(
            dplyr::group_by(otu_table_sparse, sample, seqrun) %>%
              dplyr::summarize(ingroup_nread = sum(nread)),
            by = c("sample", "seqrun")
          ) %>%
          dplyr::mutate(
            dplyr::across(
              dplyr::where(is.numeric),
              \(x) as.integer(tidyr::replace_na(x, 0L))
            )
          ) |>
          dplyr::select(sample, seqrun, raw_nread, trim_nread, filt_nread,
                        denoise_nread,  any_of("uncross_nread"),
                        nochim1_nread, nochim2_nread, nospike_nread,
                        any_of("full_length_nread"), ingroup_nread)
      )
    },

    ##### write_read_counts_{.conf_level} #####
    # character : path and file name (.rds and .tsv)
    tar_file_fast(
      write_read_counts,
      c(
        write_and_return_file(
          read_counts,
          sprintf("output/read_counts_%s.rds", .conf_level),
          "rds"
        ),
        write_and_return_file(
          read_counts,
          sprintf("output/read_counts_%s.tsv", .conf_level),
          "tsv"
        )
      )
    ),

    ##### otu_abund_table_sparse #####
    tar_fst_tbl(
      otu_abund_table_sparse,
      otu_table_sparse |>
        dplyr::left_join(read_counts, by = "sample") |>
        dplyr::left_join(
          dplyr::select(sample_table, sample, spike_weight) |>
            unique(),
          by = "sample"
        ) |>
        dplyr::group_by(sample) |>
        dplyr::transmute(
          seq_id,
          nread,
          fread = nread/sum(nread),
          w = nread/(nochim2_nread - nospike_nread + 1) * spike_weight
        ) |>
        dplyr::ungroup()

    ),

    ##### write_otu_table_sparse_{.conf_level} #####
    # character : path and file name (.tsv)
    #
    # write the otu table as a sparse tsv
    tar_file_fast(
      write_otu_table_sparse,
      c(
        write_and_return_file(
          dplyr::rename(otu_abund_table_sparse, OTU = seq_id),
          sprintf("output/otu_table_sparse_%s.tsv", .conf_level),
          type = "tsv"
        ),
        write_and_return_file(
          dplyr::rename(otu_abund_table_sparse, OTU = seq_id),
          sprintf("output/otu_table_sparse_%s.rds", .conf_level),
          type = "rds"
        )
      )
    ),

    ##### otu_unknowns_{.conf_level} #####
    # tibble:
    #  `seq_id` character : unique OTU id
    #  {UNKNOWN_RANKS} factor : for each rank, is the OTU known, novel, or uncertain
    tar_fst_tbl(
      otu_unknowns,
      dplyr::inner_join(
        asv_unknown_prob,
        asv_otu_map,
        by = c("seq_id" = "ASV")
      ) |> dplyr::summarize(
        known_prob = max(known_prob),
        novel_prob = max(novel_prob),
        .by = c(OTU, rank)
      ) |>
        dplyr::mutate(
          status = dplyr::case_when(
            known_prob > .prob_threshold ~ "known",
            novel_prob > .prob_threshold ~ "novel",
            TRUE ~ "uncertain"
          ) |>
            factor(levels = c("novel", "uncertain", "known")),
          .keep = "unused"
        ) |>
        tidyr::pivot_wider(names_from = rank, values_from = status)
    ),

    ##### write_otu_unknowns_{.conf_level} #####
    # character: path and filename
    tar_file_fast(
      write_otu_unknowns,
      write_and_return_file.data.frame(
        otu_unknowns,
        sprintf("output/otu_unknowns_%s.tsv", .conf_level),
        type = "tsv"
      )
    )
  )
)

optimotu_plan <- c(optimotu_plan, output_plan)
