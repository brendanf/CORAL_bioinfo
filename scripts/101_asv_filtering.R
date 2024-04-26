## Filter unwanted sequences (chimeras, spikes) from DADA2 output
## and produce a final ASV table
## Brendan Furneaux

seq_trim_file <- "sequences/04_denoised/seq_all_trim.fasta.gz"
if (trim_options$action == "trim") {
  seq_all_trim <- quote(seq_all)
  seq_all_trim_file <- seq_all_file
} else {
  seq_all_trim <- quote(seq_trim)
  seq_all_trim_file <- seq_trim_file
}

asv_plan <- list(
  if (trim_options$action %in% c("retain", "lowercase", "none")) {
    #### seq_trim ####
    # character filename
    # all sequences in seq_all, after trimming primers
    # this may include duplicates
    tar_file_fast(
      seq_trim,
      cutadapt_filter_trim(
        seq_all,
        primer = trim_primer_merged,
        options = cutadapt_options(
          max_err = trim_options$max_err,
          min_overlap = trim_options$min_overlap,
          action = "trim",
          discard_untrimmed = FALSE
        ),
        ncpu = local_cpus(),
        trim = "sequences/04_denoised/seq_all_trim.fasta.gz"
      )
    )
  },

  #### seq_index ####
  # character filename
  # index file for fast access to sequences in seq_all_trim
  tar_file_fast(
    seq_index,
    fastx_gz_index(!!seq_all_trim),
    deployment = "main"
  ),

  #### seqbatch ####
  # grouped tibble:
  #  `seq_idx` integer: index of this sequence in seq_all_trim
  #  `tar_group` integer: which batch this sequence is assigned to
  #  `seq_id` character: within-batch index of this sequence
  #
  # Divide the unique ASV sequences from DADA2 into batches for further
  # processing.
  # In order to save time when the pipeline is re-run after new sequences are
  # added, cache and re-use the previous batch assignments.  Then sequences
  # which already existed do not need to be re-run
  #
  # seqbatches.fst: data.frame with same columns as seqbatch
  tar_fst_tbl(
    seqbatch,
    {
      batches_file <- "data/seqbatches.fst"
      new_batchkey <- tibble::tibble(seq_idx = seq_len(sequence_size(!!seq_all_trim)))
      if (
        file.exists(batches_file) &&
        nrow(old_batchkey <- fst::read_fst(batches_file)) <= nrow(new_batchkey)
      ) {
        old_batchkey <- fst::read_fst(batches_file)
        nbatch_old <- max(old_batchkey$tar_group)
        mean_old_batchsize <- nrow(old_batchkey) / nbatch_old
        new_batchkey <- dplyr::anti_join(new_batchkey, old_batchkey, by = "seq_idx")
        if (nrow(new_batchkey) > 0L) {
          # new_batches$seq_id <- seqhash(new_batches$seq)
          min_nbatch_new <- ceiling(nrow(new_batchkey) / max_batchsize)
          if (min_nbatch_new + nbatch_old < n_workers) {
            # if the cached number of batches is less than the target number of
            # workers, and we can fit all the new sequences in the target number
            # of workers, then do that. This is the normal case when adding new
            # seqruns to the analysis (if the seqruns are small enough that the
            # batchsize is less than the maximum)
            nbatch_new <- n_workers - nbatch_old
          } else {
            # otherwise add new batches of the same average size as the old ones
            nbatch_new <- ceiling(round(nrow(new_batchkey) / mean_old_batchsize))
          }
          new_batchkey$tar_group <- rep(
            seq_len(nbatch_new) + nbatch_old,
            each = ceiling(nrow(new_batchkey) / nbatch_new),
            length.out = nrow(new_batchkey)
          )
          new_batchkey <- dplyr::mutate(
            new_batchkey,
            seq_id = as.character(seq_len(dplyr::n())),
            .by = tar_group
          )
          new_batchkey <- rbind(old_batchkey, new_batchkey)
          fst::write_fst(new_batchkey, batches_file)
        } else {
          new_batchkey <- old_batchkey
        }
      }
      if (!"tar_group" %in% names(new_batchkey)) {
        nbatch_new <- ceiling(max(
          nrow(new_batchkey) / max_batchsize, # maximum size batches
          n_workers # one batch per worker
        ))
        new_batchkey$tar_group <- rep(
          seq_len(nbatch_new),
          each = ceiling(nrow(new_batchkey) / nbatch_new),
          length.out = nrow(new_batchkey)
        )
        new_batchkey <-
          dplyr::mutate(
            new_batchkey,
            seq_id = as.character(seq_len(dplyr::n())),
            .by = tar_group
          )
        fst::write_fst(new_batchkey, batches_file)
      }
      new_batchkey
    },
    iteration = "group"
  ),

  #### seqbatch_hash ####
  tar_target(
    seqbatch_hash,
    fastx_gz_hash(
      infile = !!seq_all_trim,
      index = seq_index,
      start = min(seqbatch$seq_idx),
      n = nrow(seqbatch)
    ),
    pattern = map(seqbatch)
  ),

  #### unaligned_ref_seqs ####
  # character filename
  # sequences to use as reference for uchime
  tar_file_fast(
    unaligned_ref_seqs,
    outgroup_reference_file
  ),

  #### ref_chimeras ####
  # `integer` vector: indices in `seq_all` which are reference-based chimeras
  #
  # Find reference-based chimeras in the current seqbatch.
  tar_target(
    ref_chimeras,
    vsearch_uchime_ref(
      query = fastx_gz_extract(
        infile = seq_all_trim_file, # actual file not a dependency
        index = seq_index,
        i = seqbatch$seq_idx,
        outfile = withr::local_tempfile(fileext=".fasta.gz"),
        hash = seqbatch_hash
      ),
      ref = unaligned_ref_seqs,
      ncpu = local_cpus(),
      id_only = TRUE,
      id_is_int = TRUE
    ),
    pattern = map(seqbatch, seqbatch_hash) # per seqbatch
  ),

  #### nochim2_read_counts ####
  # tibble:
  #  `sample_key` character: as `sample_table$sample_key`
  #  `nochim2_nread` integer: number of sequences in the sample after second
  #    chimera filtering
  tar_target(
    nochim2_read_counts,
    dplyr::filter(
      seqtable_merged,
      !seq_idx %in% denovo_chimeras,
      !seq_idx %in% ref_chimeras
    ) |>
      dplyr::summarize(nochim2_nread = sum(nread), .by = sample) |>
      dplyr::rename(sample_key = sample)
  ),

  #### spikes ####
  # tibble:
  #  `seq_idx` integer: index of sequence in seq_all_trim
  #  `cluster` character: name of matching spike sequence
  #
  # find spike sequences in the current seqbatch
  tar_fst_tbl(
    spikes,
    vsearch_usearch_global(
      fastx_gz_extract(
        infile = seq_all_trim_file, # actual file not a dependency
        index = seq_index,
        i = seqbatch$seq_idx,
        outfile = withr::local_tempfile(fileext=".fasta.gz"),
        hash = seqbatch_hash
      ),
      "protaxFungi/addedmodel/amptk_synmock.udb",
      global = FALSE,
      threshold = 0.9,
      id_is_int = TRUE
    ),
    pattern = map(seqbatch, ref_chimeras) # per seqbatch
  ),

  #### nospike_read_counts ####
  # tibble:
  #  `sample_key` character: as `sample_table$sample_key`
  #  `nospike_nread` integer: number of sequences in the sample after spike
  #    removal.
  tar_fst_tbl(
    nospike_read_counts,
    dplyr::filter(
      seqtable_merged,
      !seq_idx %in% denovo_chimeras,
      !seq_idx %in% ref_chimeras
    ) |>
      dplyr::anti_join(spikes, by = "seq_idx") |>
      dplyr::summarize(nospike_nread = sum(nread), .by = sample) |>
      dplyr::rename(sample_key = sample)
  ),

  #### Amplicon models ####
  if (!identical(amplicon_model_type, "none")) {
    list(
      ##### amplicon_model_file #####
      # `character`: file name of CM or HMM file
      tar_file_fast(
        amplicon_model_file,
        model_file
      ),
      if (identical(amplicon_model_type, "CM")) {
        ##### CM #####
        list(
          ###### amplicon_model_length ######
          tar_target(
            amplicon_model_length,
            readLines(amplicon_model_file, n = 100) |>
              purrr::keep(startsWith, "CLEN") |>
              readr::parse_number()
          ),

          if (do_model_filter_only) {
            ###### do_model_filter_only ######
            ####### amplicon_model_match #######
            # `tibble`:
            #  `seq_idx` character: index of the sequence in seq_all_trim
            #  `model_from` integer: first matching position in the CM
            #  `model_to` integer: last matching postion in the CM
            #  `bit_score` numeric: bit score of the match; higher is better.
            tar_fst_tbl(
              amplicon_model_match,
              {
                sfile <- tempfile(fileext = ".dat")
                inferrnal::cmalign(
                  amplicon_model_file,
                  fastx_gz_extract(
                    infile = seq_all_trim_file, # actual file not a dependency
                    index = seq_index,
                    i = seqbatch$seq_idx,
                    outfile = withr::local_tempfile(fileext=".fasta"),
                    hash = seqbatch_hash
                  ),
                  global = TRUE,
                  notrunc = TRUE,
                  cpu = local_cpus(),
                  sfile = sfile
                )
                read_sfile(sfile) |>
                  dplyr::transmute(
                    seq_idx = as.integer(seq_id),
                    model_from = cm_from,
                    model_to = cm_to,
                    bit_score = bit_sc
                  )
              },
              pattern = map(seqbatch, seqbatch_hash)
            )
          } else if (do_model_align_only) {
            ###### do_model_align_only ######
            ####### asv_model_align #######
            tar_file_fast(
              asv_model_align,
              inferrnal::cmalign(
                amplicon_model_file,
                fastx_gz_extract(
                  infile = seq_all_trim_file, # actual file not a dependency
                  index = seq_index,
                  i = seqbatch$seq_idx,
                  outfile = withr::local_tempfile(fileext=".fasta"),
                  hash = seqbatch_hash
                ),
                global = TRUE,
                notrunc = TRUE,
                dnaout = TRUE,
                cpu = local_cpus()
              ) |>
                consensus_columns() |>
                write_sequence(
                  fname = sprintf(
                    "sequences/05_aligned/batch%05i.fasta.gz",
                    seqbatch$tar_group[1]
                  ),
                  compress = TRUE
                ),
              pattern = map(seqbatch, seqbatch_hash)
            )
          } else if (do_model_both) {
            ###### do_model_both ######
            list(
              ####### asv_cm_align #######
              # `character`: two file names, for the alignment and the alignment
              # stats
              tar_file_fast(
                asv_cm_align,
                {
                  sfile <- sprintf(
                    "sequences/05_aligned/batch%05i.sfile",
                    seqbatch$tar_group[1]
                  )
                  ensure_directory(sfile)
                  inferrnal::cmalign(
                    amplicon_model_file,
                    fastx_gz_extract(
                      infile = seq_all_trim_file, # actual file not a dependency
                      index = seq_index,
                      i = seqbatch$seq_idx,
                      outfile = withr::local_tempfile(fileext=".fasta"),
                      hash = seqbatch_hash
                    ),
                    global = TRUE,
                    notrunc = TRUE,
                    dnaout = TRUE,
                    cpu = local_cpus(),
                    sfile = sfile
                  ) |>
                    consensus_columns() |>
                    write_sequence(
                      fname = sprintf(
                        "sequences/05_aligned/batch%05i.fasta.gz",
                        seqbatch$tar_group[1]
                      ),
                      compress = TRUE
                    )|>
                    c(sfile)
                },
                pattern = map(seqbatch, seqbatch_hash)
              ),

              ####### asv_model_align #######
              tar_file_fast(
                asv_model_align,
                asv_cm_align[1],
                pattern = map(asv_cm_align),
                deployment = "main"
              ),

              ####### amplicon_model_match #######
              # `tibble`:
              #  `seq_idx` character: index of the sequence in seq_all_trim
              #  `model_from` integer: first matching position in the CM
              #  `model_to` integer: last matching postion in the CM
              #  `bit_score` numeric: bit score of the match; higher is better.
              tar_fst_tbl(
                amplicon_model_match,
                read_sfile(asv_cm_align[2]) |>
                  dplyr::transmute(
                    seq_idx = as.integer(seq_id),
                    model_from = cm_from,
                    model_to = cm_to,
                    bit_score = bit_sc
                  ),
                pattern = map(asv_cm_align),
                deployment = "main"
              )
            )
          }
        )
      } else if (identical(amplicon_model_type, "HMM")) {
        ##### HMM #####
        list(
          ###### amplicon_model_length ######
          tar_target(
            amplicon_model_length,
            readLines(amplicon_model_file, n = 100) |>
              purrr::keep(startsWith, "LENG") |>
              readr::parse_number()
          ),

          if (do_model_filter) {
            ###### amplicon_model_match ######
            # `tibble`:
            #  `seq_idx` (integer) : index of the sequence in seq_all_trim
            #  `hmm_from` (integer) : start position of the match in the HMM
            #  `hmm_to` (integer) : end position of the match in the HMM
            #  `bit_score` (numeric) : score for the match (higher is better)
            tar_fst_tbl(
              amplicon_model_match,
              nhmmer(
                seqs = fastx_gz_extract(
                  infile = seq_all_trim_file,
                  index = seq_index,
                  i = seqbatch$seq_idx,
                  outfile = withr::local_tempfile(fileext=".fasta"),
                  hash = seqbatch_hash
                ),
                hmm = amplicon_model_file
              ) |>
                dplyr::transmute(
                  seq_idx = as.integer(seq_name),
                  model_from = hmm_from,
                  model_to = hmm_to,
                  bit_score
                ),
              pattern = map(seqbatch, seqbatch_hash)
            )
          },

          if (do_model_align) {
            ###### asv_model_align ######
            tar_file_fast(
              asv_model_align,
              fastx_gz_extract(
                infile = seq_all_trim_file,
                index = seq_index,
                i = seqbatch$seq_idx,
                outfile = withr::local_tempfile(fileext=".fasta"),
                hash = seqbatch_hash
              ) |>
                fastx_split(
                  n = local_cpus(),
                  outroot = tempfile(tmpdir = withr::local_tempdir()),
                  compress = TRUE
                ) |>
                hmmalign(
                  hmm = amplicon_model_file,
                  outfile = sprintf("sequences/05_aligned/batch%05i.fasta.gz", seqbatch$tar_group[1])
                ),
              pattern = map(seqbatch, seqbatch_hash)
            )
          },

          if (do_numt_filter) {
            ###### numts ######
            tar_fst_tbl(
              numts,
              detect_numts(asv_model_align, id_is_int = TRUE),
              pattern = map(asv_model_align),
              deployment = "main"
            )
          }
        )
      } else {
        stop("invalid value for amplicon_model_type: ", amplicon_model_type)
      },

      if (do_model_filter) {
        list(
          ##### asv_full_length #####
          # `integer`: index of sequences in seqs_dedup which are full-length CM matches
          tar_target(
            asv_full_length,
            dplyr::filter(
              amplicon_model_match,
              bit_score >= model_filter$min_model_score,
              model_from <= model_filter$max_model_start,
              model_to >= model_filter$min_model_end
            )$seq_idx,
            pattern = map(amplicon_model_match)
          ),

          ##### full_length_read_counts #####
          tar_fst_tbl(
            full_length_read_counts,
            dplyr::filter(
              seqtable_merged,
              !seq_idx %in% denovo_chimeras,
              !seq_idx %in% ref_chimeras,
              seq_idx %in% asv_full_length
            ) |>
              dplyr::anti_join(spikes, by = "seq_idx") |>
              dplyr::summarize(full_length_nread = sum(nread), .by = sample) |>
              dplyr::rename(sample_key = sample)
          )
        )
      }
    )
  },

  if (protax_aligned) {
    #### aligned ####
    list(
      ##### unaligned_ref_index #####
      if (endsWith(outgroup_reference_file, ".gz")) {
        tar_file_fast(
          unaligned_ref_index,
          fastx_gz_index(unaligned_ref_seqs)
        )
      } else {
        tar_fst_tbl(
          unaligned_ref_index,
          Biostrings::fasta.index(unaligned_ref_seqs)
        )
      },

      ##### outgroup_seqbatch #####
      # tibble:
      #  `batch` integer : index of batch
      #  `batch_id` character : name of batch
      #  `from` integer : index in outgroup reference file to start this batch
      #  `to` integer : index in outgroup reference file to end this batch
      tar_fst_tbl(
        outgroup_seqbatch,
        {
          n_seq <- sequence_size(unaligned_ref_seqs)
          n_batch <- as.integer(ceiling(n_seq / max_batchsize))
          batchsize <- as.integer(ceiling(n_seq / n_batch))
          tibble::tibble(
            batch = seq_len(n_batch),
            batch_id = make_seq_names(n_batch, "ref"),
            from = (batch - 1L) * batchsize + 1L,
            to = dplyr::lead(from, 1L, default = n_seq + 1L) - 1L
          )
        }
      ),

      ##### outgroup_aligned #####
      # character: path and file name for fasta file representing a batch of
      #   aligned reference sequences
      tar_file_fast(
        outgroup_aligned,
        {
          tf <- withr::local_tempfile(fileext = ".fasta")
          !!if (endsWith(outgroup_reference_file, ".gz")) {
            quote(
              fastx_gz_extract(
                infile = unaligned_ref_seqs,
                index = unaligned_ref_index,
                i = seq(outgroup_seqbatch$from, outgroup_seqbatch$to),
                outfile = tf
              )
            )
          } else {
            quote(
              system(sprintf(
                "tail -c+%d %s | head -n%d >%s",
                unaligned_ref_index$offset[outgroup_seqbatch$from] + 1,
                unaligned_ref_seqs,
                2*(outgroup_seqbatch$to - outgroup_seqbatch$from + 1),
                tf
              ))
            )
          }
          fastx_split(
            tf,
            n = local_cpus(),
            outroot = tempfile(tmpdir = withr::local_tempdir()),
            compress = FALSE
          ) |>
            hmmalign(
              hmm = amplicon_model_file,
              outfile = sprintf("sequences/05_aligned/%s.fasta.gz", outgroup_seqbatch$batch_id)
            )
        },
        pattern = map(outgroup_seqbatch)
      ),
      ##### outgroup_taxonomy #####
      tar_fst_tbl(
        outgroup_taxonomy,
        names(Biostrings::fasta.seqlengths(outgroup_reference_file)) |>
          tibble::tibble(name = _) |>
          tidyr::separate(name, c("ref_id", "bin", "country", "taxonomy"), sep = "[|]") |>
          tidyr::separate(taxonomy, TAX_RANKS, sep = ",", extra = "drop")
      ),
      ##### best_hit_taxon #####
      # tibble:
      #  `seq_idx` integer: index in seq_all_trim
      #  `ref_id` character: reference sequence id of best hit
      #  `dist` numeric: distance to best hit
      #  `sh_id` character: species hypothesis of best hit
      #  {INGROUP_RANK} character: taxon of best hit at {INGROUP_RANK} (e.g., kingdom)
      tar_fst_tbl(
        best_hit_taxon,
        run_protax_besthit(
          aln_ref = outgroup_aligned,
          aln_query =
            fastx_split(
              asv_model_align,
              n = local_cpus(),
              outroot = tempfile(tmpdir = withr::local_tempdir()),
              compress = TRUE
            ),
          options = c(
            "-m", "300",
            "-l", amplicon_model_length,
            "-r", dplyr::last(outgroup_seqbatch$to)
          ),
          query_id_is_int = TRUE,
          ref_id_is_int = FALSE
        ) |>
          dplyr::left_join(
            dplyr::select(outgroup_taxonomy, ref_id, all_of(KNOWN_RANKS)),
            by = "ref_id"
          ),
        pattern = map(asv_model_align)
      )
    )
  } else {
    #### unaligned ####
    list(
      ##### best_hit_udb #####
      # character: path and file name for udb of reference sequences
      #
      # build a udb index for fast vsearch
      tar_file_fast(
        best_hit_udb,
        build_filtered_udb(
          infile = outgroup_reference_file,
          outfile = "sequences/outgroup_reference.udb",
          blacklist = c(
            "SH1154235.09FU", # chimeric; partial matches to two different fungi but labeled as a fern
            "SH1240531.09FU" # chimera of two fungi, labeled as a plant
          ),
          usearch = Sys.which("vsearch")
        )
      ),

      ##### best_hit_taxon #####
      # tibble:
      #  `seq_idx` integer: index in seq_all_trim
      #  `ref_id` character: reference sequence id of best hit
      #  `dist` numeric: distance to best hit
      #  `sh_id` character: species hypothesis of best hit
      #  {KNOWN_RANKS} character: taxon of best hit at {KNOWN_RANKS} (e.g., kingdom)
      if (is.null(outgroup_taxonomy_file)) {
        tar_fst_tbl(
          best_hit_taxon,
          vsearch_usearch_global(
            query = fastx_gz_extract(
              infile = seq_all_trim_file, # actual file not a dependency
              index = seq_index,
              i = seqbatch$seq_idx,
              outfile = withr::local_tempfile(fileext=".fasta"),
              hash = seqbatch_hash
            ),
            ref = best_hit_udb,
            threshold = 0.8,
            global = FALSE,
            id_is_int = TRUE
          ) |>
            dplyr::arrange(seq_idx) |>
            tidyr::separate(cluster, c("ref_id", "sh_id", "taxonomy"), sep = "[|]") |>
            tidyr::separate(taxonomy, TAX_RANKS, sep = ",", fill = "right")

        )
      } else {
        tar_fst_tbl(
          best_hit_taxon,
          vsearch_usearch_global(
            query = fastx_gz_extract(
              infile = seq_all_trim_file, # actual file not a dependency
              index = seq_index,
              i = seqbatch$seq_idx,
              outfile = withr::local_tempfile(fileext=".fasta"),
              hash = seqbatch_hash
            ),
            ref = best_hit_udb,
            threshold = 0.8,
            global = FALSE,
            id_is_int = TRUE
          ) |>
            dplyr::arrange(seq_idx) |>
            tidyr::separate(cluster, c("ref_id", "sh_id"), sep = "_") |>
            dplyr::left_join(
              readr::read_tsv(
                outgroup_taxonomy_file,
                col_names = c("sh_id", "taxonomy"),
                col_types = "cc-------"
              ),
              by = "sh_id"
            ) |>
            dplyr::mutate(
              {{INGROUP_RANK_VAR}} := sub(";.*", "", taxonomy) |> substr(4, 100),
              .keep = "unused"
            ),
          pattern = map(seqbatch, seqbatch_hash) # per seqbatch
        )
      }
    )
  },

  #### spike_table ####
  # tibble:
  #  `sample` character: sample name (as in sample_table$sample)
  #  `seqrun` character: sequencing run (as in sample_table$seqrun)
  #  `seq_id` character: unique spike ASV id, in format "Spike[0-9]+". numbers
  #    are 0-padded
  #  'seq_idx` integer: index of sequence in seqs_dedup
  #  `spike_id` character: name of best hit spike sequence
  #  `nread` integer: number of reads
  #
  # global table of ASVs which are predicted to be spikes
  # (this does include chimeras)
  tar_fst_tbl(
    spike_table,
    dplyr::arrange(spikes, seq_idx) |>
      name_seqs("Spike", "seq_id") |>
      dplyr::left_join(seqtable_merged, by = "seq_idx") |>
      dplyr::rename(spike_id = cluster, sample_key = sample) |>
      dplyr::left_join(sample_table_key, by = "sample_key") |>
      dplyr::select(sample, seqrun, seq_id, seq_idx, spike_id, nread)
  ),

  #### asv_names ####
  # tibble:
  #  `seq_idx` integer : index of a sequence in seqs_dedup
  #  `seq_id` character : unique ASV id, in format "ASV[0-9]+". numbers are
  #    0-padded
  tar_fst_tbl(
    asv_names,
    tibble::tibble(
      seq_idx =
        (!!(
          if (do_model_filter)
            quote(unname(asv_full_length))
          else
            quote(seq_len(sequence_size(!!seq_all_trim)))
        )) |>
        setdiff(denovo_chimeras) |>
        setdiff(ref_chimeras) |>
        setdiff(spikes$seq_idx) |>
        setdiff(!!(
          if (do_numt_filter)
            quote(numts$seq_idx)
          else
            quote(integer())
        )) |>
        sort()
    ) |>
      name_seqs("ASV", "seq_id")
  ),

  #### asv_table ####
  # tibble:
  #  `sample` character: sample name (as in sample_table$sample)
  #  `seqrun` character: sequencing run (as in sample_table$seqrun)
  #  `seq_id` character: unique ASV id, in format "ASV[0-9]+". numbers are
  #    0-padded
  #  `seq_idx` character: index of ASV sequence in seqs_dedup
  #  `nread` integer: number of reads
  tar_fst_tbl(
    asv_table,
    seqtable_merged |>
      dplyr::filter(
        !seq_idx %in% denovo_chimeras,
        !seq_idx %in% ref_chimeras,
        !seq_idx %in% spikes$seq_idx,
        !!(if (do_model_filter) quote(seq_idx %in% asv_full_length) else TRUE)
      ) |>
      dplyr::left_join(asv_names, by = "seq_idx") |>
      dplyr::rename(sample_key = sample) |>
      dplyr::left_join(sample_table_key, by = "sample_key") |>
      dplyr::select(sample, seqrun, seq_id, seq_idx, nread),
    deployment = "main"
  ),

  #### asv_reads ####
  # tibble:
  #  `seq_id` character: unique ASV id
  #  `nread` integer: total reads across all samples
  #
  # calculate total read counts for all ASVs (at least those present in asv_tax)
  tar_fst_tbl(
    asv_reads,
    asv_table %>%
      dplyr::group_by(seq_id) %>%
      dplyr::summarize(nread = sum(nread)) %>%
      dplyr::semi_join(asv_tax, by = "seq_id"),
    deployment = "main"
  ),

  #### asv_seq ####
  # `character` filename
  # sequence for each ASV
  tar_file_fast(
    asv_seq,
    write_sequence(
      Biostrings::readDNAStringSet(!!seq_all_trim)[asv_names$seq_idx] |>
        name_seqs(prefix = "ASV"),
      "sequences/04_denoised/asv.fasta.gz",
      compress = TRUE
    ),
    deployment = "main"
  ),

  #### asv_seq_index ####
  # `character` filename
  # sequence for each ASV
  tar_file_fast(
    asv_seq_index,
    fastx_gz_index(asv_seq),
    deployment = "main"
  ),

  #### asv_best_hit_taxon ####
  # `tibble`:
  #  `seq_id` character: unique ASV identifier
  #  `ref_id` character: reference sequence id of best hit
  #  `sh_id` character: species hypothesis of best hit
  #  `dist` numeric: distance to best hit
  #  {{INGROUP_RANK}} character: taxon at rank INGROUP_RANK (e.g. kingdom) of
  #    best hit
  tar_fst_tbl(
    asv_best_hit_taxon,
    dplyr::left_join(
      asv_names,
      best_hit_taxon,
      by = "seq_idx"
    ) |>
      dplyr::select(-seq_idx),
    deployment = "main"
  ),

  #### asv_taxsort ####
  # `tibble`:
  #  `seq_idx` integer: index of sequence in asv_taxsort_seq
  #  `seq_idx_in` integer: index of sequence in asv_seq
  tar_fst_tbl(
    asv_taxsort,
    tibble::rowid_to_column(asv_tax, "seq_idx_in") |>
      dplyr::arrange(dplyr::across(all_of(c(TAX_RANKS, "seq_id")))) |>
      tibble::rowid_to_column("seq_idx") |>
      dplyr::select(seq_idx, seq_idx_in),
    deployment = "main"
  ),

  #### asv_taxsort_seq ####
  # `character` filename
  # sequence for each ASV
  tar_file_fast(
    asv_taxsort_seq,
    write_sequence(
      Biostrings::readDNAStringSet(asv_seq)[asv_taxsort$seq_idx_in],
      "sequences/04_denoised/asv_taxsort.fasta.gz",
      compress = TRUE
    ),
    deployment = "main"
  ),

  #### asv_taxsort_seq_index ####
  # `character` filename
  # sequence for each ASV
  tar_file_fast(
    asv_taxsort_seq_index,
    fastx_gz_index(asv_taxsort_seq),
    deployment = "main"
  ),

  if (protax_aligned) {
    list(
      #### aligned_taxsort_seq ####
      # `character` filename
      # aligned sequence for eah ASV, sorted by protax taxonomy
      tar_file_fast(
        aligned_taxsort_seq,
        write_sequence(
          # Use BString instead of DNAString because it will preserve case
          Biostrings::readBStringSet(asv_model_align)[asv_taxsort$seq_idx_in],
          "sequences/05_aligned/aligned_taxsort.fasta.gz",
          compress = TRUE
        )
      ),

      #### aligned_taxsort_seq_index ####
      # `character` filename
      # sequence for each ASV
      tar_file_fast(
        aligned_taxsort_seq_index,
        fastx_gz_index(aligned_taxsort_seq),
        deployment = "main"
      )
    )
  },

  #### seqbatch_result_map ####
  tar_fst_tbl(
    seqbatch_result_map,
    seqbatch |>
      dplyr::transmute(
        seq_idx,
        result = as.raw(
          0x10 * (!seq_idx %in% denovo_chimeras) +
            0x20 * (!seq_idx %in% ref_chimeras) +
            0x40 * (!seq_idx %in% spikes$seq_idx) +
            !!(
              if (do_model_filter)
                quote(0x80 * (seq_idx %in% asv_full_length))
              else
                0
            )
        )
      ),
    pattern = !!(
      if (do_model_filter) {
        quote(map(seqbatch, ref_chimeras, spikes, asv_full_length))
      } else {
        quote(map(seqbatch, ref_chimeras, spikes))
      }
    ),
    deployment = "main"
  ),

  #### asv_map ####
  tar_fst_tbl(
    asv_map,
    seqbatch_result_map |>
      dplyr::left_join(asv_names, by = "seq_idx") |>
      dplyr::select(seq_idx, result, seq_id)
  )
)

optimotu_plan <- c(optimotu_plan, asv_plan)
