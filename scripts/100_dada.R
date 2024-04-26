# DADA2 quality filtering and denoising for Sonja's spruce log metabarcoding data
# Brendan Furneaux
# Based on DADA2 analysis for GSSP from Jenni Hultman
# edits by Sten Anslan - account for reverse complementary oriented sequences and add UNCROSS2 tag-jumps filtering per run

############################################################
## TODO: ##
# - skip RC if no rc seqs! Currently = ERROR, if no rc seqs.
# - remove empty files after cutadapt in 02_trim/
############################################################

library(magrittr)
library(targets)
library(tarchetypes)

#### DADA2 analysis based on Jenni Hultman's pipeline for GSSP

orient_meta <- tibble::tibble(
  .orient = unique(sample_table$orient)
)

#### inner_dada_plan ####
# the "inner" dada plan consists of targets which need to be done separately
# for each orientation, if not all reads are in the same orientation.

inner_dada_plan <- list(
  ##### dada2_meta_{.orient}_{.seqrun} #####
  # grouped tibble:
  #  `seqrun` character; name of sequencing run (directory in sequences/01_raw)
  #  `sample` character; name of sample, based on parsing file name
  #  `fastq_R1` character; file name with path for raw R1 file
  #  `fastq_R2` character; file name with path for raw R2 file
  #  `trim_R1` character; file name with path for trimmed R1 file
  #  `trim_R2` character; file name with path for trimmed R2 file
  #  `filt_R1` character; file name with path for filtered R1 file
  #  `filt_R2` character; file name with path for filtered R2 file
  #  `sample_key`character; common prefix of trim_R1, trim_R2, filt_R1 and
  #      filt_R2; used as sample name by dada2 functions and read counts
  #
  # `sample_table` is defined in scripts/010_load_samples.R
  # in some cases it may have list columns, so it is not safe to use fst
  tar_group_size(
    dada2_meta,
    sample_table |>
      dplyr::filter(orient == .orient, seqrun == .seqrun) |>
      dplyr::select(seqrun, sample, fastq_R1, fastq_R2, trim_R1, trim_R2,
                    filt_R1, filt_R2, sample_key, any_of(cutadapt_paired_option_names)),
    size = 96
  ),

  ##### raw_R1_{.orient}_{.seqrun} #####
  # character: path and file name
  # raw reads, for dependency tracking
  tar_file_fast(
    raw_R1,
    file.path(raw_path, dada2_meta$fastq_R1),
    pattern = map(dada2_meta)
  ),

  ##### raw_read_counts_{.orient}_{.seqrun} #####
  # tibble:
  #  `fastq_file` character: file name of raw R1 file
  #  `raw_nread` integer: number of sequences in the file
  tar_fst_tbl(
    raw_read_counts,
    tibble::tibble(
      fastq_file = raw_R1,
      raw_nread = sequence_size(fastq_file)
    ),
    pattern = map(raw_R1)
  ),


  ##### trim_{.orient}_{.seqrun} #####
  # character: file names with path of trimmed read files (fastq.gz)
  #
  # remove adapters and barcodes
  # also do some preliminary quality filtering
  tar_file_fast(
    trim,
    withr::with_connection(
      list(logfile = file(sprintf("logs/trim_%s_%s.log", .seqrun, .orient), "w")),
      dplyr::group_by(
        dada2_meta,
        dplyr::pick(any_of(c("seqrun", cutadapt_option_names)))
      ) |>
        dplyr::group_map(
          ~ purrr::pmap(
            dplyr::transmute(
              .x,
              file_R1 = file.path(raw_path, fastq_R1),
              file_R2 = file.path(raw_path, fastq_R2),
              trim_R1 = trim_R1,
              trim_R2 = trim_R2
            ),
            cutadapt_paired_filter_trim,
            primer_R1 = ifelse(.orient == "fwd", trim_primer_R1, trim_primer_R2),
            primer_R2 = ifelse(.orient == "fwd", trim_primer_R2, trim_primer_R1),
            options = update(trim_options, .x),
            ncpu = local_cpus(),
            logfile = logfile
          ) |>
            unlist(),
          .keep = TRUE
        ) |>
        unlist()
    ),
    pattern = map(dada2_meta)
  ),

  ##### trim_read_counts_{.orient}_{.seqrun} #####
  # tibble:
  #  `trim_R1` character: file name with path of trimmed R1 file
  #  `trim_nread` integer: number of sequences in the file
  #
  # count of reads per sample after adapter trimming
  tar_fst_tbl(
    trim_read_counts,
    tibble::tibble(
      trim_R1 = purrr::keep(trim, endsWith, "_R1_trim.fastq.gz"),
      trim_nread = sequence_size(trim_R1)
    ),
    pattern = map(trim)
  ),

  # DADA2 quality filtering on read-pairs

  ##### filter_pairs_{.orient}_{.seqrun} #####
  # character: file names with path of filtered read files (fastq.gz)
  #
  # additional quality filtering on read-pairs
  tar_file_fast(
    filter_pairs,
    if (nrow(dada2_meta) > 0L) {
      file.create(c(dada2_meta$filt_R1, dada2_meta$filt_R2))
      dada2::filterAndTrim(
        fwd = purrr::keep(trim, endsWith, "_R1_trim.fastq.gz"),
        filt = dada2_meta$filt_R1,
        rev = purrr::keep(trim, endsWith, "_R2_trim.fastq.gz"),
        filt.rev = dada2_meta$filt_R2,
        maxEE = dada2_maxEE, # max expected errors (fwd, rev)
        rm.phix = TRUE, #remove matches to phiX genome
        compress = TRUE, # write compressed files
        multithread = local_cpus(),
        verbose = TRUE
      )
      # return file names for samples where at least some reads passed
      c(dada2_meta$filt_R1, dada2_meta$filt_R2) %>%
        purrr::keep(file.exists)
    } else {
      character()
    },
    pattern = map(dada2_meta, trim)
  ),

  ##### filt_read_counts_{.orient}_{.seqrun} #####
  # tibble:
  #  `filt_R1` character: file name with path of filtered R1 file
  #  `filt_nread` integer: number of sequences in the file
  #
  # count of reads after filtering
  tar_fst_tbl(
    filt_read_counts,
    tibble::tibble(
      filt_R1 = purrr::keep(filter_pairs, endsWith, "_R1_filt.fastq.gz"),
      filt_nread = sequence_size(filt_R1)
    ),
    pattern = map(filter_pairs)
  ),

  ##### map over R1 and R2 #####
  # inside the tar_map, every occurrence of `read` is replaced by "R1" or "R2"
  # the read name is also appended to all target names
  # so all of this gets done separately for forward and reverse reads.
  # pattern=map() means we are also keeping the different sequencing runs
  # separate
  tar_map(
    values = list(read = c("R1", "R2")),

    ###### filtered_{read}_{.orient}_{.seqrun} ######
    # character: path and file name of filtered reads; fastq.gz
    #
    # select only the files corresponding to the read we are working on
    tar_file_fast(
      filtered,
      purrr::keep(filter_pairs, endsWith, paste0(read, "_filt.fastq.gz")),
      pattern = map(filter_pairs),
      deployment = "main"
    ),

    ###### derep_{read}_{.orient}_{.seqrun} ######
    # list of dada2 `derep` objects
    #
    # dereplicate
    tar_target(
      derep,
      dada2::derepFastq(filtered, verbose = TRUE) %>%
        set_names(file_to_sample_key(filtered)),
      pattern = map(filtered)
    ),

    ###### err_{read}_{.orient}_{.seqrun} ######
    # list: see dada2::LearnErrors
    #
    # fit error profile
    tar_target(
      err,
      if (sum(!grepl("BLANK|NEG", filtered)) > 0L) {
        dada2::learnErrors(
          purrr::discard(filtered, grepl, pattern = "BLANK|NEG"),
          multithread = local_cpus(),
          verbose = TRUE
        )
      } else {
        NULL
      }
    ),

    ###### denoise_{read}_{.orient}_{.seqrun} ######
    # list of dada2 `dada` objects
    tar_target(
      denoise,
      if (is.null(err)) {
        NULL
      } else {
        dada2::dada(derep, err = err, multithread = local_cpus(), verbose = TRUE)
      },
      pattern = map(derep)
    )
  ),

  ##### merged_{.orient}_{.seqrun} #####
  # list of data.frame; see dada2::mergePairs
  #
  # Merge paired reads and make a sequence table for each sequencing run
  tar_target(
    merged,
    if (is.null(denoise_R1) | is.null(denoise_R2)) {
      list()
    } else {
      dada2::mergePairs(
        denoise_R1,
        derep_R1,
        denoise_R2,
        derep_R2,
        minOverlap = 10,
        maxMismatch = 1,
        verbose=TRUE
      )
    },
    pattern = map(denoise_R1, derep_R1, denoise_R2, derep_R2)
  ),

  ##### seqtable_raw_{.orient}_{.seqrun} #####
  # `tibble` with columns:
  #   `sample` (character) sample name as given in sample_table$sample_key
  #   `seq_idx` (integer) index of a sequence in seq_all
  #   `nread` (integer) number of reads
  tar_fst_tbl(
    seqtable_raw,
    make_mapped_sequence_table(merged, seq_all, rc = .orient == "rev"),
    pattern = map(merged)
  ),

  ##### dada_map_{.orient}_{.seqrun} #####
  # map the raw reads to the nochim ASVs
  # indexes are per-sample
  #
  # a tibble:
  #  sample: character identifying the sample, as in sample_table
  #  raw_idx: integer index of read in the un-rarified fastq file
  #  seq_idx: integer index of ASV in seq_all
  #  flags: raw, bits give presence/absence of the read after different stages:
  #    0x01: trim
  #    0x02: filter
  #    0x04: denoise & merge
  #    0x08: tag-jump removal (if performed)
  if (isTRUE(do_uncross) && nrow(orient_meta) == 1L) {
    tar_target(
      dada_map,
      mapply(
        FUN = seq_map,
        sample = dada2_meta$sample_key,
        fq_raw = file.path(raw_path, dada2_meta$fastq_R1),
        fq_trim = dada2_meta$trim_R1,
        fq_filt = dada2_meta$filt_R1,
        dadaF = denoise_R1,
        derepF = derep_R1,
        dadaR = denoise_R2,
        derepR = derep_R2,
        merged = merged,
        MoreArgs = list(
          seq_all = seq_all,
          rc = .orient == "rev"
        ),
        SIMPLIFY = FALSE
      ) |>
        purrr::list_rbind(
          ptype = tibble::tibble(
            sample = character(),
            raw_idx = integer(),
            seq_idx = integer(),
            flags = raw()
          )
        ) |>
        add_uncross_to_seq_map(seqtable_raw, uncross),
      pattern = map(dada2_meta, denoise_R1, derep_R1, denoise_R2, derep_R2, merged)
    )
  } else {
    tar_target(
      dada_map,
      mapply(
        FUN = seq_map,
        sample = dada2_meta$sample_key,
        fq_raw = file.path(raw_path, dada2_meta$fastq_R1),
        fq_trim = dada2_meta$trim_R1,
        fq_filt = dada2_meta$filt_R1,
        dadaF = denoise_R1,
        derepF = derep_R1,
        dadaR = denoise_R2,
        derepR = derep_R2,
        merged = merged,
        MoreArgs = list(
          seq_all = seq_all,
          rc = .orient == "rev"
        ),
        SIMPLIFY = FALSE
      ) |>
        purrr::list_rbind(
          ptype = tibble::tibble(
            sample = character(),
            raw_idx = integer(),
            seq_idx = integer(),
            flags = raw()
          )
        ),
      pattern = map(dada2_meta, denoise_R1, derep_R1, denoise_R2, derep_R2, merged)
    )
  }
)

if (nrow(orient_meta) > 1L) {
  inner_dada_plan <- tar_map(
    values = orient_meta,
    inner_dada_plan
  )
} else {
  .orient <- orient_meta$.orient
}

seqrun_meta <- tibble::tibble(
  .seqrun = unique(sample_table$seqrun)
)

#### seqrun_plan ####

seqrun_plan <- tar_map(
  values = seqrun_meta,
  inner_dada_plan,

  ##### seq_merged #####
  # `character` vector
  #
  # all unique merged ASV sequences for each seqrun
  if (nrow(orient_meta) == 2) {
    tar_target(
      seq_merged,
      unique(
        as.character(
          unlist(
            c(
              lapply(merged_fwd, \(x) x$sequence),
              lapply(merged_rev, \(x) dada2::rc(x$sequence))
            )
          )
        )
      )
    )
  } else if (.orient == "fwd") {
    tar_target(
      seq_merged,
      unique(
        as.character(
          unlist(
            lapply(merged, \(x) x$sequence)
          )
        )
      )
    )
  } else if (.orient == "rev") {
    tar_target(
      seq_merged,
      unique(
        as.character(
          unlist(
            lapply(merged, \(x) dada2::rc(x$sequence))
          )
        )
      )
    )
  },

  if (nrow(orient_meta) == 2) {
    ##### seqtable_raw_{.seqrun} #####
    # `tibble`:
    #   `sample (character) - sample name as given in sample_table$sample_key
    #   `seq_idx` (integer) - index of a sequence in seq_all
    #   `nread` (integer) number of reads
    tar_target(
      seqtable_raw,
      dplyr::bind_rows(seqtable_raw_fwd, seqtable_raw_rev) |>
        dplyr::summarize(nread = sum(nread), .by = c(sample, seq_idx))
    )
  },

  ##### denoise_read_counts_{.seqrun} #####
  # tibble:
  #  `sample_key` character: as `sample_table$sample_key`
  #  `denoise_nread` integer: number of sequences in the sample after denoising
  tar_fst_tbl(
    denoise_read_counts,
    dplyr::summarize(
      seqtable_raw,
      denoise_nread = sum(nread),
      .by = sample
    ) |>
      dplyr::rename(sample_key = sample)
  ),

  if (isTRUE(do_uncross)) {
    list(
      ##### uncross_{.seqrun} #####
      # `tibble`:
      #   `sample (character) - sample name as given in sample_table$filt_key
      #   `nread` (integer) - number of reads in that sample (for this ASV)
      #   `total` (integer) - number of reads of that ASV across all samples (for
      #     this ASV)
      #   `uncross` (numeric) - UNCROSS score
      #   `is_tag_jump` (logical) - whether this occurrence is considered a likely
      #     tag jump
      #
      # `seq_idx` is not explicitly included; however the order of rows matches
      # `seqtable_raw`.
      #
      # remove tag-jumps (UNCROSS2)
      tar_fst_tbl(
        uncross,
        remove_tag_jumps(
          seqtable_raw, # raw_ASV_table
          tagjump_options$f, # f-value (expected cross-talk rate)
          tagjump_options$p, # p-value (power to rise the exponent)
          "seq_idx" # name of column which uniquely identifies the sequence
        )
      ),

      ##### uncross_summary_{.seqrun} #####
      # `tibble`:
      #   `sample (character) - sample name as given in sample_table$filt_key
      #   `Total_reads` (integer) - total reads in the sample (all ASVs)
      #   `Number_of_TagJump_Events` (integer) - number of ASVs in the sample which
      #     are considered tag jumps.
      #   `TagJump_reads` (integer) - number of reads in the sample which belong to
      #     tag-jump ASVs.
      tar_fst_tbl(
        uncross_summary,
        summarize_uncross(uncross)
      ),

      ##### uncross_read_counts_{.seqrun} #####
      # tibble:
      #  `sample_key` character: as `sample_table$sample_key`
      #  `uncross_nread` integer: number of sequences in the sample after denoising
      tar_fst_tbl(
        uncross_read_counts,
        uncross_summary |>
          dplyr::transmute(
            sample_key = sample,
            uncross_nread = Total_reads - TagJump_reads
          )
      ),
      ##### seqtable_uncross_{.seqrun} #####
      # `tibble`:
      #   `sample (character) - sample name as given in sample_table$sample_key
      #   `seq_idx` (integer) - index of a sequence in seq_all
      #   `nread` (integer) number of reads
      tar_fst_tbl(
        seqtable_uncross,
        seqtable_raw[!uncross$is_tag_jump,]
      )
    )
  },

  if (nrow(orient_meta) == 2) {
    if (isTRUE(do_uncross)) {
      tar_fst_tbl(
        dada_map,
        merge_seq_maps(dada_map_fwd, dada_map_rev) |>
          add_uncross_to_seq_map(seqtable_raw, uncross)
      )
    } else {
      tar_fst_tbl(
        dada_map,
        merge_seq_maps(dada_map_fwd, dada_map_rev)
      )
    }
  },

  ##### bimera_table_{.seqrun} #####
  # tibble:
  #  `nflag` integer: number of samples in which the sequence was considered
  #    chimeric
  #  `nsam` integer: number of samples in which the sequence occurred
  #  `seq` character: sequence
  #
  # find denovo chimeric sequences in each sample independently
  tar_fst_tbl(
    bimera_table,
    bimera_denovo_table(
      !!(if (isTRUE(do_uncross)) quote(seqtable_uncross) else quote(seqtable_raw)),
      seq_all,
      allowOneOff=TRUE,
      multithread=local_cpus()
    )
  )
)

seq_all_file <- "sequences/04_denoised/all_asv.fasta.gz"

#### dada_plan ####
dada_plan <- list(
  seqrun_plan,

  ##### seq_all #####
  # `character` vector
  #
  # all unique ASV sequences, across all seqruns
  tar_file_fast(
    seq_all,
    {
      old_seqs <-
        if (file.exists(seq_all_file)) {
          Biostrings::readDNAStringSet(seq_all_file)
        } else {
          Biostrings::DNAStringSet()
        }
      uniqs <- unique(!!tar_map_c(seqrun_plan$seq_merged))
      seqs <- c(
        old_seqs,
        Biostrings::DNAStringSet(uniqs[is.na(BiocGenerics::match(uniqs, old_seqs))])
      )
      names(seqs) <- seq_along(seqs)
      write_and_return_file(
        seqs,
        file = seq_all_file,
        compress = "gzip",
        compression_level = 9
      )
    }
  ),

  ##### denovo_chimeras #####
  # `integer` vector - index of ASV sequences which were found to be chimeric
  #
  # calculate consensus chimera calls across all seqruns
  tar_target(
    denovo_chimeras,
    combine_bimera_denovo_tables(
      !!tar_map_bind_rows(seqrun_plan$bimera_table)
    )
  ),

  ##### seqtable_merged #####
  # `tibble`:
  #   `sample (character) - sample name as given in sample_table$sample_key
  #   `seq_idx` (integer) - index of a sequence in seq_all
  #   `nread` (integer) number of reads
  tar_fst_tbl(
    seqtable_merged,
    !!tar_map_bind_rows(
      seqrun_plan,
      if (isTRUE(do_uncross)) "seqtable_uncross" else "seqtable_raw"
    )
  ),

  ##### nochim1_read_counts #####
  # tibble:
  #  `sample_key` character: as `sample_table$sample_key`
  #  `nochim1_nread` integer: number of sequences in the sample after first
  #    chimera filtering
  tar_fst_tbl(
    nochim1_read_counts,
    dplyr::filter(seqtable_merged, !seq_idx %in% denovo_chimeras) |>
    dplyr::summarize(nochim1_nread = sum(nread), .by = sample) |>
      dplyr::rename(sample_key = sample)
  )
)

optimotu_plan <- c(optimotu_plan, dada_plan)
