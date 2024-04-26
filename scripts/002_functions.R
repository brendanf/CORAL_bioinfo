bimera_denovo_table <- function(
    seqtab,
    seq = NULL,
    minFoldParentOverAbundance = 1.5,
    minParentAbundance = 2,
    allowOneOff = FALSE,
    minOneOffParentDistance = 4,
    maxShift = 16,
    multithread = FALSE,
    ...
) {
  UseMethod("bimera_denovo_table", seqtab)
}

bimera_denovo_table.matrix <- function(
  seqtab,
  seqs = colnames(seqtab),
  minFoldParentOverAbundance = 1.5,
  minParentAbundance = 2,
  allowOneOff = FALSE,
  minOneOffParentDistance = 4,
  maxShift = 16,
  multithread = FALSE
) {
  if (is.null(seqs)) seqs <- colnames(seqtab)
  if (isTRUE(multithread)) {
    RcppParallel::setThreadOptions(numThreads = "auto")
  } else if (isFALSE(multithread)) {
    RcppParallel::setThreadOptions(numThreads = 1)
  } else {
    assertthat::assert_that(
      assertthat::is.count(multithread),
      msg = "argument 'multithread' must be TRUE, FALSE, or a positive integer."
    )
    RcppParallel::setThreadOptions(numThreads = multithread)
  }
  dada2:::C_table_bimera2(
    mat = seqtab,
    seqs = seqs,
    min_fold = minFoldParentOverAbundance,
    min_abund = minParentAbundance,
    allow_one_off = allowOneOff,
    min_one_off_par_dist = minOneOffParentDistance,
    match = dada2::getDadaOpt("MATCH"),
    mismatch = dada2::getDadaOpt("MISMATCH"),
    gap_p = dada2::getDadaOpt("GAP_PENALTY"),
    max_shift = maxShift
  ) %>%
    tibble::as_tibble() %>%
    tibble::add_column(seq = seqs)
}

bimera_denovo_table.data.frame <- function(
    seqtab,
    seqs = NULL,
    minFoldParentOverAbundance = 1.5,
    minParentAbundance = 2,
    allowOneOff = FALSE,
    minOneOffParentDistance = 4,
    maxShift = 16,
    multithread = FALSE
) {
  seq_col <- intersect(c("seq", "seq_id", "seq_idx"), names(seqtab))
  if (length(seq_col) == 0) {
    stop("seqtab must have at least one of columns 'seq', 'seq_id', or 'seq_idx'")
  }
  seq_col <- seq_col[1]
  if (seq_col != "seq") {
    if (length(seqs) == 1 && file.exists(seqs)) {
      seqs <- Biostrings::readDNAStringSet(seqs)
    }
    if (methods::is(seqs, "XStringSet")) seqs <- as.character(seqs)
    if (!is.character(seqs)) {
      stop("'seqs' must be a filename of a fasta file, an XStringSet, or a character")
    }
    if (identical(seq_col, "seq_id") && !rlang::is_named(seqs)) {
      stop("'seqs' must be named if sequences are identified by 'seq_id' in 'seqtab'")
      checkmate::assert_subset(names(seqs), seqtab$seq_id)
    }
  }
  n_asv <- dplyr::n_distinct(seqtab[[seq_col]])
  n_sample <- dplyr::n_distinct(seqtab$sample)
  # max seqtable size for one partition is 1 Gb (== 2^30 bytes)
  # (not including sequences)
  # R integers are 32 bit (== 4 bytes)
  # If there are more than 250M ASVs in a single sample, then the matrix ends
  # up larger (but at that point the size of the sequences themselves is a bigger problem!)
  n_partition <- ceiling(n_asv*n_sample*4/2^30)
  sample_splits <-
    split(unique(seqtab$sample), rep(seq_len(n_partition), length.out = n_sample))
  out <- list()
  for (s in sample_splits) {
    m <- dplyr::filter(seqtab, sample %in% s) |>
      tidyr::pivot_wider(
        names_from = all_of(seq_col),
        values_from = nread,
        values_fill = list(nread = 0L)
      ) |>
      tibble::column_to_rownames("sample") |>
      as.matrix()

    switch(seq_col,
      seq_id = colnames(m) <- seqs[colnames(m)],
      seq_idx = colnames(m) <- seqs[as.integer(colnames(m))]
    )

    out_m <- bimera_denovo_table.matrix(
      seqtab = m,
      minFoldParentOverAbundance = minFoldParentOverAbundance,
      minParentAbundance = minParentAbundance,
      allowOneOff = allowOneOff,
      minOneOffParentDistance = minOneOffParentDistance,
      maxShift = maxShift,
      multithread = multithread
    )
    switch(
      seq_col,
      seq_idx = out_m$seq <- match(out_m$seq, seqs),
      seq_id = out_m$seq <- names(seqs)[match(out_m$seq, seqs)]
    )
    names(out_m)[3] <- "seq_idx"

    out <- c(
      out,
      list(out_m)
    )
  }

  dplyr::bind_rows(out) |>
    dplyr::summarize(nflag = sum(nflag), nsam = sum(nsam), .by = all_of(seq_col))
}

combine_bimera_denovo_tables <- function(
  bimdf,
  minSampleFraction = 0.9,
  ignoreNNegatives = 1L,
  verbose = FALSE
) {
  seq_col <- intersect(c("seq", "seq_id", "seq_idx"), names(bimdf))
  if (length(seq_col) == 0) {
    stop("bimdf must have at least one of columns 'seq', 'seq_id', or 'seq_idx'")
  }
  seq_col <- seq_col[1]

  bimdf <- dplyr::summarize(
    bimdf,
    dplyr::across(everything(), sum),
    .by = any_of(seq_col)
  )
  ## This snippet modified from DADA2
  bims.out <- with(
    bimdf,
    nflag >= nsam | (nflag > 0 & nflag >= (nsam - ignoreNNegatives) * minSampleFraction)
  )
  if (verbose)
    message("Identified ", sum(bims.out), " bimeras out of ",
            length(bims.out), " input sequences.")
  ## end snippet from DADA2
  bimdf[[seq_col]][bims.out]
}

remove_bimera_denovo_tables <- function(
  seqtabs,
  bimdf,
  minSampleFraction = 0.9,
  ignoreNNegatives = 1L,
  verbose = FALSE
) {
  bims.out <- combine_bimera_denovo_tables(
    bimdf,
    minSampleFraction = minSampleFraction,
    ignoreNNegatives = ignoreNNegatives,
    verbose = verbose
  )
  dplyr::bind_rows(seqtabs) |>
    dplyr::filter(!bims.out[seq])
}

#' Map the fate of indivdual reads through trimming and filtering steps
#'
#' @param fq_raw (character) raw fastq file
#' @param fq_trim (character) trimmed fastq file
#' @param fq_filt (character) filtered fastq file
#'
#' @return `tibble` with columns:
#'   `raw_idx` (integer) index of the sequence in fastq_raw; or if sequence names
#'     in the fastq files are hex encoded integers (e.g., during subsampling to
#'     indicate the original index) then `seq_id` comes from the sequence names.
#'   `trim_idx` (integer) index of the sequence in fastq_trim
#'   `filt_idx` (integer) index of the sequence in fastq_filt
fastq_seq_map <- function(fq_raw, fq_trim, fq_filt) {
  out <- tibble::tibble(
    raw_idx = fastq_names(fq_raw)
  ) |>
    dplyr::left_join(
      tibble::enframe(fastq_names(fq_trim), value = "raw_idx", name = "trim_idx"),
      by = "raw_idx"
    ) |>
    dplyr::left_join(
      tibble::enframe(fastq_names(fq_filt), value = "raw_idx", name = "filt_idx"),
      by = "raw_idx"
    )
  if (nrow(out) > 100) {
    if (!all(grepl("^[1-9a-f]+$", out$raw_idx[1:100]))) {
      out$raw_idx <- seq_along(out$raw_idx)
      return(out)
    }
  }
  if (all(grepl("^[1-9a-f]+$", out$raw_idx))) {
    out$raw_idx <- as.integer(paste0("0x", out$raw_idx))
  } else {
    out$raw_idx <- seq_along(out$raw_idx)
  }
  out
}

#' Map the fate of individual reads through dada2 dereplication, denoising, and merge.
#'
#' @param dadaF (`dada2::dada-class` object or list of such objects) denoised
#' forward reads
#' @param derepF (`dada2::derep-class` object or list of such objects)
#' dereplicated forward reads
#' @param dadaR (`dada2::dada-class` object or list of such objects) denoised
#' reverse reads
#' @param derepR (`dada2::derep-class` object or list of such objects)
#' dereplicated reverse reads
#' @param merged (`data.frame` returned by `dada2::mergePairs()` or list of
#' such objects) results of merginf the denoised reads in dadaF and dadaR
#'
#' @return a `data.frame` with three columns:
#'   - `fwd_idx` (integer) index of forward ASV in `dadaF`
#'   - `rev_idx` (integer) index of reverse ASV in `dadaR`
#'   - `merge_idx` (integer) row index of merged ASV `merged`
#' Each row of this `data.frame` represents a single read in the fastq files
#' originally passed to `dada2::derepFastq()`, and the rows are in the same
#' order as the reads.
#' If the inputs were lists, then the output is a list of `data.frame`s as
#' described above.
dada_merge_map <- function(dadaF, derepF, dadaR, derepR, merged) {
  if (all(
    methods::is(dadaF, "dada"),
    methods::is(dadaR, "dada"),
    methods::is(derepF, "derep"),
    methods::is(derepR, "derep"),
    methods::is(merged, "data.frame")
  )) {
    tibble::tibble(
      fwd_idx = dadaF$map[derepF$map],
      rev_idx = dadaR$map[derepR$map]
    ) |>
      dplyr::left_join(
        tibble::rowid_to_column(merged[c("forward", "reverse")], "merge_idx"),
        by = c("fwd_idx" = "forward", "rev_idx" = "reverse")
      )
  } else if (all(
    rlang::is_bare_list(dadaF),
    rlang::is_bare_list(dadaR),
    rlang::is_bare_list(derepF),
    rlang::is_bare_list(derepR),
    rlang::is_bare_list(merged)
  )) {
    purrr::pmap(list(dadaF, derepF, dadaR, derepR, merged), dada_merge_map)
  }
}

#' Map the fate of individual reads through merging to find unique reads
#'
#' @param sample (character) name of the sample
#' @param fq_raw (character) name of the raw fastq R1 file
#' @param fq_trim (character) name of the trimmed fastq R1 file
#' @param fq_file (character) name of the filtered fastq R1 file
#' @param dadaF (dada-object) denoised R1
#' @param derepF (derep-object) dereplicated R1
#' @param dadaR (dada-object) denoised R2
#' @param derepR (dada-object) dereplicated R2
#' @param merged (data.frame as returned by `dada2::mergePairs`) result of merging dadaF and dadaR
#' @param seq_all (character) unique ASV sequences
#' @param rc (logical) if TRUE, sequences in `merged` are reverse-complemented
#'  relative to seq_all.
#'
#' @return `data.frame` with columns:
#'   `sample` (character) the sample name
#'   `raw_idx` (integer) the index of the sequence in the raw file; see `seq_map()`
#'   `seq_idx` (integer) the index of the sequence in seq_all
#'   `flags` (raw) bitset indicating the presence of the sequence at different stages:
#'    0x01 = trimmed
#'    0x02 = filtered
#'    0x04 = denoised & merged
seq_map <- function(sample, fq_raw, fq_trim, fq_filt, dadaF, derepF, dadaR, derepR, merged, seq_all, rc = FALSE) {
  seq_map <- fastq_seq_map(fq_raw, fq_trim, fq_filt)
  dada_map <- dada_merge_map(dadaF, derepF, dadaR, derepR, merged)
  seq_map$dada_idx <-
  seq_map$seq_idx <- match(merged$sequence, seq_all)[dada_map$merge_idx[seq_map$filt_idx]]
  dplyr::transmute(
    seq_map,
    sample = sample,
    raw_idx,
    seq_idx,
    flags = as.raw(
      ifelse(is.na(trim_idx), 0, 0x01) +
        ifelse(is.na(filt_idx), 0, 0x02) +
        ifelse(is.na(dada_idx), 0, 0x04)
    )
  )
}

merge_seq_maps <- function(seqmap_fwd, seqmap_rev) {
  dplyr::full_join(
    seqmap_fwd,
    seqmap_rev,
    by = c("sample", "raw_idx"),
    suffix = c("_fwd", "_rev")
  ) |>
    dplyr::transmute(
      sample,
      raw_idx,
      seq_idx = dplyr::coalesce(seq_idx_fwd, seq_idx_rev),
      flags = flags_fwd | flags_rev
    )
}

add_uncross_to_seq_map <- function(seqmap, seqtable_raw, uncross) {
  dplyr::left_join(
    seqmap,
    tibble::tibble(
      sample = seqtable_raw$sample,
      seq_idx = seqtable_raw$seq_idx,
      is_tag_jump = uncross$is_tag_jump
    ),
    by = c("sample", "seq_idx")
  ) |>
    dplyr::transmute(
      sample,
      raw_idx,
      seq_idx,
      flags = flags | as.raw(ifelse(is.na(is_tag_jump) | is_tag_jump, 0, 0x08))
    )
}

sort_seq_table <- function(seqtable, ...) {
  UseMethod("sort_seq_table", seqtable)
}

sort_seq_table.matrix <- function(seqtable, ...) {
  colorder <- order(
    -colSums(seqtable > 0), # prevalence, highest to lowest
    -colSums(seqtable), # abundance, highest to lowest
    -apply(seqtable, 2, var), # variance, highest to lowest
    seqhash(colnames(seqtable)) # hash of sequence (pseudorandom but stable)
  )
  if (is.null(attr(seqtable, "map"))) {
    seqtable[order(rownames(seqtable)), colorder]
  } else {
    structure(
      seqtable[order(rownames(seqtable)), colorder],
      map = dplyr::mutate(attr(seqtable, "map"), seq_id_out = order(colorder)[seq_id_out])
    )
  }
}

# if possible, returns an ordering permutation over the sequences
sort_seq_table.data.frame <- function(seqtable, seqs = NULL, abund_col = "nread", ...) {
  # TODO: add some verification here
  abund <- as.symbol(abund_col)
  seqorder <- dplyr::summarize(
    seqtable,
    prevalence = dplyr::n(),
    abundance = sum(!!abund),
    variance = if (prevalence == 1) 0 else var(!!abund),
    .by = any_of(c("seq", "seq_idx", "seq_id"))
  )
  seqorder$hash <-
    if ("seq" %in% names(seqorder)) {
      seqhash(seqorder$seq)
    } else if ("seq_idx" %in% names(seqorder)) {
      hash_sequences(seqs, use_names = FALSE)[as.integer(seqorder$seq_idx)]
    } else if ("seq_id" %in% names(seqorder)) {
      unname(hash_sequences(seqs, use_names = TRUE)[seqorder$seq_id])
    }
  out <- order(
    -seqorder$prevalence,
    -seqorder$abundance,
    -seqorder$variance,
    seqorder$hash
  )
  if ("seq" %in% names(seqorder)) {
    seqorder$seq[out]
  } else if ("seq_idx" %in% names(seqorder)) {
    out
  } else if ("seq_id" %in% names(seqorder)) {
    seqorder$seq_id[out]
  }
}

#' Calculate clustering thresholds for each taxon, falling back to its ancestor
#' taxa as necessary
#'
#' @param rank (`character` string) the rank within which clustering will be
#' performed
#' @param conf_level (`character` string) as in `fmeasure_optima`
#' @param taxon_table (`data.frame`) taxonomy table; column "seq_id" gives the
#' sequence ID, and columns {ROOT_RANK} to {TIP_RANK} (e.g., "kingdom" to
#' "species") give the taxonomy at each rank.
#' @param fmeasure_optima (`data.frame`) optimum clustering thresholds within
#' various taxa; column "rank" gives the rank which is approximated by
#' clustering; "superrank" gives the rank of the taxon within which the
#' clustering threshold was optimized; "supertaxon" gives that taxon name;
#' "conf_level" gives a string description of the confidence level threshold for
#' taxonomic assignments; "threshold" gives the optimum clustering threshold;
#' "f_measure" gives the F measure at the optimum threshold.
#' @param default (`character string`) default taxon to define threshold to use
#' when taxonomy is unknown. default: INGROUP_TAXON
#'
#'
#' @return (named `numeric`, where names are taxa and values are clustering
#' thresholds)
#' unknown sequences, grouped by taxonomy at the parent rank. Sequences where
#' the parent rank is also unknown are in the item named `"_NA_"`.
calc_taxon_thresholds <- function(rank, conf_level, taxon_table,
                                  fmeasure_optima, default = INGROUP_TAXON) {
  rank_name <- rlang::sym(rank)
  dplyr::select(taxon_table, !!ROOT_RANK:!!rank_name) %>%
    dplyr::filter(!is.na(!!rank_name)) %>%
    unique() %>%
    purrr::reduce(
      c(superranks(rank), rank),
      function(thresholds, r) {
        dplyr::left_join(
          thresholds,
          dplyr::filter(
            fmeasure_optima,
            rank == subranks(!!rank)[1],
            superrank == r,
            conf_level == !!conf_level
          ) %>%
            dplyr::select(
              !!r := supertaxon,
              !!paste0("threshold_", r) := threshold
            ),
          by = r
        )
      },
      .init = .
    ) %>%
    dplyr::transmute(
      !!rank_name := !!rank_name,
      threshold = {.} %>%
        dplyr::select(dplyr::starts_with("threshold")) %>%
        rev() %>%
        do.call(dplyr::coalesce, .)
    ) %>%
    tibble::deframe() %>%
    c("_NA_" = dplyr::filter(
      fmeasure_optima,
      rank == subranks(!!rank)[1],
      supertaxon == default,
      conf_level == !!conf_level
    )$threshold)
}

threshold_as_dist <- function(thresholds) {
  if (any(thresholds > 50)) {
    1 - 0.01 * thresholds
  } else if (any(thresholds > 1)) {
    0.01 * thresholds
  } else if (any(thresholds > 0.5) ) {
    1 - thresholds
  } else {
    thresholds
  }
}

#' Calculate clustering thresholds for each taxon, falling back to its ancestor
#' taxa as necessary
#'
#' @param rank (`character` string) the rank within which clustering will be
#' performed
#' @param conf_level (`character` string) as in `fmeasure_optima`
#' @param taxon_table (`data.frame`) taxonomy table; column "seq_id" gives the
#' sequence ID, and columns {ROOT_RANK} to {TIP_RANK} (e.g., "kingdom" to
#' "species") give the taxonomy at each rank.
#' @param fmeasure_optima (`data.frame`) optimum clustering thresholds within
#' various taxa; column "rank" gives the rank which is approximated by
#' clustering; "superrank" gives the rank of the taxon within which the
#' clustering threshold was optimized; "supertaxon" gives that taxon name;
#' "conf_level" gives a string description of the confidence level threshold for
#' taxonomic assignments; "threshold" gives the optimum clustering threshold;
#' "f_measure" gives the F measure at the optimum threshold.
#' @param default (`character string`) default taxon to define threshold to use
#' when taxonomy is unknown. default: INGROUP_TAXON
#'
#'
#' @return (named `list` of `double` vectors)
calc_subtaxon_thresholds <- function(rank, conf_level, taxon_table,
                                  fmeasure_optima, default = INGROUP_TAXON) {
  rank_name <- rlang::sym(rank)
  dplyr::select(taxon_table, {{ROOT_RANK}}:{{rank_name}}) %>%
    tidyr::crossing(subrank = subranks(rank)) %>%
    dplyr::filter(!is.na(!!rank_name)) %>%
    unique() %>%
    purrr::reduce(
      c(superranks(rank), rank),
      function(thresholds, r) {
        dplyr::left_join(
          thresholds,
          dplyr::filter(
            fmeasure_optima,
            rank %in% subranks(!!rank),
            superrank == r,
            conf_level == !!conf_level
          ) %>%
            dplyr::select(
              subrank = rank,
              !!r := supertaxon,
              !!paste0("threshold_", r) := threshold
            ),
          by = c("subrank", r)
        )
      },
      .init = .
    ) %>%
    dplyr::transmute(
      subrank = rank2factor(subrank),
      !!rank_name := !!rank_name,
      threshold = {.} %>%
        dplyr::select(dplyr::starts_with("threshold")) %>%
        rev() %>%
        do.call(dplyr::coalesce, .) %>%
        threshold_as_dist()
    ) %>%
    dplyr::arrange(subrank) %>%
    split(.[[rank]]) %>%
    lapply(dplyr::select, !any_of(rank)) %>%
    lapply(tibble::deframe) %>%
    lapply(cummax) %>%
    c(
      "_NA_" = dplyr::filter(
        fmeasure_optima,
        rank %in% subranks(!!rank),
        supertaxon == default,
        conf_level == !!conf_level
      ) %>% dplyr::transmute(rank = rank2factor(rank), threshold_as_dist(threshold)) %>%
        dplyr::arrange(rank) %>%
        tibble::deframe() %>%
        cummax() %>%
        list()
    )
}

parse_protax_nameprob <- function(nameprob, id_is_int = FALSE) {
  checkmate::assert_flag(id_is_int)
  id_col <- if (isTRUE(id_is_int)) "seq_idx" else "seq_id"
  id_col_name <- as.symbol(id_col)
  set_names(nameprob, basename(nameprob)) |>
    lapply(readLines) |>
    tibble::enframe() |>
    tidyr::extract(
      name,
      into = "rank",
      regex = "query(\\d+)\\.nameprob",
      convert = TRUE
    ) |>
    tidyr::unchop(value) |>
    dplyr::mutate(
      rank = rank2factor(TAX_RANKS[rank]),
      value = gsub("([^\t]+)\t([0-9.]+)", "\\1:\\2", value) %>%
        gsub("(:[0-9.]+)\t", "\\1;", .)
    ) |>
    tidyr::separate(value, into = c(id_col, "nameprob"), sep = "\t", fill = "right") |>
    tidyr::separate_rows(nameprob, sep = ";") |>
    tidyr::separate(nameprob, into = c("name", "prob"), sep = ":", convert = TRUE) |>
    tidyr::extract(name, into = c("parent_taxonomy", "taxon"), regex = "(.+),([^,]+)$") |>
    dplyr::mutate(
      !!id_col_name := if (id_is_int) as.integer(!!id_col_name) else !!id_col_name,
      prob = ifelse(is.na(taxon), 0, prob)
    ) |>
    dplyr::arrange(!!id_col_name, dplyr::desc(rank), dplyr::desc(prob))
}

# combine tip classifications to build a full PROTAX taxonomy

build_taxonomy <- function(...) {
  tax <- tibble::tibble(
    classification = union(...),
    rank = ifelse(
      classification == "root",
      0L,
      stringr::str_count(classification, stringr::fixed(",")) + 1L),
    parent = ifelse(
      rank <= 1,
      "root",
      sub(",[^,]+$", "", classification)
    )
  ) %>%
    split(.$rank)
  tax[[1]]$taxon_id = 0L
  tax[[1]]$prior = 1
  tax[[1]]$parent_id = 0L
  tax[[2]]$taxon_id = 1L:2L
  tax[[2]]$prior = c(0.99, 0.01)

  for (r in length(tax):3L) {
    if (r == length(tax)) {
      tax[[r]]$prior = 0.99/nrow(tax[[r]])
    } else {
      tax[[r]]$prior <- NULL
      tax[[r]] <- dplyr::left_join(
        tax[[r]],
        dplyr::group_by(tax[[r+1]], parent) %>%
          dplyr::summarise(prior = sum(prior)) %>%
          dplyr::rename(classification = parent),
        by = "classification"
      )
    }
  }
  for (r in 2L:length(tax)) {
    tax[[r]]$taxon_id <- seq_len(nrow(tax[[r]])) + max(tax[[r-1L]]$taxon_id)
    tax[[r]]$parent_id <- NULL
    tax[[r]] <- dplyr::left_join(
      tax[[r]],
      dplyr::select(
        tax[[r-1]],
        parent = classification,
        parent_id = taxon_id
      ),
      by = "parent"
    )
  }
  dplyr::bind_rows(tax) %>%
    dplyr::select(taxon_id, parent_id, rank, classification, prior)
}

# Format a classification for use as a Sintax reference DB
sintax_format <- function(s) {
  s <- sub(",", ";p:", s, fixed = TRUE)
  s <- sub(",", ";c:", s, fixed = TRUE)
  s <- sub(",", ";o:", s, fixed = TRUE)
  s <- sub(",", ";f:", s, fixed = TRUE)
  s <- sub(",", ";g:", s, fixed = TRUE)
  s <- sub(",", ";s:", s, fixed = TRUE)
  s <- chartr(";", ",", s)
  paste0("tax=d:", s, ";")
}

build_taxonomy_new <- function(...) {
  tax <- tibble::tibble(classification = unique(c(...))) |>
    tidyr::separate_wider_delim(
      1,
      delim = ",",
      names = TAX_RANKS,
      too_few = "align_start"
    ) |>
    dplyr::mutate({{ROOT_RANK}} := ifelse({{ROOT_RANK}}=="root", NA_character_, {{ROOT_RANK}}))

  # remove all taxa with children
  for (rank in 6:1) {
    rankname <- rlang::sym(TAX_RANKS[rank + 1])
    tax <- dplyr::filter(
      tax,
      if (any(!is.na({{rankname}} ))) !is.na({{rankname}}) else TRUE,
      .by = all_of(TAX_RANKS[1:rank])
    )
  }
  tax$n <- 1L
  tax_out <- list()
  tax_out[[1]] <- tibble::tibble(
    taxon_id = 0L,
    parent_id = 0L,
    rank = 0L,
    classification = "root",
    prior = 1.0
  )
  for (i in seq_along(TAX_RANKS)) {
    rankname <- rlang::sym(TAX_RANKS[i])
    tax_i <- dplyr::select(tax, one_of(TAX_RANKS[1:i]), n)
    tax_i <- tax_i[!is.na(tax_i[[TAX_RANKS[i]]]),]
    tax_i <- dplyr::summarize(tax_i, n = sum(n), .by = one_of(TAX_RANKS[1:i]))
    if (i == 1) {
      tax_i$parent_classification <- "root"
    } else {
      tax_i$parent_classification <-
        apply(tax_i[seq_len(i - 1)], 1, paste, collapse = ",", simplify = TRUE)
    }
    tax_i <- dplyr::left_join(
      tax_i,
      dplyr::select(
        tax_out[[i]],
        parent_classification = classification,
        parent_id = taxon_id,
        prior
      ),
      by = "parent_classification"
    )
    # todo: allow to specify priors
    # todo: unknown species?
    tax_i <- dplyr::mutate(tax_i, prior = prior * n / sum(n), .by = parent_id)
    if (i > 1L) {
      tax_i <- dplyr::mutate(
        tax_i,
        classification = paste(parent_classification, {{rankname}}, sep = ",")
      )
    } else {
      tax_i <- dplyr::mutate(
        tax_i,
        classification = {{rankname}}
      ) |>
        dplyr::arrange(classification == "nonFungi")
    }
    tax_i$taxon_id <- seq_len(nrow(tax_i)) + max(tax_out[[i]]$taxon_id)
    tax_i$rank <- i
    tax_out[[i + 1]] <-
      dplyr::select(tax_i, taxon_id, parent_id, rank, classification, prior)
  }
  dplyr::bind_rows(tax_out)
}

# Truncate classification(s) at a given rank
truncate_taxonomy <- function(s, rank) {
  regex <- paste0("(^([^,]+,){", rank-1L, "}[^,]+).*")
  out <- gsub(regex, "\\1", s)
  out[!grepl(regex, s)] <- NA_character_
  out
}

# Remove mycobank numbers from genus and species names
remove_mycobank_number <- function(taxon) {
  dplyr::if_else(startsWith(taxon, "pseudo"), taxon, sub("_[0-9]+$", "", taxon))
}


# Find OTUs which contain sequences with _any_ probability
# of being a target taxon
find_target_taxa <- function(target_taxa, asv_all_tax_prob, asv_taxonomy, otu_taxonomy) {
  asv_otu_key <-
    dplyr::inner_join(
      dplyr::select(asv_taxonomy, asv_seq_id = seq_id, {{TIP_RANK_VAR}}),
      dplyr::select(otu_taxonomy, seq_id, {{TIP_RANK_VAR}}),
      by = TIP_RANK
    ) %>%
    dplyr::select(-{{TIP_RANK_VAR}})
  otu_long_taxonomy <- tidyr::pivot_longer(
    otu_taxonomy,
    all_of(TAX_RANKS),
    names_to = "rank",
    values_to = "otu_taxon",
    names_transform = list(rank = rank2factor)
  ) %>%
    dplyr::select(seq_id, rank, otu_taxon)
  dplyr::select(
    asv_all_tax_prob,
    asv_seq_id = seq_id,
    rank,
    protax_taxon = taxon,
    protax_prob = prob
  ) %>%
    dplyr::inner_join(asv_otu_key, by = "asv_seq_id") %>%
    dplyr::group_by(seq_id) %>%
    dplyr::filter(any(protax_taxon %in% target_taxa)) %>%
    dplyr::left_join(otu_long_taxonomy, by = c("seq_id", "rank")) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(seq_id, asv_seq_id, rank) %>%
    dplyr::select(seq_id, asv_seq_id, otu_taxon, rank, everything())
}

summarize_uncross <- function(uncross) {
  uncross |>
  dplyr::summarize(
    Total_reads = sum(nread),
    Number_of_TagJump_Events = sum(is_tag_jump),
    TagJump_reads = sum(nread[is_tag_jump], na.rm = TRUE),
    ReadPercent_removed <- TagJump_reads / Total_reads * 100,
    .by = sample
  )
}

detect_numts <- function(a2m, id_is_int = FALSE) {
  checkmate::assert_file(a2m, access = "r")
  checkmate::assert_flag(id_is_int)
  a2m <- Biostrings::readBStringSet(a2m) |>
    sub(pattern = "^[acgt]+", replacement = "") |>
    sub(pattern = "[acgt]+$", replacement = "")

  a2m_frameshifts <- gregexpr("-+|[actg]+", a2m) |>
    purrr::map_dfr(
      ~data.frame(pos = .x, len = attr(.x, "match.length")),
      .id = "i"
    ) |>
    dplyr::filter(pos != 1L, pos + len != 659L, len %% 3 != 0)
  a2m_stops <- gregexpr("TA[GA]", a2m) |>
    purrr::map_dfr(
      ~data.frame(pos = .x, len = attr(.x, "match.length")),
      .id = "i"
    ) |>
    dplyr::filter(pos %% 3 == 2, pos > 0)
  numts = dplyr::bind_rows(
    frameshift = a2m_frameshifts,
    stop_codon = a2m_stops,
    .id = "numt_indicator"
  ) |>
    dplyr::mutate(seq_id = names(a2m)[as.integer(i)], .keep = "unused")

  if (id_is_int) {
    dplyr::transmute(
      numts,
      seq_idx = as.integer(seq_id),
      numt_indicator,
      pos,
      len
    )
  } else (
    dplyr::select(numts, seq_id, numt_indicator, pos, len)
  )
}

# convert a list of data to the XML format to be sent to KronaTools
xml_format <- function(data_format) {
  lapply(data_format, vapply, sprintf, "", fmt = "<val>{%s}</val>") |>
    vapply(paste, "", collapse = ",") |>
    purrr::imap_chr(sprintf, fmt="<%2$s>%1$s</%2$s>") |>
    paste(collapse = "\n")
}

krona_xml_nodes <- function(
    data,
    .rank,
    maxrank = rank2factor(TIP_RANK),
    outfile,
    pre = NULL,
    post = NULL,
    taxonomy = paste(KNOWN_TAXA, collapse = ","),
    node_data_format = NULL,
    node_xml_format = xml_format(node_data_format),
    ...
) {
  if (is.character(.rank)) .rank <- rank2factor(.rank)
  con <- outfile
  if (!methods::is(con, "connection")) {
    con <- file(con, open = "w")
    on.exit(close(con))
  }
  my_data <- data
  if (!is.null(taxonomy)) {
    my_data <- dplyr::filter(data, startsWith(parent_taxonomy, taxonomy))
  }
  xml <- dplyr::filter(my_data, rank == .rank) |>
    dplyr::transmute(
      taxon = taxon,
      taxonomy = ifelse(is.na(parent_taxonomy), taxon, paste(parent_taxonomy, taxon, sep = ",")),
      pre = glue::glue(
        '<node name="{taxon}">',
        node_xml_format,
        .sep = "\n"
      ),
      post = "</node>"
    )
  if (!is.null(pre)) {
    writeLines(pre, con)
  }
  if (.rank == maxrank) {
    writeLines(paste(xml$pre, xml$post, sep = "\n"), con)
  } else {
    purrr::pwalk(
      xml,
      krona_xml_nodes,
      data = my_data,
      .rank = subranks(.rank)[1],
      maxrank = maxrank,
      outfile = con,
      ...,
      node_xml_format = node_xml_format
    )
  }
  if (!is.null(post)) writeLines(post, con)
  outfile
}

read_sfile <- function(file) {
  tibble::tibble(
    text = readLines(file),
    is_widths = grepl("^#[- ]+$", text),
    part = cumsum(is_widths)
  ) |>
    dplyr::filter(is_widths | !startsWith(text, "#")) |>
    dplyr::group_split(part, .keep = FALSE) |>
    purrr::discard(\(x) nrow(x) == 1) |>
    purrr::map_dfr(
      \(x) {
        paste(x$text, collapse = "\n") |>
        readr::read_fwf(
          col_positions =  stringr::str_locate_all(x$text[1], "-+")[[1]] |>
            tibble::as_tibble() |>
            tibble::add_column(
              col_names = c("idx", "seq_id", "match_len", "cm_from", "cm_to",
                            "trunc", "bit_sc", "avg_pp", "time_band_calc",
                            "time_alignment", "time_total", "mem_mb")
            ) |>
            do.call(readr::fwf_positions, args = _),
          skip = 1,
          col_types = "iciiicdddddd"
        )
      }
    )
}

consensus_columns <- function(aln) {
  checkmate::assert_names(names(aln), must.include = c("alignment", "GC"))
  checkmate::assert_names(names(aln$GC), must.include = "RF")
  checkmate::assert_class(aln$alignment, "MultipleAlignment")
  checkmate::assert_class(aln$GC$RF, "BString")
  dots <- gregexpr("[.]+", aln$GC$RF)[[1]]
  Biostrings::colmask(aln$alignment) <-
    IRanges::IRanges(start = dots, width = attr(dots, "match.length"))
  methods::as(aln$alignment, "DNAStringSet")
}

file_to_sample_key <- function(filename) {
  sub("_(fwd|rev)_R[12]_(filt|trim)\\.fastq\\.gz", "", basename(filename))
}

# remove potential tag-jump from DADA2 ASVs table
#   core by Vladimir Mikryukov,
#   edited for 'targets' by Sten Anslan
#   modified to match OptimOTU style by Brendan Furneaux
remove_tag_jumps <- function(seqtable, f, p, id_col = "seq") {
  checkmate::assert_data_frame(seqtable)
  checkmate::assert_names(names(seqtable), must.include = c(id_col, "sample", "nread"))
  ## Load ASV table
  cat("...Number of ASVs: ", dplyr::n_distinct(seqtable[[id_col]]), "\n")
  n <- dplyr::n_distinct(seqtable$sample)
  cat("...Number of samples: ", n, "\n")

  ## UNCROSS score (with original parameter - take a root from the exp in denominator, to make curves more steep)
  uncross_score <- function(x, N, n, f = 0.01, tmin = 0.1, p = 1){
    # x = ASV abundance in a sample
    # N = total ASV abundance
    # n = number of samples
    # f = expected cross-talk rate, e.g. 0.01
    # tmin = min score to be considered as cross-talk
    # p = power to rise the exponent (default, 1; use 1/2 or 1/3 to make cureves more stepp)

    z <- f * N / n               # Expected treshold
    sc <- 2 / (1 + exp(x/z)^p)   # t-score
    data.frame(uncross = sc, is_tag_jump = sc >= tmin)
  }

  ## Estimate total abundance of sequence per plate
  out <- seqtable |>
    dplyr::mutate(total = sum(nread, na.rm = TRUE), .by = dplyr::all_of(id_col)) |>
    dplyr::select(-dplyr::all_of(id_col))



  ## Esimate UNCROSS score
  out <- cbind(
    out,
    uncross_score(
      x = out$nread,
      N = out$total,
      n = n,
      f = as.numeric(f),
      p = as.numeric(p)
    )
  )
  cat("...Number of tag-jumps: ", sum(out$is_tag_jump, na.rm = TRUE), "\n")
  # fwrite(x = TJ, file = "TagJump_stats.txt", sep = "\t")

  ## Remove detected tag-jumps from the ASV table
  out
}
