#### Functions which call external software from R
# Brendan Furneaux 2022

find_executable <- function(executable) {
  checkmate::assert_character(executable)
  out <- Sys.getenv(executable)
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- Sys.getenv(toupper(executable))
  }
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- Sys.which(executable)
  }
  if (nchar(out) == 0 || !file.exists(out)) {
    out <- list.files(path = "bin", pattern = executable, recursive = TRUE, full.names = TRUE)
  }
  checkmate::assert_file_exists(out, access = "x", .var.name = executable)
  out
}

# try to find the vsearch executable
find_vsearch <- function() {
  find_executable("vsearch")
}

# try to find the cutadapt executable
find_cutadapt <- function() {
  find_executable("cutadapt")
}

# try to find the hmmalign executable
find_hmmalign <- function() {
  find_executable("hmmalign")
}

# try to find the hmmsearch executable
find_hmmsearch <- function() {
  find_executable("hmmsearch")
}

# try to find the nhmmer executable
find_nhmmer <- function() {
  find_executable("nhmmer")
}

#' "usearch_global" function of vsearch
#'
#' @param query (`data.frame`, `Biostrings::DNAStringSet`, or `character` vector) query
#' sequences
#' @param ref (`data.frame`, `Biostrings::DNAStringSet`, `character` vector, or
#' file name) reference sequences
#' @param threshold (`numeric` scalar) identity threshold, in range 0.0-1.0
#' @param global (`logical` flag) if `TRUE`, end gaps and internal gaps are
#' penalized equally.  Otherwise end gaps are not penalized.
#' @param ncpu (`integer` count) number of threads to use
#'
#' @return `tibble::tibble` with columns `seq_id`, `clust`, and `dist`, where `seq_id` is
#' the name of a sequence from `query`, `clust` is the closest match to that
#' sequence in `ref`, and `dist` is the distance between them
#' @export
vsearch_usearch_global <- function(query, ref, threshold, global = TRUE,
                                   ncpu = local_cpus(), id_is_int = FALSE) {
  checkmate::check_flag(id_is_int)
  if (is.character(query) && length(query) == 1 && file.exists(query)) {
    tquery <- query
  } else {
    tquery <- withr::local_tempfile(pattern = "query", fileext = ".fasta")
    write_sequence(query, tquery)
  }
  if (is.character(ref) && length(ref) == 1 && file.exists(ref)) {
    tref <- ref
  } else {
    tref <- withr::local_tempfile(pattern = "ref", fileext = ".fasta")
    write_sequence(ref, tref)
  }
  assertthat::assert_that(assertthat::is.flag(global))
  gap <- if (global) "1" else "1I/0E"
  uc = system(
    paste(
      find_vsearch(),
      "--usearch_global", tquery,
      "--db", tref,
      "--id", threshold,
      "--uc", "-",
      "--maxaccepts", "100",
      "--top_hits_only",
      "--threads", ncpu,
      "--gapopen", gap,
      "--gapext", gap,
      "--match", "1",
      "--mismatch", "-1",
      "| awk '$1==\"H\" {print $9,$10,$4}'"
    ),
    intern = TRUE
  )
  stopifnot(attr(uc, "status") == 0)
  if (length(uc) > 0) {
    readr::read_delim(
      I(uc),
      col_names = c(if (id_is_int) "seq_idx" else "seq_id", "cluster", "dist"),
      delim = " ",
      col_types = if (id_is_int) "icd" else"ccd"
    )
  } else if (id_is_int) {
    tibble::tibble(
      seq_idx = integer(),
      cluster = character(),
      dist = numeric()
    )
  } else {
    tibble::tibble(
      seq_id = character(),
      cluster = character(),
      dist = numeric()
    )
  }
}

vsearch_uchime_ref <- function(query, ref, ncpu = local_cpus(), id_only = FALSE, id_is_int = FALSE) {
  checkmate::assert_flag(id_is_int)
  if (checkmate::test_file_exists(query, "r")) {
    tquery <- query
  } else {
    tquery <- tempfile("query", fileext = ".fasta")
    on.exit(unlink(c(tquery), force = TRUE))
    write_sequence(query, tquery)
  }
  if (checkmate::test_file_exists(ref, "r")) {
    tref <- ref
  } else {
    tref <- tempfile("ref", fileext = ".fasta")
    on.exit(unlink(c(tref), force = TRUE), add = TRUE)
    write_sequence(ref, tref)
  }
  tchimeras <- tempfile("chimeras", fileext = ".fasta")
  on.exit(unlink(tchimeras), TRUE)
  vs <- system2(
    find_vsearch(),
    args = c(
      "--uchime_ref", tquery,
      "--db", tref,
      "--chimeras", tchimeras,
      "--threads", ncpu
    )
  )
  stopifnot(vs == 0L)
  if (id_only) {
    out <- names(Biostrings::fasta.seqlengths(tchimeras))
    if (id_is_int) {
      as.integer(out)
    } else {
      out
    }
  } else {
    out <- Biostrings::readDNAStringSet(tchimeras) %>%
    as.character() %>%
    tibble::enframe(name = "seq_id", value = "seq")
    if (id_is_int) {
      dplyr::transmute(out, seq_idx = as.integer(seq_id), seq)
    } else {
      out
    }
  }
}

vsearch_usearch_global_closed_ref <- function(query, ref, threshold, ...) {
  out <- tibble::tibble(seq_id = character(0), cluster = character(0))
  while(sequence_size(query) > 0 && sequence_size(ref) > 0) {
    result <- vsearch_usearch_global(query, ref, threshold, ...)
    if (nrow(out) > 0) {
      result <- dplyr::left_join(
        result,
        out,
        by = c("cluster" = "seq_id"),
        suffix = c(".orig", "")
      ) %>%
        dplyr::select(seq_id, cluster)
    }
    out <- dplyr::bind_rows(out, result)
    ref <- select_sequence(query, result$seq_id)
    query <- select_sequence(query, result$seq_id, negate = TRUE)
  }
  out
}

# build a usearch database (UDB) file using USEARCH

build_udb <- function(infile, outfile, type = c("usearch", "sintax", "ublast"),
                      usearch = Sys.which("usearch")) {
  type <- match.arg(type)
  command <- paste0("-makeudb_", type)
  args <- c(
    command, infile,
    "-output", outfile
  )
  result <- system2(usearch, args)
  stopifnot(result == 0)
  outfile
}

build_filtered_udb <- function(
  infile,
  outfile,
  type = c("usearch", "sintax", "ublast"),
  blacklist,
  usearch = Sys.which("usearch")
) {
  # make sure we have valid arguments
  type <- match.arg(type)
  command <- paste0("-makeudb_", type)
  stopifnot(system2(usearch, "--version")==0)

  # make a temp file and a temp fifo
  blf <- tempfile(fileext = ".txt")
  tf <- tempfile(fileext = ".fasta")
  on.exit(unlink(c(tf, blf), force = TRUE))
  writeLines(blacklist, blf)
  stopifnot(system2("mkfifo", tf) == 0)

  # first usearch call removes the blacklisted sequences
  system2(
    usearch,
    args = c(
      "--fastx_getseqs", infile,
      "--labels", blf,
      "--label_substr_match",
      "--notmatched", tf
    ),
    wait = FALSE
  )
  # second usearch call creates the udb file
  result = system2(
    usearch,
    args = c(
      command, tf,
      "--output", outfile
    )
  )
  stopifnot(result == 0)
  outfile
}

vsearch_cluster_smallmem <- function(seq, threshold = 1, ncpu = local_cpus()) {
  if (is.character(seq) && length(seq) == 1 && file.exists(seq)) {
    tout <- seq
  } else {
    tout <- withr::local_tempfile(pattern = "data", fileext = ".fasta")
    write_sequence(seq, tout)
  }
  tin <- tempfile("data", fileext = ".uc")
  uc = system(
    paste(
      find_vsearch(),
      "--cluster_smallmem", tout,
      "--usersort",
      "--id", threshold,
      "--uc -",
      "--threads", ncpu,
      "| awk '$1==\"H\" {print $9,$10}'"
    ),
    intern = TRUE
  )
  stopifnot(attr(uc, "status") == 0)
  if (length(uc) > 0) {
    readr::read_delim(
      I(uc),
      col_names = c("query", "hit"),
      delim = " ",
      col_types = "cc"
    )
  } else {
    tibble::tibble(query = character(), hit = character())
  }
}

collapseNoMismatch_vsearch <- function(seqtab, ..., ncpu = local_cpus()) {
  seqs <- colnames(seqtab)
  names(seqs) <- seq_along(seqs)
  matches <- vsearch_cluster_smallmem(seqs, ncpu = ncpu)
  map <- tibble::tibble(
    seq_idx_in = seq_len(ncol(seqtab)),
    seq_idx_out = seq_len(ncol(seqtab))
  )
  if (nrow(matches) > 0) {
    matches$query <- as.integer(matches$query)
    matches$hit <- as.integer(matches$hit)
    matches <- matches[order(matches$query),]
    for (i in unique(matches$hit)) {
      seqtab[,i] <- seqtab[,i] +
        as.integer(rowSums(seqtab[,matches$query[matches$hit == i], drop = FALSE]))
    }
    seqtab <- seqtab[,-matches$query]
    map$seq_idx_out[matches$query] <- matches$hit
    map$seq_idx_out = map$seq_idx_out - findInterval(map$seq_idx_out, matches$query)
  }
  attr(seqtab, "map") <- map
  return(seqtab)
}

nomismatch_hits_vsearch <- function(seqtab, seqs = NULL,
                                    abund_col = "nread",
                                    fastx_index = NULL,
                                    ...,
                                    ncpu = local_cpus()) {
  if (is.null(seqs)) {
    checkmate::assert_names(
      names(seqtab),
      must.include = c("seq", abund_col),
      disjunct.from = "seq_idx"
    )
    seqs <- sort_seq_table(seqtab, abund_col = abund_col, ...)
  } else {
    checkmate::assert(
      checkmate::check_names(
        names(seqtab),
        must.include = c("seq_idx", abund_col),
        disjunct.from = "seq"
      ),
      checkmate::check_names(
        names(seqtab),
        must.include = c("seq_id", abund_col),
        disjunct.from = "seq"
      )
    )
    o <- sort_seq_table(seqtab, seqs = seqs, abund_col = abund_col, ...)
    if (checkmate::check_file_exists(seqs)) {
      # no easy way to re-order without reading it all into memory
      seqs <- Biostrings::readDNAStringSet(seqs)
    }
    if (is.character(o)) o <- match(o, names(seqs))
    seqs <- seqs[o]
    names(seqs) <- as.character(o)
  }
  nseqs <- sequence_size(seqs)
  vsearch_cluster_smallmem(seqs, ncpu = ncpu) |>
    dplyr::mutate(dplyr::across(everything(), as.integer))
}

deduplicate_seqtable <- function(seqtable, hits, abund_col = "nread",
                                 sample_cols = "sample", merge = TRUE) {
  checkmate::assert_character(abund_col)
  checkmate::assert_character(sample_cols)
  checkmate::assert_data_frame(seqtable)
  checkmate::assert_names(
    names(seqtable),
    must.include = c("seq_idx", abund_col, sample_cols)
  )
  checkmate::assert_data_frame(hits)
  checkmate::assert_names(
    names(hits),
    must.include = c("query", "hit")
  )
  checkmate::assert_integer(hits$query, lower = 0L, any.missing = FALSE)
  checkmate::assert_integer(hits$hit, lower = 0L, any.missing = FALSE)
  checkmate::assert_integer(seqtable$seq_idx, lower = 0L, any.missing = FALSE)
  checkmate::assert_flag(merge)

  hits <- dplyr::arrange(hits, query)
  seqtable <- dplyr::left_join(
    seqtable,
    hits,
    by = c("seq_idx" = "query")
  ) |>
    dplyr::mutate(
      seq_idx = dplyr::coalesce(hit, seq_idx),
      .keep = "unused"
    )
  if (isTRUE(merge)) {
    seqtable <- dplyr::summarize(
      seqtable,
      dplyr::across(all_of(abund_col), sum),
      .by = any_of(c("seq_idx", sample_cols))
    )
  }
  seqtable$seq_idx <- seqtable$seq_idx - findInterval(seqtable$seq_idx, hits$query)
  seqtable
}

# TODO: allow sequence list which is not a file, not named with integers, etc.
# TODO: do this without reading the whole file into memory
deduplicate_seqs <- function(seqs, hits, outfile) {
  checkmate::assert_file_exists(seqs, "r")
  checkmate::assert_data_frame(hits)
  checkmate::assert_names(
    names(hits),
    must.include = c("query", "hit")
  )
  checkmate::assert_integer(hits$query, lower = 0L, any.missing = FALSE)
  checkmate::assert_integer(hits$hit, lower = 0L, any.missing = FALSE)
  if (nrow(hits) > 0) {
    out <- Biostrings::readBStringSet(seqs)[-hits$query]
    names(out) <- as.character(seq_along(out))
    write_sequence(out, outfile, compress = endsWith(outfile, ".gz"), compression_level = 9)
  } else {
    file.copy(seqs, outfile, overwrite = TRUE)
    outfile
  }
}

deduplicate_seq_idx <- function(seq_idx, hits, merge = TRUE) {
  checkmate::assert_integerish(seq_idx, lower = 1)
  deduplicate_seqtable(
    seqtable = tibble::tibble(
      seq_idx = seq_idx
    ),
    hits = hits,
    sample_cols = character(),
    abund_col = character(),
    merge = merge
  )$seq_idx
}

cutadapt_paired_option_names <- c(
  "max_err",
  "min_overlap",
  "action",
  "discard_untrimmed",
  "max_n",
  "max_ee",
  "min_length",
  "max_length",
  "truncQ_R1",
  "truncQ_R2",
  "cut_R1",
  "cut_R2"
)

cutadapt_paired_options <- function(
    max_err = 0.2,
    min_overlap = 10L,
    action = "retain",
    discard_untrimmed = TRUE,
    max_n = 0L,
    max_ee = NULL,
    min_length = NULL, max_length = NULL,
    truncQ_R1 = 2, truncQ_R2 = 2,
    cut_R1 = NULL, cut_R2 = NULL
) {
  checkmate::assert_number(max_err, lower = 0, null.ok = TRUE)
  checkmate::assert_count(min_overlap, positive = TRUE, null.ok = TRUE)
  checkmate::assert_choice(
    action,
    c("trim", "retain", "mask", "lowercase", "none"),
    null.ok = TRUE
  )
  checkmate::assert_flag(discard_untrimmed)
  checkmate::assert_count(max_n, null.ok = TRUE)
  checkmate::assert_number(max_ee, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_count(min_length, positive = TRUE, null.ok = TRUE)
  checkmate::assert_integerish(truncQ_R1, min.len = 1, max.len = 2, null.ok = TRUE)
  checkmate::assert_integerish(truncQ_R2, min.len = 1, max.len = 2, null.ok = TRUE)
  checkmate::assert_integerish(cut_R1, min.len = 1, max.len = 2, null.ok = TRUE)
  checkmate::assert_integerish(cut_R2, min.len = 1, max.len = 2, null.ok = TRUE)
  structure(
    list(
      max_err = max_err,
      min_overlap = min_overlap,
      action = action,
      discard_untrimmed = discard_untrimmed,
      max_n = max_n,
      max_ee = max_ee,
      min_length = min_length,
      truncQ_R1 = truncQ_R1,
      truncQ_R2 = truncQ_R2,
      cut_R1 = cut_R1,
      cut_R2 = cut_R2
    ),
    class = "cutadapt_paired_options"
  )
}

cutadapt_paired_filter_trim <- function(
  file_R1, file_R2,
  primer_R1, primer_R2,
  trim_R1, trim_R2,
  options = cutadapt_paired_options(),
  ncpu = local_cpus(),
  cutadapt = find_cutadapt(),
  logfile = NULL,
  ...
) {
  checkmate::assert_class(options, "cutadapt_paired_options")
  args <- c(
    "-g", primer_R1,
    "-G", primer_R2,
    "-o", trim_R1,
    "-p", trim_R2
  )
  if (!is.null(options$max_err)) {
    args <- c(args, "-e", options$max_err)
  }
  if (!is.null(options$min_overlap)) {
    args <- c(args, "-O", options$min_overlap)
  }
  if (!is.null(options$action)) {
    args <- c(args, paste0("--action=", options$action))
  }
  if (isTRUE(options$discard_untrimmed)) {
    args <- c(args, "--discard-untrimmed")
  }
  if (!is.null(options$max_n)) {
    args <- c(args, "--max-n", options$max_n)
  }
  if (!is.null(options$max_ee)) {
    args <- c(args, "--max-ee", options$max_ee)
  }
  if (!is.null(options$min_length)) {
    args <- c(args, "-m", options$min_length)
  }
  if (!is.null(options$max_length)) {
    args <- c(args, "-M", options$max_length)
  }
  if (!is.null(options$truncQ_R1)) {
    args <- c(args, "-q", paste(round(options$truncQ_R1), collapse = ","))
  }
  if (!is.null(options$truncQ_R2)) {
    args <- c(args, "-Q", paste(round(options$truncQ_R2), collapse = ","))
  }
  if (!is.null(options$cut_R1)) {
    args <- c(args, "-u", paste(round(options$cut_R1), collapse = ","))
  }
  if (!is.null(options$cut_R2)) {
    args <- c(args, "-U", paste(round(options$cut_R2), collapse = ","))
  }
  if (!is.null(ncpu)) {
    assertthat::assert_that(assertthat::is.count(ncpu))
    args <- c(args, "-j", ncpu)
  }
  args <- c(args, file_R1, file_R2)
  out <- processx::run(
    cutadapt,
    args = args,
    error_on_status = TRUE
  )
  if (!is.null(logfile)) writeLines(out$stdout, logfile)
  c(trim_R1, trim_R2)
}

cutadapt_option_names <- c(
  "max_err",
  "min_overlap",
  "action",
  "discard_untrimmed",
  "max_n",
  "max_ee",
  "min_length",
  "max_length",
  "truncQ",
  "cut"
)

cutadapt_options <- function(
    max_err = 0.2,
    min_overlap = 10L,
    action = "retain",
    discard_untrimmed = TRUE,
    max_n = 0L,
    max_ee = NULL,
    min_length = NULL, max_length = NULL,
    truncQ = NULL,
    cut = NULL
) {
  checkmate::assert_number(max_err, lower = 0, null.ok = TRUE)
  checkmate::assert_count(min_overlap, positive = TRUE, null.ok = TRUE)
  checkmate::assert_choice(
    action,
    c("trim", "retain", "mask", "lowercase", "none"),
    null.ok = TRUE
  )
  checkmate::assert_flag(discard_untrimmed)
  checkmate::assert_count(max_n, null.ok = TRUE)
  checkmate::assert_number(max_ee, lower = 0, finite = TRUE, null.ok = TRUE)
  checkmate::assert_count(min_length, positive = TRUE, null.ok = TRUE)
  checkmate::assert_count(max_length, positive = TRUE, null.ok = TRUE)
  checkmate::assert_integerish(truncQ, min.len = 1, max.len = 2, null.ok = TRUE)
  checkmate::assert_integerish(cut, min.len = 1, max.len = 2, null.ok = TRUE)
  structure(
    list(
      max_err = max_err,
      min_overlap = min_overlap,
      action = action,
      discard_untrimmed = discard_untrimmed,
      max_n = max_n,
      max_ee = max_ee,
      min_length = min_length,
      max_length = max_length,
      truncQ = truncQ,
      cut = cut
    ),
    class = "cutadapt_options"
  )
}

# replace existing (default?) options with values from a data.frame or list.
# if a data.frame, then the values in the data.frame should be all the same.

update.cutadapt_paired_options <- function(options, new_options) {
  checkmate::assert(
    checkmate::check_list(new_options, null.ok = TRUE),
    checkmate::check_data_frame(new_options, null.ok = TRUE),
    checkmate::check_character(new_options, null.ok = TRUE)
  )
  new_options <-
    new_options[intersect(names(new_options), cutadapt_paired_option_names)]
  if (is.data.frame(new_options) && ncol(new_options) > 0) {
    new_options <- unique(new_options)
    if (nrow(new_options) > 1L)
      stop(
        "'new_options' must be the same for all samples in a batch. \n",
        "Current batch has ", nrow(new_options),
        "unique combinations of options."
      )
  }
  new_options <- lapply(unclass(new_options), unlist)
  if (length(new_options) > 0) {
    options[names(new_options)] <- new_options
    do.call(cutadapt_paired_options, options)
  } else {
    options
  }
}

update.cutadapt_options <- function(options, new_options) {
  checkmate::assert(
    checkmate::check_list(new_options, null.ok = TRUE),
    checkmate::check_data_frame(new_options, null.ok = TRUE),
    checkmate::check_character(new_options, null.ok = TRUE)
  )
  new_options <-
    new_options[intersect(names(new_options), cutadapt_option_names)]
  if (is.data.frame(new_options) && ncol(new_options) > 0) {
    new_options <- unique(new_options)
    if (nrow(new_options) > 1L)
      stop(
        "'new_options' must be the same for all samples in a batch. \n",
        "Current batch has ", nrow(new_options),
        "unique combinations of options."
      )
  }
  new_options <- lapply(unclass(new_options), unlist)
  if (length(new_options) > 0) {
    options[names(new_options)] <- new_options
    do.call(cutadapt_options, options)
  } else {
    options
  }
}

cutadapt_filter_trim <- function(
  file,
  primer,
  trim,
  options = cutadapt_options(),
  ncpu = local_cpus(),
  cutadapt = find_cutadapt(),
  logfile = NULL,
  ...
) {
  args <- c(
    "-g", primer,
    "-o", trim
  )
  if (!is.null(options$max_err)) {
    args <- c(args, "-e", options$max_err)
  }
  if (!is.null(options$min_overlap)) {
    args <- c(args, "-O", round(options$min_overlap))
  }
  if (!is.null(options$action)) {
    args <- c(args, paste0("--action=", options$action))
  }
  if (isTRUE(options$discard_untrimmed)) {
    args <- c(args, "--discard-untrimmed")
  }
  if (!is.null(options$max_n)) {
    args <- c(args, "--max-n", options$max_n)
  }
  if (!is.null(options$max_ee)) {
    args <- c(args, "--max-ee", options$max_ee)
  }
  if (!is.null(options$min_length)) {
    args <- c(args, "-m", options$min_length)
  }
  if (!is.null(options$max_length)) {
    args <- c(args, "-M", options$max_length)
  }
  if (!is.null(options$truncQ)) {
    args <- c(args, "-q", paste(round(options$truncQ), collapse = ","))
  }
  if (!is.null(options$cut)) {
    args <- c(args, "-u", paste(round(options$cut), collapse = ","))
  }
  if (!is.null(ncpu)) {
    checkmate::assert_count(ncpu, positive = TRUE)
    args <- c(args, "-j", ncpu)
  }
  args <- c(args, file)
  out <- processx::run(
    cutadapt,
    args = args,
    error_on_status = TRUE
  )
  if (!is.null(logfile)) writeLines(out$stdout, logfile)
  trim
}

trim_primer <- function(seqs, primer, ...) {
  tempseqs <- tempfile(fileext = ".fasta")
  write_sequence(seqs, tempseqs)
  temptrimmed <- tempfile(fileext = ".fasta")
  on.exit(unlink(c(tempseqs, temptrimmed), force = TRUE))
  cutadapt_filter_trim(
    file = tempseqs,
    primer = primer,
    trim = temptrimmed,
    ...
  )
  Biostrings::readDNAStringSet(temptrimmed) %>%
    as.character() %>%
    tibble::enframe(name = "seq_id", value = "seq")
}


hmmalign <- function(seqs, hmm, outfile, outformat = "A2M",
                     compress = endsWith(outfile, ".gz")) {
  checkmate::assert_string(hmm)
  checkmate::assert_file_exists(hmm, access = "r")
  ensure_directory(outfile)
  checkmate::assert_path_for_output(outfile, overwrite = TRUE)
  checkmate::assert_choice(outformat, c("A2M", "a2m", "afa", "AFA"))
  checkmate::assert_flag(compress)
  exec <- find_hmmalign()
  checkmate::assert_file_exists(exec, access = "x")
  if (checkmate::test_file_exists(seqs, "r")) {
    tseqs <- seqs
    n <- length(seqs)
  } else if (checkmate::test_list(seqs, types = c("character", "XStringSet", "data.frame"))) {
    n <- length(seqs)
    tseqs <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
    purrr::pwalk(list(seq = seqs, fname = tseqs), write_sequence)
  } else {
    checkmate::assert_multi_class(seqs, c("data.frame", "character", "XStringSet"))
    n <- 1
    tseqs <- withr::local_tempfile(fileext = ".fasta")
    write_sequence(seqs, tseqs)
  }
  checkmate::assert(
    length(outfile) == 1,
    length(outfile) == n
  )
  if (length(outfile) == 1 && n > 1) {
    tout <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
  } else {
    tout <- outfile
  }
  if (compress) {
    if (identical(tout, outfile)) {
      tout <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
    }
  }
  mout <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
  for (i in seq_len(n)) {
    processx::run("mkfifo", mout[i])
  }

  args <- data.frame(
    "--outformat", outformat,
    "--trim",
    "-o", mout,
    hmm,
    tseqs
  )
  args <- as.matrix(args)
  hmmer <- vector("list", n)
  deline <- vector("list", n)
  for (i in seq_len(n)) {
    hmmer[[i]] <- processx::process$new(
      command = exec,
      args = args[i,],
      supervise = TRUE
    )
    deline[[i]] <- processx::process$new(
      command = "awk",
      args = 'BEGIN{ORS=""};NR>1&&/^>/{print "\\n"};{print};/^>/{print "\\n"};END{print "\\n"}',
      stdin = mout[i],
      stdout = tout[i],
      supervise = TRUE
    )
  }
  hmmer_return <- integer()
  for (i in seq_len(n)) {
    hmmer[[i]]$wait()
    hmmer_return <- union(hmmer[[i]]$get_exit_status(), hmmer_return)
  }
  stopifnot(identical(hmmer_return, 0L))
  for (i in seq_len(n)) {
    deline[[i]]$wait()
  }

  if (compress && length(outfile) == n) {
    gzip <- vector("list", n)
    for (i in seq_len(n)) {
      gzip[[i]] <- processx::process$new(
        command = "gzip",
        error_on_status = TRUE,
        args = c("-c", tout[i]),
        stdout = outfile[i]
      )
    }
    gzip_return <- integer()
    for (i in seq_len(n)) {
      gzip_return <- union(gzip[[i]]$wait()$status, gzip_return)
    }
    stopifnot(identical(gzip_return, 0L))
  } else if (length(outfile) < n) {
    fastx_combine(tout, outfile)
  }
  outfile
}

read_hmmer_tblout <- function(file, col_names, col_types) {
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
            col_positions =  stringr::str_locate_all(x$text[1], "#?-+")[[1]] |>
              tibble::as_tibble() |>
              tibble::add_column(col_names = col_names) |>
              do.call(readr::fwf_positions, args = _),
            skip = 1,
            col_types = col_types
          )
      }
    )

}

read_domtblout <- function(file) {
  read_hmmer_tblout(
    file,
    col_names = c("seq_name", "seq_accno", "seq_length", "hmm_name",
                  "hmm_accno", "hmm_length", "Evalue", "full_score",
                  "full_bias", "hit_num", "total_hits", "c_Evalue",
                  "i_Evalue", "hit_score", "hit_bias", "hmm_from", "hmm_to",
                  "seq_from", "seq_to", "env_from", "env_to", "acc",
                  "description"),
    col_types = "cciccidddiiddddiiiiiidc"
  )
}

read_dna_tblout <- function(file) {
  read_hmmer_tblout(
    file,
    col_names = c(
      "seq_name", "seq_accno", "hmm_name", "hmm_accno",
      "hmm_from", "hmm_to", "seq_from", "seq_to", "env_from", "env_to",
      "seq_len", "strand", "Evalue", "bit_score", "bias", "description"
    ),
    col_types = "cccciiiiiiicnnnc"
  )

}

hmmsearch <- function(seqs, hmm) {
  checkmate::assert_string(hmm)
  checkmate::assert_file_exists(hmm, access = "r")
  exec <- find_hmmsearch()
  checkmate::assert_file_exists(exec, access = "x")
  if (checkmate::test_file_exists(seqs, "r")) {
    tseqs <- seqs
    n <- length(seqs)
  } else if (checkmate::test_list(seqs, types = c("character", "XStringSet", "data.frame"))) {
    n <- length(seqs)
    tseqs <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
    purrr::pwalk(list(seq = seqs, fname = tseqs), write_sequence)
  } else {
    checkmate::assert_multi_class(seqs, c("data.frame", "character", "XStringSet"))
    n <- 1
    tseqs <- withr::local_tempfile(fileext = ".fasta")
    write_sequence(seqs, tseqs)
  }
  outfile <- replicate(n, tempfile(fileext = ".hmmout"))
  args <- data.frame(
    "--noali",
    "--notextw",
    "--domtblout", outfile,
    hmm,
    tseqs
  )
  args <- as.matrix(args)

  hmmer <- vector("list", n)
  for (i in seq_len(n)) {
    hmmer[[i]] <- processx::process$new(
      command = exec,
      args = args[i,],
      supervise = TRUE
    )
  }
  hmmer_return <- integer()
  for (i in seq_len(n)) {
    hmmer[[i]]$wait()
    hmmer_return <- union(hmmer[[i]]$get_exit_status(), hmmer_return)
    stopifnot(identical(hmmer_return, 0L))
  }
  purrr::map_dfr(
    outfile,
    read_domtblout
  )
}

nhmmer <- function(seqs, hmm, ncpu = local_cpus()) {
  checkmate::assert_string(hmm)
  checkmate::assert_file_exists(hmm, access = "r")
  checkmate::assert_count(ncpu)
  exec <- find_nhmmer()
  checkmate::assert_file_exists(exec, access = "x")
  if (length(seqs) == 1 && checkmate::test_file_exists(seqs, "r")) {
    tseqs <- seqs
  } else {
    checkmate::assert_multi_class(seqs, c("data.frame", "character", "XStringSet"))
    tseqs <- withr::local_tempfile(fileext = ".fasta")
    write_sequence(seqs, tseqs)
  }
  outfile <- withr::local_tempfile(fileext = ".hmmout")
  args <- c(
    "--noali",
    "--notextw",
    "--tblout", outfile,
    "--watson",
    "--cpu", ncpu,
    hmm,
    tseqs
  )
  processx::run(
      command = exec,
      args = args,
      error_on_status = TRUE
  )
  read_dna_tblout(outfile)
}

run_protax <- function(seqs, outdir, modeldir, ncpu = local_cpus()) {
  if (dir.exists(outdir)) unlink(outdir, recursive = TRUE)
  dir.create(outdir)
  if (length(seqs) == 1 && file.exists(seqs)) {
    if (seqs != file.path(outdir, "all.fa"))
      file.copy(seqs, file.path(outdir, "all.fa"))
  } else {
    write_sequence(seqs, file.path(outdir, "all.fa"))
  }
  status <- system2(
    "scripts/runprotax",
    c(outdir, modeldir, ncpu)
  )
  stopifnot(status == 0)
  list.files(outdir, full.names = TRUE)
}

# vectorized on aln_seqs
# attempts to run them ALL in parallel, be careful!
run_protax_animal <- function(aln_seqs, modeldir, min_p = 0.1, rep_p = 0.01,
                              strip_inserts = FALSE, id_is_int = FALSE,
                              info = FALSE, options = character()) {
  checkmate::assert_file_exists(aln_seqs, access = "r")
  checkmate::assert_directory_exists(modeldir)
  priors <- file.path(modeldir, "taxonomy.priors")
  checkmate::assert_file_exists(priors, access = "r")
  refs <- file.path(modeldir, "refs.aln")
  checkmate::assert_file_exists(refs, access = "r")
  rseqs <- file.path(modeldir, "model.rseqs.numeric")
  checkmate::assert_file_exists(rseqs, access = "r")
  pars <- file.path(modeldir, "model.pars")
  checkmate::assert_file_exists(pars, access = "r")
  scs <- file.path(modeldir, "model.scs")
  checkmate::assert_file_exists(scs, access = "r")
  checkmate::assert_flag(info)
  executable <- find_executable(if (info) "classify_info" else "classify_v2")
  checkmate::check_number(min_p, lower = 0, upper = 1, finite = TRUE)
  checkmate::check_number(rep_p, lower = 0, upper = min_p, finite = TRUE)
  checkmate::assert_flag(strip_inserts)
  checkmate::assert_flag(id_is_int)
  checkmate::assert_character(options)
  args <- c("-t", rep_p, options, priors, refs, rseqs, pars, scs, as.character(min_p))

  is_gz <-endsWith(aln_seqs, ".gz")
  stopifnot(all(is_gz) | all(!is_gz))
  is_gz <- all(is_gz)

  n <- length(aln_seqs)
  if (strip_inserts | is_gz) {
    prepipe <- vector("list", n)
    protax_in <- replicate(n, withr::local_tempfile(fileext = ".fasta"))
    pipecommand <- if (is_gz) "zcat" else "cat"
    pipecommand <- paste(pipecommand, aln_seqs)
    if (strip_inserts) pipecommand <- paste(pipecommand, "| tr -d 'acgt'")
    pipecommand <- paste(pipecommand, ">", protax_in)
  } else {
    protax_in <- aln_seqs
  }

  protax <- vector("list", n)
  outfiles <- replicate(n, withr::local_tempfile())
  for (i in seq_len(n)) {
    if (strip_inserts | is_gz) {
      system(pipecommand[i])
    }
    protax[[i]] <- processx::process$new(
      command = executable,
      args = c(args, protax_in[i]),
      stdout = outfiles[i]
    )
  }
  protax_exit_status = 0L
  output <- vector("list", n)
  for (i in seq_len(n)) {
    protax[[i]]$wait()
    protax_exit_status <- max(protax_exit_status, protax[[i]]$get_exit_status())
    stopifnot(protax_exit_status == 0L)
    output[[i]] <- readr::read_delim(
      outfiles[i],
      col_names = c(
        if (id_is_int) "seq_idx" else "seq_id",
        "rank",
        "taxonomy",
        "prob",
        if (info) c("best_id", "best_dist", "second_id", "second_dist") else NULL
      ),
      col_types = paste0(
        if (id_is_int) "i" else "c",
        "icn",
        if (info) "-in-in-" else ""
      )
    )
  }
  dplyr::bind_rows(output)
}

run_protax_besthit <- function(aln_query, aln_ref, options = character(),
                               query_id_is_int = TRUE, ref_id_is_int = TRUE,
                               command = "dist_best") {
  checkmate::assert_file_exists(aln_query, access = "r")
  n_query <- length(aln_query)
  checkmate::assert_file_exists(aln_ref, access = "r")
  n_ref <- length(aln_ref)
  executable <- find_executable(command)
  checkmate::assert_character(options)
  checkmate::assert_flag(query_id_is_int)
  checkmate::assert_flag(ref_id_is_int)

  is_gz_ref <- endsWith(aln_ref, ".gz")
  stopifnot(all(is_gz_ref) | all(!is_gz_ref))
  is_gz_ref <- all(is_gz_ref)

  if (n_ref > 1) {
    ref <- replicate(n_query, withr::local_tempfile(fileext = ".fasta"))
    refcommand <- if (is_gz_ref) "zcat" else "cat"
    for (i in seq_len(n_query)) {
      system2("mkfifo", args = ref[i])
      system2(
        command = refcommand,
        args = aln_ref,
        stdout = ref[i],
        wait = FALSE
      )
    }
  } else {
    ref <- rep_len(aln_ref, n_query)
  }

  besthit <- vector("list", n_query)
  outfiles <- replicate(n_query, withr::local_tempfile())
  for (i in seq_len(n_query)) {
    besthit[[i]] <- processx::process$new(
      command = executable,
      args = c(options, ref[i], aln_query[i]),
      stdout = outfiles[i]
    )
  }
  besthit_exit_status = 0L
  output <- vector("list", n_query)
  for (i in seq_len(n_query)) {
    besthit[[i]]$wait()
    besthit_exit_status <- max(besthit_exit_status, besthit[[i]]$get_exit_status())
    stopifnot(besthit_exit_status == 0L)
    besthit[[i]] <- readr::read_delim(
      outfiles[i],
      col_names = c(
        if (query_id_is_int) "seq_idx" else "seq_id",
        if (ref_id_is_int) "ref_idx" else "ref_id",
        "dist"
      ),
      col_types = paste0(
        if (query_id_is_int) "i" else "c",
        if (ref_id_is_int) "i" else "c",
        "d"
      ),
      delim = " "
    )
  }
  dplyr::bind_rows(besthit)
}

run_protax_bipart <- function(aln_query, aln_ref, max_d = 0.2,
                              query_id_is_int = TRUE, ref_id_is_int = TRUE) {
  checkmate::assert_number(max_d, lower = 0, upper = 1)
  run_protax_besthit(
    aln_query = aln_ref,
    aln_ref = aln_query,
    options = as.character(max_d),
    query_id_is_int = query_id_is_int,
    ref_id_is_int = ref_id_is_int,
    command = "dist_bipart"
  )
}

protax_besthit_closedref <- function(infile, index, i, unknowns, thresh,
                                     seq_width, ncpu = local_cpus(),
                                     max_gap = 100L) {
  checkmate::assert_file(infile, "r")
  checkmate::assert_file(index, "r")
  checkmate::assert_integerish(i)
  checkmate::assert_logical(unknowns, len = length(i))

  seqs <- fastx_gz_random_access_extract(
    infile = infile,
    index = index,
    i = i,
    ncpu = ncpu,
    max_gap = max_gap
  )

  queries <- names(seqs)[unknowns]
  refs <- names(seqs)[!unknowns]

  out <- NULL
  while (length(queries) > 0 && length(refs) > 0) {
    query_file <- write_sequence(
      seqs[queries],
      withr::local_tempfile(fileext=".fasta")
    ) |>
      fastx_split(
        # only parallelize if it is likely to be worth it.
        n = if (length(queries) > 1000) local_cpus() else 1L,
        outroot = tempfile(tmpdir = withr::local_tempdir()),
        compress = FALSE
      )

    ref_file <- write_sequence(
      seqs[refs],
      withr::local_tempfile(fileext=".fasta")
    )

    newout <-
      run_protax_besthit(
        aln_query = query_file,
        aln_ref = ref_file,
        options = c(
          "-m", "100",
          "-l", seq_width,
          "-r", as.character(length(refs))
        ),
        query_id_is_int = FALSE,
        ref_id_is_int = FALSE
      ) |>
      dplyr::filter(dist <= 1 - thresh/100)

    refs <- unique(newout$seq_id)
    queries <- setdiff(queries, refs)
    if (is.null(out)) {
      out <- dplyr::rename(newout, cluster = ref_id)
    } else {
      newout <- dplyr::left_join(newout, out, by = c("ref_id" = "seq_id")) |>
        dplyr::select(-ref_id)
      out <- dplyr::bind_rows(out, newout)
    }
  }
  out
}

seq_cluster_protax <- function(aln_seq, aln_index, which, thresh, aln_len) {
  nslice <- floor(sqrt(local_cpus()-1))
  allseq <- fastx_gz_random_access_extract(
    infile = aln_seq,
    index = aln_index,
    i = which
  )
  names(allseq) <- as.character(seq_along(allseq) - 1L)
  if (length(which) > 1000) {
    seq <- character(nslice)
    spl <- sort(rep_len(seq_len(nslice), length(which)))
    i <- split(seq_len(length(which)), spl)
    for (j in seq_len(nslice)) {
      seq[j] <- write_sequence(
        allseq[i[[j]]],
        fname = withr::local_tempfile(fileext = ".fasta")
      )
    }
    i_half <- lapply(i, \(x) split(x, rep(c(1, 2), length.out = length(x))))
    i_half <- do.call(c, args = i_half)
    seq_half <- character(nslice*2)
    for (j in seq(3, nslice*2)) {
      seq_half[j] <- write_sequence(
        allseq[i_half[[j]]],
        fname = tempfile(tmpdir = withr::local_tempdir())
      )
    }
  } else {
    nslice <- 1
    i <- list(seq_along(allseq))
    seq <- write_sequence(
      allseq,
      fname = withr::local_tempfile(fileext = ".fasta")
    )
  }
  distmx <- withr::local_tempfile()
  system2("mkfifo", distmx)
  for (j in seq_len(nslice)) {
    system2(
      find_executable("dist_matrix"),
      args = c(
        "-l", aln_len,
        "-i", length(i[[j]]),
        "-m", 100,
        max(thresh),
        seq[j]
      ),
      stdout = distmx,
      wait = FALSE
    )
    if (j < nslice) {
      for (k in (2*j+1):(2*nslice)) {
        system2(
          find_executable("dist_bipart"),
          args = c(
            "-l", aln_len,
            "-i", length(i[[j]]),
            "-m", 100,
            max(thresh),
            seq_half[k],
            seq[j]
          ),
          stdout = distmx,
          wait = FALSE
        )
      }
    }
  }
  optimotu::distmx_cluster(
    distmx = distmx,
    names = as.character(which),
    threshold_config = optimotu::threshold_set(thresh),
    clust_config = optimotu::clust_tree(),
    parallel_config = optimotu::parallel_concurrent(max(1, local_cpus() - nslice))
  )

}

parse_protaxAnimal_output <- function(x) {
  out <-
    tibble::tibble(output = x) |>
    tidyr::separate(
      output,
      into = c("seq_id", "assignment"),
      extra = "merge",
      fill = "right",
      convert = TRUE
    ) |>
    dplyr::mutate(
      assignment = gsub("([^ ]+) ([0-9.]+)", "\\1\x1f\\2", assignment)
    ) |>
    tidyr::separate_longer_delim(assignment, " ") |>
    tidyr::separate(
      assignment,
      into = c("taxonomy", "prob"),
      sep = "\x1f",
      convert = TRUE
    ) |>
    dplyr::mutate(
      rank = int2rankfactor(
        stringr::str_count(taxonomy, ",") + 1L + length(KNOWN_RANKS)
      )
    ) |>
    tidyr::extract(taxonomy, into = c("parent_taxonomy", "taxon"), regex = "(?:(.+),)?([^,]+)$") |>
    dplyr::mutate(
      parent_taxonomy = paste(
        paste(KNOWN_TAXA, collapse = ","),
        parent_taxonomy,
        sep = ","
      ) |>
        trimws(whitespace = ","),
      taxon = dplyr::na_if(taxon, "unk")
    ) |>
    dplyr::select(seq_id, rank, parent_taxonomy, taxon, prob) |>
    dplyr::arrange(seq_id, rank)

  if (is.integer(out$seq_id)) out <- dplyr::rename(out, seq_idx = seq_id)
  out
}

fastq_names <- function(fq) {
  if (!file.exists(fq)) return(character())
  if (endsWith(fq, ".gz")) {
    system(paste("zcat", fq, "| awk 'NR%4==1{print substr($1, 2)}'"), intern = TRUE)
  } else {
    system(paste(" awk 'NR%4==1{print substr($1, 2)}'", fq), intern = TRUE)
  }
}


#' Quickly retrieve sequences from a gzipped file using an index
#'
#' @param file (`character` filename) file to create an index for
#'
#' @return file name of the created index
#' @rdname fastx_gz
fastx_gz_index <- function(file) {
  index <- sprintf("%s.fqi", file)
  args <- c(
    "index",
    sprintf("-f=%s", file),
    sprintf("-i=%s", index),
    "-w"
  )
  out = system2("bin/fastqindex_0.9.0b", args)
  stopifnot(out == 0)
  checkmate::assert_file_exists(index, "r")
  index
}

#' @param infile (`character` filename) gzipped fasta or fastq file
#' @param index (`character` filename) index file for `infile`
#' @param i (`integer` vector) indices to extract
#' @param outfile (`character` filename) file to write the extracted sequences to
#' @param renumber (`logical` flag) if `TRUE`, replace the sequence names with
#'   integers, starting at 0.
#' @param append (`logical` flag) if `TRUE`, append to `outfile` if it already
#'   exists, rather than overwriting.
#'
#' @return filename of the output file
#' @rdname fastx_gz
fastx_gz_extract <- function(infile, index, i, outfile, renumber = FALSE, append = FALSE, hash = NULL) {
  checkmate::assert_file_exists(infile, "r")
  checkmate::assert_file_exists(index, "r")
  checkmate::assert_integerish(i, lower = 1)
  checkmate::assert_string(outfile)
  checkmate::assert_flag(renumber)
  checkmate::assert_flag(append)
  if (file.exists(outfile) && !append) unlink(outfile)
  ensure_directory(outfile)
  if (!file.exists(outfile)) file.create(outfile)
  start <- which(i != dplyr::lag(i, 1, -1) + 1L)
  end <- c(start[-1] - 1L, length(i))
  is_fastq <- endsWith(infile, "fastq.gz") || endsWith(infile, "fq.gz")
  command <- sprintf(
    "bin/fastqindex_0.9.0b extract -s=%i -n=%i -e=%i -f=%s -i=%s | tail -n+8",
    i[start] - 1L,
    end - start + 1L,
    if (is_fastq) 4 else 2,
    infile,
    index
  )
  if (renumber) {
    command = paste(
      command,
      sprintf(
        "| awk -v n=%i 'NR%%%i==1{print \"%c\" n; n++; next}; {print}'",
        cumsum(dplyr::lag(end - start + 1L, 1L, 0L)),
        if (is_fastq) 4 else 2,
        if (is_fastq) "@" else ">"
      )
    )
  }
  if (endsWith(outfile, ".gz")) {
    command = paste(command, "| gzip -c -")
  }
  command = paste(command, ">>", outfile)
  result <- vapply(command, system, 0L)
  stopifnot(all(result == 0))
  outfile
}

#' @param infile (`character` filename) gzipped fasta or fastq file
#' @param index (`character` filename) index file for `infile`
#' @param i (`integer` vector) indices to extract
#' @param outfile (`character` filename) file to write the extracted sequences to
#' @param renumber (`logical` flag) if `TRUE`, replace the sequence names with
#'   integers, starting at 0.
#' @param append (`logical` flag) if `TRUE`, append to `outfile` if it already
#'   exists, rather than overwriting.
#'
#' @return filename of the output file
#' @rdname fastx_gz
fastx_gz_random_access_extract <- function(infile, index, i, outfile = NULL, renumber = FALSE, append = FALSE, hash = NULL, max_gap = 100L, ncpu = local_cpus()) {
  checkmate::assert_file_exists(infile, "r")
  checkmate::assert_file_exists(index, "r")
  checkmate::assert_integerish(i, lower = 1)
  checkmate::assert_string(outfile, null.ok = TRUE)
  checkmate::assert_flag(renumber)
  checkmate::assert_flag(append)
  checkmate::assert_integerish(max_gap, lower = 1)
  isort <- sort(unique(i))
  start <- which(isort > dplyr::lag(isort, 1, -max_gap) + max_gap)
  end <- c(start[-1] - 1L, length(i))
  is_fastq <- endsWith(infile, "fastq.gz") || endsWith(infile, "fq.gz")
  tmpfile <- replicate(length(start), withr::local_tempfile(fileext = ".fasta"))
  processes <- vector("list", length(start))
  for (j in seq_len(min(ncpu, length(start)))) {
    processes[[j]] <- processx::process$new(
      command = "bin/fastqindex_0.9.0b",
      arg = c(
        "extract",
        sprintf("-s=%d", isort[start[j]] - 1L),
        sprintf("-n=%d", isort[end[j]] - isort[start[j]] + 1L),
        sprintf("-e=%d", if (is_fastq) 4 else 2),
        sprintf("-f=%s", infile),
        sprintf("-i=%s", index),
        sprintf("-o=%s", tmpfile[j])
      ),
      stderr = ""
    )
  }
  j <- 1
  fastqindex_return <- 0
  while (j <= length(start) && !is.null(processes[[j]])) {
    fastqindex_return <- fastqindex_return + processes[[j]]$wait()$get_exit_status()
    if (fastqindex_return == 0 && (k <- j + ncpu) <= length(start)) {
      processes[[k]] <- processx::process$new(
        command = "bin/fastqindex_0.9.0b",
        arg = c(
          "extract",
          sprintf("-s=%d", isort[start[k]] - 1L),
          sprintf("-n=%d", isort[end[k]] - isort[start[k]] + 1L),
          sprintf("-e=%d", if (is_fastq) 4 else 2),
          sprintf("-f=%s", infile),
          sprintf("-i=%s", index),
          sprintf("-o=%s", tmpfile[k])
        ),
        stderr = ""
      )
    }
    j <- j + 1
  }
  stopifnot(fastqindex_return == 0)

  included <- unlist(mapply("seq", isort[start], isort[end]))
  selected <- match(i, included)
  seqs <- Biostrings::readBStringSet(tmpfile, seek.first.rec = TRUE)[selected]
  if (renumber) {
    names(seqs) <- as.character(seq_along(seqs))
  }
  if (is.null(outfile)) {
    seqs
  } else {
    write_sequence(seqs, outfile, compress = endsWith(outfile, ".gz"))
  }
}

fastx_gz_multi_extract <- function(infile, index, ilist, outfiles, renumber = FALSE, append = FALSE) {
  checkmate::assert_file_exists(infile, "r")
  checkmate::assert_file_exists(index, "r")
  checkmate::assert_list(i, types = "integerish")
  checkmate::assert_path_for_output(outfiles)
  stopifnot(length(ilist) == length(outfiles))
  checkmate::assert_flag(renumber)
  checkmate::assert_flag(append)
  for (of in outfiles) {
    if (file.exists(of) && !append) unlink(of)
    ensure_directory(of)
    if (!file.exists(of)) file.create(of)
  }
  start <- which(i != dplyr::lag(i, 1, -1) + 1L)
  end <- c(start[-1] - 1L, length(i))
  is_fastq <- endsWith(infile, "fastq.gz") || endsWith(infile, "fq.gz")
  command <- sprintf(
    "bin/fastqindex_0.9.0b extract -s=%i -n=%i -e=%i -f=%s -i=%s | tail -n+8",
    i[start] - 1L,
    end - start + 1L,
    if (is_fastq) 4 else 2,
    infile,
    index
  )
  if (renumber) {
    command = paste(
      command,
      sprintf(
        "| awk -v n=%i 'NR%%%i==1{print \"%c\" n; n++; next}; {print}'",
        cumsum(dplyr::lag(end - start + 1L, 1L, 0L)),
        if (is_fastq) 4 else 2,
        if (is_fastq) "@" else ">"
      )
    )
  }
  stopifnot(all(endsWith(outfiles, ".gz")) || all(!endsWith(outfiles, ".gz")))
  if (all(endsWith(outfiles, ".gz"))) {
    command = paste(command, "| gzip -c -")
  }
  command = paste(command, ">>", outfiles)
  result <- vapply(command, system, 0L)
  stopifnot(all(result == 0))
  outfiles
}


#' @param start (`integer` scalar) one-based index to start hashing
#' @param n (`integer` scalar) number of sequences to hash
#'
#' @return (`character`) md5 hash
fastx_gz_hash <- function(infile, index, start, n) {
  checkmate::assert_file_exists(infile, "r")
  checkmate::assert_file_exists(index, "r")
  checkmate::assert_integerish(start, lower = 1)
  checkmate::assert_integerish(n, lower = 1)
  is_fastq <- endsWith(infile, "fastq.gz") || endsWith(infile, "fq.gz")
  command <- sprintf(
    "bin/fastqindex_0.9.0b extract -s=%i -n=%i -e=%i -f=%s -i=%s | tail -n+8 | md5sum",
    start - 1L,
    n,
    if (is_fastq) 4 else 2,
    infile,
    index
  )
  result <- system(command, intern = TRUE)
  stopifnot(attr(result, "status") == 0)
  c(strtrim(result, 32))
}

fastx_rename <- function(infile, names, outfile) {
  checkmate::assert_file_exists(infile, "r")
  checkmate::assert_character(names)
  if (length(names) != 1 || !file.exists(names)) {
    stopifnot(sequence_size(infile) == length(names))
    names <- write_and_return_file(
      names,
      withr::local_tempfile(fileext = ".txt")
    )
  }
  if (grepl(fastq_regex, infile)) {
    header_condition <- "NR%4==1"
    header_token <- "@"
  } else if (grepl(fasta_regex, infile)) {
    header_condition <- "/^>/"
    header_token <- ">"
  } else {
    stop("Cannot determine file type for ", infile)
  }
  command <- sprintf(
    "awk '%s{getline name < \"%s\"; print \"%s\" name; next}; {print}'",
    header_condition,
    names,
    header_token
  )
  if (endsWith(infile, ".gz")) {
    command <- paste("zcat", infile, "|", command)
  } else {
    command <- paste(command, "<", infile)
  }
  if (endsWith(outfile, ".gz")) {
    command <- paste(command, "| gzip -c - >", outfile)
  } else {
    command <- paste(command, ">", outfile)
  }
  result <- system(command)
  stopifnot(result == 0L)
  outfile
}

# splits a (possibly gzipped) fastx file into n subfiles
# never loads anything into memory, if subfiles are compressed it happens in
# parallel
# this would be more portable with a custom c function (I think?)
fastx_split <- function(infile, n, outroot = tempfile(), compress = FALSE) {
  checkmate::assert_string(infile)
  checkmate::assert_file(infile, access = "r")
  checkmate::assert_int(n, lower = 1, upper = 64)
  checkmate::assert_path_for_output(outroot)
  checkmate::assert_flag(compress)

  is_fastq <- grepl(fastq_regex, infile)
  is_gz <- endsWith(infile, ".gz")

  if (n == 1 && is_gz == compress) return(infile)

  suffix <- if(is_fastq) ".fastq" else ".fasta"
  if (compress) suffix <- paste0(suffix, ".gz")
  command <- if (is_gz) "zcat" else "cat"
  command <- paste(command, infile, "| paste -d'\u1f' - -")
  if (is_fastq) command <- paste(command, "- -")

  digits <- floor(log10(n)) + 1
  command <- paste(
    command,
    " | split",
    sprintf("-nr/%i", n),
    "--numeric-suffixes=1",
    "-a", digits,
    "--additional-suffix", suffix,
    "--filter='tr \"\u1f\" \"\\n\""
    )
  if (compress) command <- paste(command, "| gzip -c -")
  command <- paste(
    command, ">$FILE'",
    "-", outroot
  )
  command <- paste("bash -c", shQuote(command))
  result <- system(command)
  stopifnot(result == 0)
  outfiles <- sprintf(paste0("%s%0", digits, "i%s"), outroot, seq_len(n), suffix)
  checkmate::assert_file(outfiles, "r")
  outfiles
}

# rejoins some split fastx files
# of course this is most useful if something happened to them in between.
# If no lines are missing or reordered in the infiles, then the outfile
# will end up in the same order as the original file which was split.
# Otherwise this will almost certainly not happen.
fastx_combine <- function(infiles, outfile) {
  checkmate::assert_file(infiles, "r")
  checkmate::assert_path_for_output(outfile, overwrite = TRUE)
  is_fastq <- grepl(fastq_regex, infiles)
  stopifnot(all(is_fastq) | all(!is_fastq))
  is_fastq <- all(is_fastq)

  is_gz <-endsWith(infiles, ".gz")
  stopifnot(all(is_gz) | all(!is_gz))
  is_gz <- all(is_gz)

  compress <- endsWith(outfile, ".gz")

  command <- 'paste -d"\\n"'
  subcommand <- if (is_gz) "<( zcat" else "<( cat"
  subcommand <- paste(subcommand, infiles, "| paste -d'\u1f' - -")
  if (is_fastq) subcommand <- paste(subcommand, "- -")
  subcommand <- paste(subcommand, ")", collapse = " ")
  command <- paste(command, subcommand, "| tr '\u1f' \"\\n\" | sed -n /^$/!p")
  if (compress) command <- paste(command, "| gzip -c -")
  command <- paste(command, ">", outfile)
  command <- paste("bash -c", shQuote(command))
  result <- system(command)
  stopifnot(result == 0)
  checkmate::assert_file(outfile, "r")
  outfile
}
