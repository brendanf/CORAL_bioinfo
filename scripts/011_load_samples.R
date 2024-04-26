library(magrittr)

#### find all the sequencing files and load the sequencing metadata

#### Get the samples organized ####

#define paths
path <- "."
meta_path <- file.path(path, "metadata")
seq_path <- file.path(path, "sequences")
raw_path <- file.path(seq_path, "01_raw")
trim_path <- file.path(seq_path, "02_trim")
filt_path <- file.path(seq_path, "03_filter")
asv_path <- file.path(seq_path, "04_denoised")
protax_path <- file.path(seq_path, "05_protax")
output_path <- file.path(path, "output")
log_path <- file.path(path, "logs")

# create paths if missing
if (!dir.exists(filt_path)) dir.create(filt_path, recursive = TRUE)
if (!dir.exists(trim_path)) dir.create(trim_path, recursive = TRUE)
if (!dir.exists(asv_path)) dir.create(asv_path, recursive = TRUE)
if (!dir.exists(protax_path)) dir.create(protax_path, recursive = TRUE)
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)
if (!dir.exists(log_path)) dir.create(log_path, recursive = TRUE)

true_vals <- c("1", "y", "Y", "yes", "Yes", "YES", "t", "T", "true", "True", "TRUE")
false_vals <- c("0", "n", "N", "no", "No", "NO", "f", "F", "false", "False", "FALSE")

if (!isFALSE(pipeline_options$custom_sample_table)) {
  #todo: support sample table in other formats (csv, excel, ...)
  sample_table <- suppressWarnings(
    readr::read_tsv(
      pipeline_options$custom_sample_table,
      col_types = readr::cols(
        seqrun = readr::col_character(),
        sample = readr::col_character(),
        neg_control = readr::col_character(),
        pos_control = readr::col_character(),
        fastq_R1 = readr::col_character(),
        fastq_R2 = readr::col_character(),
        orient = readr::col_character(),
        .default = readr::col_guess()
      )
    )
  )
  checkmate::check_data_frame(
    sample_table,
    col.names = "named"
  )
  checkmate::check_names(
    names(sample_table),
    must.include = c("sample", "seqrun", "fastq_R1", "fastq_R2")
  )
  checkmate::assert_character(sample_table$sample, any.missing = FALSE)
  checkmate::assert_character(sample_table$seqrun, any.missing = FALSE)
  checkmate::assert_file_exists(file.path(raw_path, sample_table$fastq_R1), access = "r")
  checkmate::assert_file_exists(file.path(raw_path, sample_table$fastq_R2), access = "r")
  if ("pos_control" %in% names(sample_table)) {
    checkmate::assert_subset(sample_table$pos_control, c(true_vals, false_vals))
    sample_table$pos_control <- sample_table$pos_control %in% true_vals
  } else {
    sample_table$pos_control <- FALSE
  }
  if ("neg_control" %in% names(sample_table)) {
    checkmate::assert_subset(sample_table$neg_control, c(true_vals, false_vals))
    sample_table$neg_control <- sample_table$neg_control %in% true_vals
  } else {
    sample_table$neg_control <- FALSE
  }
  if ("orient" %in% names(sample_table)) {
    checkmate::assert_subset(sample_table$orient, c("fwd", "rev", "mixed"))
    if (pipeline_options$orient != "custom") {
      warning(
        "custom sample table '", pipeline_options$custom_sample_table,
        "' includes an 'orient' column, but option 'orient' is '",
        pipeline_options$orient, "'. The 'orient' column will be ignored."
      )
      sample_table$orient <- NULL
    } else {
      if (any(sample_table$orient == "mixed")) {
        sample_table <-
          dplyr::left_join(
            sample_table,
            tibble::tibble(
              orient = c("fwd", "rev", "mixed", "mixed"),
              new_orient = c("fwd", "rev", "fwd", "rev")
            ),
            by = "orient",
            multiple = "all"
          ) |>
          dplyr::mutate(orient = new_orient, .keep = "unused")
      }
    }
  }
  sample_table <- dplyr::mutate(
    sample_table,
    dplyr::across(
      any_of(c("truncQ_R1", "truncQ_R2", "cut_R1", "cut_R2")),
      \(x) lapply(strsplit(as.character(x), ","), as.numeric)
    )
  )
} else {
  # find files
  sample_table <- tibble::tibble(
    fastq_R1 = sort(list.files(raw_path, paste0(".*R1(_001)?[.]", pipeline_options$file_extension), recursive = TRUE)),
    fastq_R2 = sort(list.files(raw_path, paste0(".*R2(_001)?[.]", pipeline_options$file_extension), recursive = TRUE))
  ) %>%
    # parse filenames
    tidyr::extract(
      fastq_R1,
      into = c("seqrun", "sample"),
      regex = paste0("([^/]+)/(?:.*/)?(.+?)[._](?:S\\d+_L001_)?R1(?:_001)?[.]", pipeline_options$file_extension),
      remove = FALSE
    ) %>%
    dplyr::mutate(
      sample = dplyr::if_else(
        startsWith(sample, "BLANK"),
        paste(seqrun, sample, sep = "_"),
      sample
    )
  )
}

switch(
  pipeline_options$orient,
  fwd = sample_table$orient <- "fwd",
  rev = sample_table$orient <- "rev",
  mixed = sample_table <- tidyr::crossing(sample_table, orient = c("fwd", "rev")),
  custom = if (isFALSE(pipeline_options$custom_sample_table)) {
    stop("option 'orient: custom' requires a custom sample table is given.")
  } else if (!"orient" %in% names(sample_table)) {
    stop("option 'orient: custom' required a column named 'orient' in the",
    " custom sample table, with values consisting of 'fwd', 'rev', and 'mixed'")
  },
  stop("unknown value for option 'orient'; should be 'fwd', 'rev', 'mixed', or 'custom'")
)

sample_table <- sample_table %>%
  # generate filenames for trimmed and filtered reads
  dplyr::mutate(
    sample_key = paste(seqrun, sample, sep = "_"),
    trim_R1 = file.path(
      trim_path,
      paste(sample_key, orient, "R1_trim.fastq.gz", sep = "_")
    ),
    trim_R2 = file.path(
      trim_path,
      paste(sample_key, orient, "R2_trim.fastq.gz", sep = "_")
    ),
    filt_R1 = file.path(
      filt_path,
      paste(sample_key, orient, "R1_filt.fastq.gz", sep = "_")
    ),
    filt_R2 = file.path(
      filt_path,
      paste(sample_key, orient, "R2_filt.fastq.gz", sep = "_")
    ),
    sample_key = file_to_sample_key(filt_R1) # to be sure
  )

# spike_strength is used along with the nonspike/spike ratio to convert from
# read number to "weight"
if (!("spike_weight") %in% names(sample_table))
  sample_table$spike_weight <- 1

assertthat::assert_that(
  !any(is.na(sample_table$seqrun)),
  !any(is.na(sample_table$sample)),
  is.numeric(sample_table$spike_weight),
  !any(duplicated(sample_table[c("fastq_R1", "orient")])),
  !any(duplicated(sample_table[c("fastq_R2", "orient")])),
  !any(duplicated(sample_table$trim_R1)),
  !any(duplicated(sample_table$trim_R2)),
  !any(duplicated(sample_table$filt_R1)),
  !any(duplicated(sample_table$filt_R2))
)

n_seqrun <- dplyr::n_distinct(sample_table$seqrun)
n_orient_seqrun <- dplyr::n_distinct(sample_table$seqrun, sample_table$orient)
n_workers <- max(min_workers, min(max_workers, n_orient_seqrun*workers_per_seqrun))

sample_table_key <- dplyr::select(
  sample_table,
  sample,
  seqrun,
  sample_key,
) |>
  unique()

cat("Found", dplyr::n_distinct(sample_table$sample, sample_table$seqrun),
    "samples in", n_seqrun, "runs.\n",
    "sample_table targets hash is:", targets:::digest_obj64(sample_table), "\n"
)
for (n in colnames(sample_table)) {
  cat(sprintf("sample_table$%s hash: %s\n", n, targets:::digest_obj64(sample_table[[n]])))
}
