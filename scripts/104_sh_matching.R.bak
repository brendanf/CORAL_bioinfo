# parse Unite SH assignments
# copied from GSSP 
library(magrittr)
library(targets)
library(tarchetypes)

unite_thresholds <- tibble::tibble(
  sim_percent = seq(97, 99.5, 0.5),
  dist_percent = 100-sim_percent,
  sim_frac = sim_percent/100,
  dist_frac = dist_percent/100,
  sim_name = sub(".", "", as.character(sim_percent), fixed = TRUE),
  dist_name = sub(".", "", as.character(10*dist_frac), fixed = TRUE)
)

sh_datafile <- "sh_matching_data_0_5.zip"
sh_dataurl <- "https://files.plutof.ut.ee/public/orig/9C/FD/9CFD7C58956E5331F1497853359E874DEB639B17B04DB264C8828D04FA964A8F.zip"

# do we need to do ITSx?
# typically this means that the clustering thresholds are calculated for ITS2,
# but our sequences also include flanking regions.
do_itsx <- Sys.getenv("DO_ITSX", names = FALSE) != ""
# if this pipeline does ITSx, then it doesn't need to also be done in sh_matching.
do_itsx_in_sh_matching <- if (isTRUE(do_itsx)) "no" else "yes"

SH_plan <- list(
  #### sh_meta ####
  # split the data up into the same number of groups as we have sequencing
  # runs, because that is the number of workers we have most likely chosen.
  tar_fst_tbl(
    sh_meta,
    tibble::tibble(
      jobnumber = seq.int(n_seqrun),
      sh_infile = sprintf("indata/source_%d", jobnumber),
      sh_outfile = sprintf("outdata/source_%d.zip", jobnumber)
    )
  ),

  if (do_itsx) {
    list(
      #### asv_its2_seq ####
      # extract only the ITS2 region
      tar_fst_tbl(
        asv_its2_seq,
        dplyr::select(grouped_asv_seq, ASV, seq) %>%
          tibble::deframe() %>%
          Biostrings::DNAStringSet() %>%
          rITSx::itsx(
            taxon = "all",
            complement = FALSE,
            cpu = local_cpus(),
            summary = FALSE,
            graphical = FALSE,
            fasta = FALSE,
            preserve = TRUE,
            save_regions = "ITS2",
            positions = FALSE,
            not_found = FALSE,
            read_function = Biostrings::readDNAStringSet
          ) %>%
          magrittr::extract2("ITS2") %>%
          tibble::enframe(value = "seq"),
        pattern = map(grouped_asv_seq)
      ),
      #### asvs_to_unite ####
      # write a fasta file of non-spike ASV sequences for SH matching
      tar_file(
        asvs_to_unite,
        write_sequence(asv_its2_seq, sh_meta$sh_infile),
        priority = 1,
        deployment = "main",
        pattern = map(asv_its2_seq, sh_meta)
      )
    )
  } else {
    #### asvs_to_unite ####
    # write a fasta file of non-spike ASV sequences for SH matching
    tar_file(
      asvs_to_unite,
      write_sequence(grouped_asv_seq, sh_meta$sh_infile),
      priority = 1,
      deployment = "main",
      pattern = map(grouped_asv_seq, sh_meta)
    )
  },
  #### sh_matching_script ####
  tar_file(
    sh_matching_script,
    "sh_matching_pub/sh_matching_analysis/run_pipeline.sh"
  ),
  #### sh_matching_analysis ####
  tar_file(
    sh_matching_analysis,
    list.files(path="sh_matching_pub/sh_matching_analysis/scripts", full.names = TRUE)
  ),
  #### sh_matching_data ####
  tar_file(
    sh_matching_data,
    {
      download.file(
        url = sh_dataurl,
        destfile = sh_datafile,
        quiet = TRUE
      )
      out <- unzip(
        zipfile = sh_datafile,
        exdir = "data/sh_matching_data",
        junkpaths = TRUE,
        overwrite = TRUE
      )
      unlink(sh_datafile)
      out
    }
  ),
  #### sh_matching ####
  # run the SH matching pipeline locally
  tar_file(
    sh_matching,
    {
      # mention these here just for dependency tracking
      asvs_to_unite
      sh_matching_analysis
      sh_matching_data
      # remove the output file if it exists
      unlink(sh_meta$sh_outfile, force = TRUE)
      ensure_directory(sh_meta$sh_outfile)
      # make sure the temp directory is there
      # running on csc this should be mapped to the local scratch drive
      if (!dir.exists("userdir")) dir.create("userdir")
      # run the pipeline
      system2(sh_matching_script, c(sh_meta$jobnumber, "its2", do_itsx_in_sh_matching))
      # return the output file (for dependency tracking)
      sh_meta$sh_outfile
    },
    priority = 1,
    pattern = map(asvs_to_unite, sh_meta)
  ),
  
  #### unite_excluded ####
  # look up the sequences which Unite excluded, and why.
  # typically these will be based on chimeras or length
  tar_fst_tbl(
    unite_excluded,
    dplyr::left_join(
      readr::read_tsv(
        unz(sh_matching, sprintf("excluded_%d.txt", sh_meta$jobnumber)),
        col_types = "ccc",
        col_names = c("seq_id_tmp", "step", "reason")
      ),
      readr::read_tsv(
        unz(sh_matching, sprintf("source_%s_names", sh_meta$jobnumber)),
        col_types = "cci",
        col_names = c("seq_accno", "seq_id_tmp", "n")
      ),
      by = "seq_id_tmp"
    ),
    pattern = map(sh_matching, sh_meta)
  ),
  tar_map(
    values = unite_thresholds,
    names = sim_name,
    #### unite_matches_out ####
    # parse the Unite matches_out* files
    # matches_out is for sequences with at least 75% similarity to a sequence
    # already in Unite.
    # matches_out_1 is for sequences which do NOT have a match at 75% similarity
    # to a sequence already in Unite.
    tar_fst_tbl(
      unite_matches_out,
      purrr::map_dfr(
        c("matches", "matches_1"),
        # for some reason readr::read_tsv had parsing errors
        ~ readLines(
          unz(sh_matching, sprintf("matches/%s_out_%s.csv", .x, dist_name))
        ) %>%
          tibble::tibble(data = .) %>%
          tidyr::separate(
            col = data,
            into = unlist(strsplit(.$data[1], "\t")),
            sep = "\t",
            extra = "drop"
          ) %>%
          dplyr::slice(-1)
      ) %>%
        # the Unite column names include the threshold (as dissimilarity),
        # but we are taking care of that by keeping the data separate.
        dplyr::rename_with(.fn = sub, pattern = " \\(\\d\\.\\d\\)",
                           replacement = "") %>%
        tidyr::separate(
          `SH/compound taxonomy`,
          into = c("kingdom", "phylum", "class", "order", "family", "genus", "species"),
          sep = ";",
          fill = "right",
          remove = FALSE
        ) %>%
        dplyr::mutate(
          dplyr::across(kingdom:species, sub, pattern = "^[kpcofgs]__", replacement = ""),
          dplyr::across(kingdom:species, dplyr::na_if, "unidentified"),
          species = ifelse(endsWith(species, "_sp"), NA_character_, species)
        ),
      pattern = map(sh_matching)
    )
  )
)
