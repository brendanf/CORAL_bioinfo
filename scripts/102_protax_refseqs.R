# Add additional reference sequences to the Protax reference data (if provided)

if (protax_aligned) {
  default_model_dir <- protax_root
  taxonomy_filename <- "taxonomy.priors"
} else {
  default_model_dir <- file.path(protax_root, "addedmodel")
  taxonomy_filename <- "taxonomy"
}


protax_usearch <- file.path(protax_root, "scripts", "usearch10.0.240_i86linux32")

refseq_plan <- list(
  #### taxonomy_default_file ####
  # character: path and file name (tsv file, no extension)
  #
  # The default Protax taxonomy file
  tar_file_fast(
    taxonomy_default_file,
    file.path(default_model_dir, taxonomy_filename)
  ),

  #### taxonomy_default ####
  # tibble:
  #  `taxon_id` integer : taxon index
  #  `parent_id` integer : index of parent taxon
  #  `rank` integer : 0 = root, 1=ROOT_RANK..n=TIP_RANK
  #  `classification` character : full comma-separated classification of this
  #    taxon, not including "root"
  #  `prior` numeric : prior for a new sequence to belong to this taxon; by
  #    default equal to the number of species belonging to this taxon divided
  #    by the total number of species
  tar_fst_tbl(
    taxonomy_default,
    readr::read_tsv(
      taxonomy_default_file,
      col_names = c("taxon_id", "parent_id", "rank", "classification", "prior"),
      col_types = "iiicd",
      col_select = 1:5
    )
  )
)

if (checkmate::test_file_exists(pipeline_options$added_reference$fasta) &&
    checkmate::test_file_exists(pipeline_options$added_reference$table)) {
  custom_protax_dir <- "custom_protax"
  if (!dir.exists(custom_protax_dir)) dir.create(custom_protax_dir)

  refseq_plan <- list(
    refseq_plan,

    #### taxonomy_ascii7_default_file ####
    # character: path and file name (tsv format, no extension)
    #
    # The default Protax taxonomy file, modified to use only ascii characters
    tar_file_fast(
      taxonomy_ascii7_default_file,
      file.path(default_model_dir, "taxonomy.ascii7"),
      deployment = "main"
    ),

    #### taxonomy_ascii7_default ####
    # tibble:
    #  `taxon_id` integer : taxon index
    #  `parent_id` integer : index of parent taxon
    #  `rank` integer : 0 = root, 1=ROOT_RANK..n=TIP_RANK
    #  `classification` character : full comma-separated classification of this
    #    taxon, not including "root"; modified to include only ascii characters
    #  `prior` numeric : prior for a new sequence to belong to this taxon; by
    #    default equal to the number of species belonging to this taxon divided
    #    by the total number of species
    tar_fst_tbl(
      taxonomy_ascii7_default,
      readr::read_tsv(
        taxonomy_ascii7_default_file,
        col_names = c("taxon_id", "parent_id", "rank", "classification", "prior"),
        col_types = "iiicd"
      ),
      deployment = "main"
    ),

    #### new_refseq_file ####
    # character : path and file name (fasta format)
    #
    # user-provided reference sequences to be added to Protax
    tar_file_fast(
      new_refseq_file,
      pipeline_options$added_reference$fasta,
      deployment = "main"
    ),

    #### new_refseq ####
    # Biostrings::DNAStringSet : user-provided reference sequences to be added
    #   to Protax
    tar_target(
      new_refseq,
      Biostrings::readDNAStringSet(new_refseq_file),
      deployment = "main"
    ),

    #### new_refseq_metadata_file ####
    # character : path and file name (.xlsx)
    tar_file_fast(
      new_refseq_metadata_file,
      pipeline_options$added_reference$table,
      deployment = "main"
    ),

    #### new_refseq_metadata ####
    # tibble:
    #  `Culture_ID` character : ID which should match the name of the sequence
    #    in `new_refseq`
    #  `Protax_synonym` character : full comma-separated classification of this
    #    taxon, not including "root"
    #  ... : there may be other columns in the provided file, but they are not
    #    used.
    tar_fst_tbl(
      new_refseq_metadata,
      # TODO: accept csv/tsv/ods here too
      readxl::read_excel(new_refseq_metadata_file) %>%
        dplyr::filter(Culture_ID %in% names(new_refseq)),
      deployment = "main"
    ),

    #### taxonomy_new ####
    # tibble:
    #  `taxon_id` integer : taxon index
    #  `parent_id` integer : index of parent taxon
    #  `rank` integer : 0 = root, 1=ROOT_RANK..n=TIP_RANK
    #  `classification` character : full comma-separated classification of this
    #    taxon, not including "root"
    #  `prior` numeric : prior for a new sequence to belong to this taxon; by
    #    default equal to the number of species belonging to this taxon divided
    #    by the total number of species
    #
    # Add user-provided taxa to Protax's taxonomy
    tar_fst_tbl(
      taxonomy_new,
      build_taxonomy_new(
        taxonomy_default$classification,
        new_refseq_metadata$Protax_synonym
      ),
      deployment = "main"
    ),

    #### write_protax_taxonomy_new ####
    # character : path and file name (tsv format, no extension)
    #
    # write taxonomy_new to the new Protax model directory
    tar_file_fast(
      write_protax_taxonomy_new,
      write_and_return_file(
        taxonomy_new,
        file.path(custom_protax_dir, "taxonomy"),
        "tsv",
        col_names = FALSE
      ),
      deployment = "main"
    ),

    #### write_its2_new ####
    # character : path and file name (.fa, fasta format)
    #
    # generate the new Protax sequence reference by appending the user-provided
    # sequences to the existing references
    tar_file_fast(
      write_its2_new,
      {
        outfile <- file.path(custom_protax_dir, "its2.fa")
        file.copy(file.path(default_model_dir, "its2.fa"), outfile, overwrite = TRUE)
        file.append(outfile, new_refseq_file)
        outfile
      },
      deployment = "main"
    ),

    #### write_sintaxits2_new ####
    # character: path and file name (.fa, fasta format)
    #
    # generate the new reference sequence file for Sintax by appending the user-
    # provided file to the Protax default file, after formatting the FASTA
    # headers to include the taxonomy in Sintax format
    tar_file_fast(
      write_sintaxits2_new,
      {
        outfile <- file.path(custom_protax_dir, "sintaxits2train.fa")
        file.copy(file_path(default_model_dir, "sintaxits2train.fa"), outfile, overwrite = TRUE)
        tibble::enframe(as.character(new_refseq), name = "Culture_ID") %>%
          dplyr::left_join(
            dplyr::select(new_refseq_metadata, Culture_ID, Protax_synonym),
            by = "Culture_ID"
          ) %>%
          dplyr::transmute(
            name = paste(Culture_ID, sintax_format(Protax_synonym), sep = ";"),
            value = value
          ) %>%
          tibble::deframe() %>%
          write_sequence(outfile, append = TRUE)
        outfile
      },
      deployment = "main"
    ),
    #### write_its2udb_new ####
    # character : path and file name (*.udb, usearch database format)
    #
    # convert the new reference sequences to UDB for fast searching
    tar_file_fast(
      write_its2udb_new,
      build_udb(
        write_its2_new,
        file.path(custom_protax_dir, "its2.udb"),
        type = "usearch",
        usearch = protax_usearch
      )
    ),
    #### write_sintaxits2udb_new ####
    # character : path and file name (*.udb, usearch database format)
    #
    # convert the new reference sequences to UDB for fast taxonomic assignment
    # by Sintax
    tar_file_fast(
      write_sintaxits2udb_new,
      build_udb(
        write_sintaxits2_new,
        file.path(custom_protax_dir, "sintaxits2.udb"),
        type = "sintax",
        usearch = protax_usearch
      )
    ),
    #### write amptksynmockudb ####
    # character : path and file name (*.udb, usearch databse format)
    #
    # nothing needs to be changed about the old synmock file, so just copy it
    # over.
    tar_file_fast(
      write_amptksynmockudb,
      {
        outfile <- file.path(custom_protax_dir, "amptk_synmock.udb")
        file.symlink(
          fs::path_rel(
            path = file.path(default_model_dir, "amptk_synmock.udb"),
            start = custom_protax_dir),
          outfile
        )
        outfile
      },
      deployment = "main"
    ),

    #### write_protax_taxonomy.ascii7_new ####
    # character : path and file name (no extension, TSV format)
    #
    # write the ascii-cleaned version of the new taxonomy file
    tar_file_fast(
      write_protax_taxonomy.ascii7_new,
      write_and_return_file(
        dplyr::mutate(
          taxonomy_new,
          classification = ascii_clean(classification)
        ),
        file.path(custom_protax_dir, "taxonomy.ascii7"),
        "tsv",
        col_names = FALSE
      ),
      deployment = "main"
    ),
    #### repeat for ranks from 2 (phylum) to 7 (species) ####
    tar_map(
      values = list(.rank = 2:7),

      ##### write_protax_tax_{.rank} #####
      # character : path and file name (no extension, TSV format)
      #
      # taxonomy file, truncated to the chosen rank
      tar_file_fast(
        write_protax_tax,
        write_and_return_file(
          dplyr::filter(taxonomy_new, rank <= .rank),
          file.path(custom_protax_dir, paste0("tax", .rank)),
          "tsv",
          col_names = FALSE
        ),
        deployment = "main"
      ),

      ##### write_protax_ref.tax_{.rank} #####
      # character : path and file name (no extension, TSV format)
      #
      # append classification of user-provided sequences to Protax files;
      # this is "for each reference sequence, what is its classification"
      tar_file_fast(
        write_protax_ref.tax,
        {
          outfile <- paste0("ref.tax", .rank)
          file.copy(
            from = file.path(default_model_dir, outfile),
            to = file.path(custom_protax_dir, outfile),
            overwrite = TRUE
          )
          dplyr::transmute(
            new_refseq_metadata,
            Culture_ID = Culture_ID,
            Protax_synonym = truncate_taxonomy(Protax_synonym, .rank)
          ) |>
            dplyr::filter(!is.na(Protax_synonym)) |>
            write_and_return_file(
              file.path(custom_protax_dir, outfile),
              "tsv",
              append = TRUE
            )
        },
        deployment = "main"
      ),

      ##### write_protax_rseqs #####
      # character : path and file name (no extension, TSV format)
      #
      # add classification of user-provided sequences to Protax files;
      # this is "for each classification, what are its reference sequences"
      tar_file_fast(
        write_protax_rseqs,
        dplyr::bind_rows(
          readr::read_tsv(
            file.path(default_model_dir, sprintf("rseqs%d", .rank)),
            col_names = c("taxon_id", "accno"),
            col_types = "ic"
          ) %>%
            tidyr::separate_rows(accno, sep = ","),
          dplyr::transmute(
            new_refseq_metadata,
            accno = Culture_ID,
            classification = truncate_taxonomy(Protax_synonym, .rank)
          ) %>%
            dplyr::left_join(
              dplyr::filter(taxonomy_new, rank == .rank),
              by = "classification"
            ) %>%
            dplyr::select(taxon_id, accno)
        ) %>%
          dplyr::group_by(taxon_id) %>%
          dplyr::summarize(accno = paste(accno, collapse = ",")) %>%
          write_and_return_file(
            file.path(custom_protax_dir, sprintf("rseqs%d", .rank)),
            type = "tsv",
            col_names = FALSE
          ),
        deployment = "main"
      )
    )
  )

  # programatically find all of the output files from the refseqplan so far,
  # and make a target to require all of them
  custom_protax_files <- purrr::keep(get_target_names(refseq_plan), startsWith, "write_")

  paste("{", paste(custom_protax_files, collapse = "\n"), shQuote(custom_protax_dir), "}", sep = "\n")

  refseq_plan <- c(
    refseq_plan,
    #### custom_protax ####
    # character : directory name
    #
    # requires all of the modified files, and returns the directory name
    tar_target_raw(
      name = "custom_protax",
      command = parse(
        text = paste(
          "{",
          paste(custom_protax_files, collapse = "\n"),
          shQuote(custom_protax_dir),
          "}",
          sep = "\n"
        )
      ),
      format = "file",
      deployment = "main"
    ),
    #### protax_model ####
    # character : directory name
    #
    # The directory containing the protax references to use; in this case use
    # the version with user-supplied references
    tar_file_fast(
      protax_model,
      custom_protax,
      deployment = "main"
    )
  )
} else {
  refseq_plan <- list(
    refseq_plan,
    tar_fst_tbl(
      taxonomy_new,
      taxonomy_default,
      deployment = "main"
    ),

    #### protax_model ####
    # character : directory name
    #
    # The directory containing the protax references to use; in this case use
    # the default version
    tar_file_fast(
      protax_model,
      default_model_dir,
      deployment = "main"
    )
  )
}

optimotu_plan <- c(optimotu_plan, refseq_plan)
