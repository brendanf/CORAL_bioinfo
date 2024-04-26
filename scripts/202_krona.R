krona_plan <- list(

  #### krona_script ####
  # character: KronaTools script for embedding in Krona plots
  tar_target(
    krona_script,
    readLines(
      withr::local_connection(
        url("http://marbl.github.io/Krona/src/krona-2.0.js")
      )
    ),
    deployment = "main"
  ),

  #### krona_shortcut_icon ####
  # character: base64-encoded KronaTools shortcut icon for embedding in Krona
  #   plots
  tar_target(
    krona_shortcut_icon,
    base64enc::dataURI(
      withr::local_connection(
        url("http://marbl.github.io/Krona//img/favicon.ico", open = "rb")
      ),
      mime = "image/x-icon"
    ),
    deployment = "main"
  ),

  #### krona_hiddenimage ####
  # character: base64-encoded KronaTools "hidden" icon for embedding in Krona
  #   plots
  tar_target(
    krona_hiddenimage,
    base64enc::dataURI(
      withr::local_connection(
        url("http://marbl.github.io/Krona//img/hidden.png", open = "rb")
      ),
      mime = "image/png"
    ),
    deployment = "main"
  ),

  #### krona_loadingimage ####
  # character: base64-encoded KronaTools "loading" icon for embedding in Krona
  #   plots
  tar_target(
    krona_loadingimage,
    base64enc::dataURI(
      withr::local_connection(
        url("http://marbl.github.io/Krona//img/loading.gif", open = "rb")
      ),
      mime = "image/gif"
    ),
    deployment = "main"
  ),

  #### krona_logo ####
  # character: base64-encoded KronaTools logo for embedding in Krona plots
  tar_target(
    krona_logo,
    base64enc::dataURI(
      withr::local_connection(
        url("http://marbl.github.io/Krona//img/logo-small.png", open = "rb")
      ),
      mime = "image/png"
    ),
    deployment = "main"
  ),

  tar_map(
    values = tibble::tibble(
      .conf_level = c("plausible", "reliable"),
      otu_taxonomy = paste0("otu_taxonomy_", .conf_level) %>%
        rlang::syms()
    ),
    names = .conf_level,
    #### otu_krona_data_{.conf_level} ####
    # tibble:
    #  `rank` ordered factor: species:kingdom, rank of the taxon for this line
    #  `taxon` character: name of taxon for this line
    #  `parent_taxonomy` character: comma-delimited parent classification
    #  `phylum_unknown_fread` numeric: fraction of reads in this taxon with
    #    unknown phylum
    #  `phylum_unknown_fotu` numeric: fraction of OTUs in this taxon with
    #    unknown phylum
    #  `phylum_unknown_focc` numeric: fraction of occurrences in this taxon with
    #    unknown phylum
    #  `class_unknown_fread` numeric: fraction of reads in this taxon with
    #    unknown class
    #  `class_unknown_fotu` numeric: fraction of OTUs in this taxon with
    #    unknown class
    #  `class_unknown_focc` numeric: fraction of occurrences in this taxon with
    #    unknown class
    #  `order_unknown_fread` numeric: fraction of reads in this taxon with
    #    unknown order
    #  `order_unknown_fotu` numeric: fraction of OTUs in this taxon with
    #    unknown order
    #  `order_unknown_focc` numeric: fraction of occurrences in this taxon with
    #    unknown order
    #  `family_unknown_fread` numeric: fraction of reads in this taxon with
    #    unknown family
    #  `family_unknown_fotu` numeric: fraction of OTUs in this taxon with
    #    unknown family
    #  `family_unknown_focc` numeric: fraction of occurrences in this taxon with
    #    unknown family
    #  `genus_unknown_fread` numeric: fraction of reads in this taxon with
    #    unknown genus
    #  `genus_unknown_fotu` numeric: fraction of OTUs in this taxon with
    #    unknown genus
    #  `genus_unknown_focc` numeric: fraction of occurrences in this taxon with
    #    unknown genus
    #  `species_unknown_fread` numeric: fraction of reads in this taxon with
    #    unknown species
    #  `species_unknown_fotu` numeric: fraction of OTUs in this taxon with
    #    unknown species
    #  `species_unknown_focc` numeric: fraction of occurrences in this taxon with
    #    unknown species
    #  `nread` integer: total number of reads for this taxon
    #  `nocc` integer: total number of occurrences for this taxon
    #  `notu` integer: total number of otus in this taxon
    #  `child_unknown_fread` numeric: fraction of reads in this taxon which are
    #    unidentified at the child rank
    #  `child_unknown_fotu` numeric: fraction of OTUs in this taxon which are
    #    unidentified at the child rank
    #  `child_unknown_focc` numeric: fraction of occurrences in this taxon which are
    #    unidentified at the child rank
    #  `fread` numeric: fraction of all reads belonging to this taxon
    #  `focc` numeric: raction of all occurrences belonging to this taxon
    #  `fotu` numeric: raction of all otus belonging to this taxon
    #
    tar_fst_tbl(
      otu_krona_data,
      if (nrow(otu_taxonomy) == 0) {
        tibble::tibble(
          rank = rank2factor(character()),
          taxon = character(),
          parent_taxonomy = character(),
          phylum_unknown_fread = numeric(),
          phylum_unknown_fotu = numeric(),
          phylum_unknown_focc = numeric(),
          class_unknown_fread = numeric(),
          class_unknown_fotu = numeric(),
          class_unknown_focc = numeric(),
          order_unknown_fread = numeric(),
          order_unknown_fotu = numeric(),
          order_unknown_focc = numeric(),
          family_unknown_fread = numeric(),
          family_unknown_fotu = numeric(),
          family_unknown_focc = numeric(),
          genus_unknown_fread = numeric(),
          genus_unknown_fotu = numeric(),
          genus_unknown_focc = numeric(),
          species_unknown_fread = numeric(),
          species_unknown_fotu = numeric(),
          species_unknown_focc = numeric(),
          nread = integer(),
          nocc = integer(),
          notu = integer(),
          child_unknown_fread = numeric(),
          child_unknown_fotu = numeric(),
          child_unknown_focc = numeric(),
          fread = numeric(),
          focc = numeric(),
          fotu = numeric()
        )
      } else {
        otu_taxonomy |>
          dplyr::mutate(
            genus = remove_mycobank_number(genus),
            species = remove_mycobank_number(species),
            phylum_parent = kingdom,
            class_parent = paste(phylum_parent, phylum, sep = ","),
            order_parent = paste(class_parent, class, sep = ","),
            family_parent = paste(order_parent, order, sep = ","),
            genus_parent = paste(family_parent, family, sep = ","),
            species_parent = paste(genus_parent, genus, sep = ","),
            dplyr::across(
              .cols = phylum:species,
              .fns = \(x) startsWith(x, "pseudo"),
              .names = "{.col}_unknown"
            )
          ) |>
          dplyr::rename_with(.fn = paste0, .cols = kingdom:species, "_taxon") |>
          tidyr::pivot_longer(
            kingdom_taxon:species_parent,
            names_to = c("rank", ".value"),
            names_sep = "_",
            names_transform = list(rank = rank2factor)
          ) |>
          dplyr::mutate(taxon = chartr("_", " ", taxon)) |>
          dplyr::group_by(rank, taxon, parent) |>
          dplyr::summarize(
            dplyr::across(
              phylum_unknown:species_unknown,
              list(
                fread = ~sum(nread * .)/sum(nread),
                fotu = ~sum(.)/dplyr::n(),
                focc = ~sum(nsample*.)/(sum(nsample))
              ),
              .names = "{.col}_{.fn}"
            ),
            nread = sum(nread),
            nocc = sum(nsample),
            notu = dplyr::n(),
            .groups = "drop"
          ) |>
          dplyr::mutate(
            child_unknown_fread = dplyr::case_when(
              rank == "kingdom" ~ phylum_unknown_fread,
              rank == "phylum" ~ class_unknown_fread,
              rank == "class" ~ order_unknown_fread,
              rank == "order" ~ family_unknown_fread,
              rank == "family" ~ genus_unknown_fread,
              TRUE ~ species_unknown_fread
            ),
            child_unknown_focc = dplyr::case_when(
              rank == "kingdom" ~ phylum_unknown_focc,
              rank == "phylum" ~ class_unknown_focc,
              rank == "class" ~ order_unknown_focc,
              rank == "order" ~ family_unknown_focc,
              rank == "family" ~ genus_unknown_focc,
              TRUE ~ species_unknown_focc
            ),
            child_unknown_fotu = dplyr::case_when(
              rank == "kingdom" ~ phylum_unknown_fotu,
              rank == "phylum" ~ class_unknown_fotu,
              rank == "class" ~ order_unknown_fotu,
              rank == "order" ~ family_unknown_fotu,
              rank == "family" ~ genus_unknown_fotu,
              TRUE ~ species_unknown_fotu
            )
          ) |>
          dplyr::group_by(rank) |>
          dplyr::mutate(
            fread = nread/sum(nread),
            focc = nocc/sum(nocc),
            fotu = notu/sum(notu)
          ) |>
          dplyr::rename(parent_taxonomy = parent)
      },
      deployment = "main"
    ),

    #### write_otu_krona_{.conf_level} ####
    # character (output filename)
    #
    # write a stand-alone HTML file containing the Krona plot
    tar_file_fast(
      write_otu_krona,
      sprintf("output/otu_krona_%s.html", .conf_level) |>
        krona_xml_nodes(
          data = dplyr::filter(otu_krona_data, (nocc>=5)|(notu>=5)|(nread>1000)),
          .rank = ROOT_RANK,
          maxrank = TIP_RANK,
          outfile = _,
          node_data_format = list(
            f = c("focc", "fread", "fotu"),
            nocc = rep("nocc", 3),
            nread = rep("nread", 3),
            notu = rep("notu", 3),
            sp = c("species_unknown_focc", "species_unknown_fread",
                   "species_unknown_fotu"),
            gen = c("genus_unknown_focc", "genus_unknown_fread",
                    "genus_unknown_fotu"),
            fam = c("family_unknown_focc", "family_unknown_fread",
                    "family_unknown_fotu")
          ),
          taxonomy = NULL,
          pre = c(
            '<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">',
            '<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">',
            ' <head>',
            '  <meta charset="utf-8"/>',
            paste0('  <link rel="shortcut icon" href="', krona_shortcut_icon, '"/>'),
            '  <script id="notfound" type="text/javascript">window.onload=function(){document.body.innerHTML=""}</script>',
            '  <script language="javascript" type="text/javascript">',
            krona_script,
            '  </script>',
            ' </head>',
            ' <body>',
            paste0('  <img id="hiddenImage" src="', krona_hiddenimage, '" style="display:none" alt="Hidden Image"/>'),
            paste0('  <img id="loadingImage" src="', krona_loadingimage, '" style="display:none" alt="Loading Indicator"/>'),
            paste0('  <img id="logo" src="', krona_logo, '" style="display:none" alt="Logo of Krona"/>'),
            '  <noscript>Javascript must be enabled to view this page.</noscript>',
            '  <div style="display:none">',
            '<krona>',
            '<attributes magnitude="f">',
            '<attribute display="Weighted fraction">f</attribute>',
            '<attribute display="Total occurences">nocc</attribute>',
            '<attribute display="Total reads">nread</attribute>',
            '<attribute display="Total OTUs">notu</attribute>',
            '<attribute display="Weighted fraction belonging to unknown species">sp</attribute>',
            '<attribute display="Weighted fraction belonging to unknown genera">gen</attribute>',
            '<attribute display="Weighted fraction belonging to unknown families">fam</attribute>',
            '</attributes>',
            '<color attribute="sp" valueStart="0" valueEnd="1" hueStart="120" hueEnd="0" default="true"></color>',
            '<datasets>',
            '<dataset>Occurence weighting</dataset>',
            '<dataset>Read abundance weighting</dataset>',
            '<dataset>OTU richness weighting</dataset>',
            '</datasets>'
          ),
          post = "</krona>"
        ),
      deployment = "main"
    )
  )
)

optimotu_plan <- c(optimotu_plan, krona_plan)
