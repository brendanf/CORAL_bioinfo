output_targets <- rlang::syms(
  tarchetypes::tar_select_names(optimotu_plan, starts_with("write_"))
  )

optimotu_plan <- c(
  optimotu_plan,
  list(
    tar_file_fast(
      zip_output,
      {
        outfile <- sprintf(
          "%s/%s_%s.zip",
          "output",
          project_name,
          strftime(Sys.Date(), "%Y%m%d")
        )
        i <- 1
        while(file.exists(outfile)) {
          outfile <- sprintf(
            "%s/%s_%s_%i.zip",
            "output",
            project_name,
            strftime(Sys.Date(), "%Y%m%d"),
            i
          )
          i <- i + 1
        }
        result <- zip(outfile, c(!!!output_targets), zip = "zip", flags = "-j9X")
        stopifnot(result == 0)
        outfile
      }
    )
  )
)
