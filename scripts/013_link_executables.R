# put links to programs in the place where sh_matching will look for them
# (./bin will be bound to /sh_matching/programs)

link_executable <- function(filename) {
  unlink(filename)
  ensure_directory(filename)
  file.symlink(Sys.which(basename(filename)), filename)
  filename
}

withr::with_dir(
  "bin",
  {
    link_executable("krona/bin/ktImportText")
    link_executable("vsearch/bin/vsearch")
  }
)
