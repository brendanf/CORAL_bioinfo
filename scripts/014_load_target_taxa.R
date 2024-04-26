target_taxa <- if (file.exists("data/target_taxa.txt")) {
  readLines("data/target_taxa.txt")
} else {
  character()
}
