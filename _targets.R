# Bioinformatics pipeline for Sonja Saine's spruce log metabarcoding project
#   by Brendan Furneaux
# Based on:
#   DADA2 analysis for GSSP from Jenni Hultman
#   PROTAX taxonomy assignment for GSSP by Panu Somervuo
#   Clustering for GSSP by Brendan Furneaux

library(targets)
library(tarchetypes)
library(magrittr)
library(qs)
library(fst)

tar_option_set(
  format = "qs",
  memory = "transient",
  garbage_collection = TRUE,
  priority = 0.5,
  workspace_on_error = TRUE
)

optimotu_plan <- list()

# Numbered R scripts define the targets plan.
# They are numbered in the order they are used.

# certain runners (all except run_node) run the 0** scripts, because they need
# information from the configuration file.  Only run them now if they have not
# been run before.
if (!exists("pipeline_options")) {
  for (f in list.files("scripts", "^0[[:digit:]]{2}_.+[.]R$", full.names = TRUE)) {
    source(f)
  }
}

# the 1** and higher scritps are the ones which actually define the plan.
for (f in list.files("scripts", "^[1-9][[:digit:]]{2}_.+[.]R$", full.names = TRUE)) {
  source(f)
}

cat("Detected", local_cpus(), "cores for main process.\n" )

# Make sure log directory exists
if (!dir.exists("logs")) dir.create("logs")

# End this file with a list of target objects.
invisible(optimotu_plan)
