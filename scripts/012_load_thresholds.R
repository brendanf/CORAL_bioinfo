# test the performance of different clustering thresholds for replicating
# "true" taxonomic groups.

fmeasure_optima <- readr::read_tsv(threshold_file, col_types = "ccccdd")
