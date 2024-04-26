# define functions and metadata for the plan
for (f in list.files("scripts", "^0[[:digit:]]{2}_.+[.]R$", full.names = TRUE)) {
  source(f)
}

targets::tar_option_set(
  # by default, workers run the targets, retrieve their prerequisites, and
  # store their results
  deployment = "worker",
  storage = "worker",
  retrieval = "worker",
  error = "abridge",
  resources = targets::tar_resources(
    clustermq = targets::tar_resources_clustermq(
      template = list(
        job_name = "worker_its2_taxonomy_first",
        cores = 8,
        time = "12:00:00",
        temp_space=10,
        memory=16384
      )
    )
  )
)
options(
  clustermq.scheduler = "slurm",
  clustermq.template = file.path(getwd(), "slurm", "puhti_clustermq.tmpl")
)

cat("Running pipeline with a pool of", n_workers, "clustermq workers.\n")

target <- strsplit(Sys.getenv("OPTIMOTU_TARGET"), "[, ;]")[[1]]
if (length(target) > 0) {
  targets::tar_make_clustermq(
    names = any_of(target),
    callr_function=NULL,
    workers = n_workers,
    reporter = "timestamp"
  )
} else {
  targets::tar_make_clustermq(
    callr_function=NULL,
    workers = n_workers,
    reporter = "timestamp"
  )
}
