

#### define alternate version of crew_controller_slurm ####
# this is a workaround for crew worker failures when too many try to dial in to
# the master at the same time.

# code modified from crew.cluster
if (requireNamespace("crew.cluster", quietly = TRUE)) {

  crew_launcher_slurm2 <- function(
    name = NULL,
    seconds_interval = 0.5,
    seconds_timeout = 60,
    seconds_launch = 86400,
    seconds_idle = Inf,
    seconds_wall = Inf,
    tasks_max = Inf,
    tasks_timers = 0L,
    reset_globals = TRUE,
    reset_packages = FALSE,
    reset_options = FALSE,
    garbage_collection = FALSE,
    launch_max = 5L,
    tls = crew::crew_tls(mode = "automatic"),
    verbose = FALSE,
    command_submit = as.character(Sys.which("sbatch")),
    command_terminate = as.character(Sys.which("scancel")),
    command_delete = NULL,
    script_directory = tempdir(),
    script_lines = character(0L),
    slurm_log_output = "/dev/null",
    slurm_log_error = "/dev/null",
    slurm_memory_gigabytes_per_cpu = NULL,
    slurm_cpus_per_task = NULL,
    slurm_time_minutes = 1440,
    slurm_partition = NULL
  ) {
    name <- as.character(name %||% crew::crew_random_name())
    if (!is.null(command_delete)) {
      crew::crew_deprecate(
        name = "command_delete",
        date = "2023-01-08",
        version = "0.1.4.9001",
        alternative = "command_terminate"
      )
      command_terminate <- command_delete
    }
    launcher <- crew_class_launcher_slurm2$new(
      name = name,
      seconds_interval = seconds_interval,
      seconds_timeout = seconds_timeout,
      seconds_launch = seconds_launch,
      seconds_idle = seconds_idle,
      seconds_wall = seconds_wall,
      tasks_max = tasks_max,
      tasks_timers = tasks_timers,
      reset_globals = reset_globals,
      reset_packages = reset_packages,
      reset_options = reset_options,
      garbage_collection = garbage_collection,
      launch_max = launch_max,
      tls = tls,
      verbose = verbose,
      command_submit = command_submit,
      command_terminate = command_terminate,
      script_directory = script_directory,
      script_lines = script_lines,
      slurm_log_output = slurm_log_output,
      slurm_log_error = slurm_log_error,
      slurm_memory_gigabytes_per_cpu = slurm_memory_gigabytes_per_cpu,
      slurm_cpus_per_task = slurm_cpus_per_task,
      slurm_time_minutes = slurm_time_minutes,
      slurm_partition = slurm_partition
    )
    launcher$validate()
    launcher
  }

  crew_class_launcher_slurm2 <- R6::R6Class(
    classname = "crew_class_launcher_slurm2",
    inherit = crew.cluster::crew_class_launcher_slurm,
    private = list(
      .prefix = NULL,
      #### added ####
      .last_job = NULL,
      .args_launch = function(script) {
        shQuote(
          c(
            "--parsable",
            if (!is.null(private$.last_job))
              paste0("--dependency=after:", private$.last_job)
            else
              NULL,
            script
          )
        )
      }
      #### end added ####
    ),
    public = list(
      launch_worker = function(call, name, launcher, worker, instance) {
        lines <- c(self$script(name = name), paste("R -e", shQuote(call)))
        if (is.null(private$.prefix)) {
          if (!file.exists(super$script_directory)) {
            dir.create(super$script_directory, recursive = TRUE)
          }
          private$.prefix <- crew::crew_random_name()
        }
        script <- crew.cluster:::path_script(
          dir = super$script_directory,
          prefix = private$.prefix,
          launcher = launcher,
          worker = worker
        )
        writeLines(text = lines, con = script)
        #### modified ####
        launch_result <- system2(
          command = super$command_submit,
          args = private$.args_launch(script = script),
          stdout = TRUE,
          stderr = crew.cluster:::if_any(super$verbose, "", FALSE),
          wait = TRUE
        )
        if (is.null(attr(launch_result, "status"))) {
          # `sbatch --parsable` returns JobID, optionally followed by semicolon and cluster name.
          # add more sanity checks here
          job <- strsplit(launch_result, ";")[[1]][1]
          private$.last_job <- job
        }
        if (any(super$verbose)) {
          ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%OS2", tz = "GMT")
          targets:::cli_blue_bullet(paste(ts, "submitted job", launch_result))
        }
        #### end modified ####
        list(name = name, script = script)
      }
    )
  )

  crew_controller_slurm2 <- function(
    name = NULL,
    workers = 1L,
    host = NULL,
    port = NULL,
    tls = crew::crew_tls(mode = "automatic"),
    tls_enable = NULL,
    tls_config = NULL,
    seconds_interval = 0.25,
    seconds_timeout = 60,
    seconds_launch = 86400,
    seconds_idle = Inf,
    seconds_wall = Inf,
    seconds_exit = NULL,
    tasks_max = Inf,
    tasks_timers = 0L,
    reset_globals = TRUE,
    reset_packages = FALSE,
    reset_options = FALSE,
    garbage_collection = FALSE,
    launch_max = 5L,
    verbose = FALSE,
    command_submit = as.character(Sys.which("sbatch")),
    command_terminate = as.character(Sys.which("scancel")),
    command_delete = NULL,
    script_directory = tempdir(),
    script_lines = character(0L),
    slurm_log_output = "/dev/null",
    slurm_log_error = "/dev/null",
    slurm_memory_gigabytes_per_cpu = NULL,
    slurm_cpus_per_task = NULL,
    slurm_time_minutes = 1440,
    slurm_partition = NULL
  ) {
    if (!is.null(seconds_exit)) {
      crew::crew_deprecate(
        name = "seconds_exit",
        date = "2023-09-21",
        version = "0.5.0.9002",
        alternative = "none (no longer necessary)"
      )
    }
    client <- crew::crew_client(
      name = name,
      workers = workers,
      host = host,
      port = port,
      tls = tls,
      tls_enable = tls_enable,
      tls_config = tls_config,
      seconds_interval = seconds_interval,
      seconds_timeout = seconds_timeout
    )
    launcher <- crew_launcher_slurm2(
      name = name,
      seconds_interval = seconds_interval,
      seconds_timeout = seconds_timeout,
      seconds_launch = seconds_launch,
      seconds_idle = seconds_idle,
      seconds_wall = seconds_wall,
      tasks_max = tasks_max,
      tasks_timers = tasks_timers,
      reset_globals = reset_globals,
      reset_packages = reset_packages,
      reset_options = reset_options,
      garbage_collection = garbage_collection,
      launch_max = launch_max,
      tls = tls,
      verbose = verbose,
      command_submit = command_submit,
      command_terminate = command_terminate,
      command_delete = command_delete,
      script_directory = script_directory,
      script_lines = script_lines,
      slurm_log_output = slurm_log_output,
      slurm_log_error = slurm_log_error,
      slurm_memory_gigabytes_per_cpu = slurm_memory_gigabytes_per_cpu,
      slurm_cpus_per_task = slurm_cpus_per_task,
      slurm_time_minutes = slurm_time_minutes,
      slurm_partition = slurm_partition
    )
    controller <- crew::crew_controller(client = client, launcher = launcher)
    controller$validate()
    controller
  }
}
