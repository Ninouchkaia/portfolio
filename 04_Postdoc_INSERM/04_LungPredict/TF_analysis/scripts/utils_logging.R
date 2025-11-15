# scripts/utils_logging.R

init_logger <- function(log_path = "results/logs/pipeline.log") {
  dir.create(dirname(log_path), showWarnings = FALSE, recursive = TRUE)
  assign(".lungpredict_logfile", log_path, envir = .GlobalEnv)
}

log_msg <- function(..., level = "INFO") {
  logfile <- get(".lungpredict_logfile", envir = .GlobalEnv, inherits = FALSE)
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  msg <- paste0("[", ts, "][", level, "] ", paste(..., collapse = " "))
  cat(msg, "\n", file = logfile, append = TRUE)
  invisible(msg)
}

log_info  <- function(...) log_msg(..., level = "INFO")
log_debug <- function(...) log_msg(..., level = "DEBUG")
