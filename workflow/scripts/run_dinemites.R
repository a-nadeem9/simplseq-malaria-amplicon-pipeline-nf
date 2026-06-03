#!/usr/bin/env Rscript
# run_dinemites.R
# ---------------------------------------------------------------------------
# Runs the DINEMITES package on prepared longitudinal input data.
# Supports three model types: simple, clustering, bayesian.
# Produces allele probabilities, molFOI, new infections, and per-subject plots.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

# ── CLI arguments ──────────────────────────────────────────────────────────
option_list <- list(
  make_option("--input", type = "character",
              help = "Path to DINEMITES input TSV (from simplseq_to_dinemites.R)"),
  make_option("--model", type = "character", default = "simple",
              help = "Model type: simple, clustering, or bayesian [default %default]"),
  make_option("--outdir", type = "character",
              help = "Output directory for DINEMITES results"),
  make_option("--n_lags", type = "integer", default = 3,
              help = "Simple model: number of recent samples to check [default %default]"),
  make_option("--t_lag", type = "character", default = "Inf",
              help = "Simple model: time window in days, or Inf for no day cutoff [default %default]"),
  make_option("--seed", type = "integer", default = 1,
              help = "Random seed for stochastic DINEMITES models [default %default]"),
  make_option("--refresh", type = "integer", default = 100,
              help = "Stan/clustering progress refresh interval [default %default]"),
  make_option("--bayesian_lag_days", type = "integer", default = 30,
              help = "Lag window in days for generated Bayesian covariates [default %default]"),
  make_option("--bayesian_chains", type = "integer", default = 1,
              help = "Bayesian model Stan chains [default %default]"),
  make_option("--bayesian_parallel_chains", type = "integer", default = 1,
              help = "Bayesian model parallel Stan chains [default %default]"),
  make_option("--bayesian_iter_warmup", type = "integer", default = 500,
              help = "Bayesian model warmup iterations [default %default]"),
  make_option("--bayesian_iter_sampling", type = "integer", default = 500,
              help = "Bayesian model sampling iterations [default %default]"),
  make_option("--bayesian_adapt_delta", type = "double", default = 0.99,
              help = "Bayesian model Stan adapt_delta [default %default]"),
  make_option("--bayesian_drop_out", type = "character", default = "false",
              help = "Use the Bayesian drop-out model: true/false [default %default]"),
  make_option("--plot_width", type = "double", default = 18,
              help = "DINEMITES plot width in inches [default %default]"),
  make_option("--plot_height", type = "double", default = 0,
              help = "DINEMITES plot height in inches; 0 auto-scales [default %default]")
)
args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args$input) || is.null(args$outdir)) {
  stop("Required arguments: --input, --outdir")
}

model_type <- tolower(args$model)
if (!model_type %in% c("simple", "clustering", "bayesian")) {
  stop("Model type must be one of: simple, clustering, bayesian. Got: ", model_type)
}

parse_t_lag <- function(value) {
  raw_value <- trimws(as.character(value))
  if (tolower(raw_value) %in% c("", "inf", "infinity", "none")) {
    return(Inf)
  }
  parsed <- suppressWarnings(as.numeric(raw_value))
  if (is.na(parsed) || parsed < 0) {
    stop("--t_lag must be a non-negative number of days or Inf. Got: ", raw_value)
  }
  parsed
}

parse_bool <- function(value) {
  raw_value <- tolower(trimws(as.character(value)))
  if (raw_value %in% c("1", "true", "yes", "y", "on")) {
    return(TRUE)
  }
  if (raw_value %in% c("0", "false", "no", "n", "off", "")) {
    return(FALSE)
  }
  stop("Boolean value must be true or false. Got: ", value)
}

if (is.na(args$n_lags) || args$n_lags < 1) {
  stop("--n_lags must be at least 1. Got: ", args$n_lags)
}
args$t_lag <- parse_t_lag(args$t_lag)
args$bayesian_drop_out <- parse_bool(args$bayesian_drop_out)

if (is.na(args$seed) || args$seed < 1) {
  stop("--seed must be at least 1. Got: ", args$seed)
}
if (is.na(args$refresh) || args$refresh < 0) {
  stop("--refresh must be at least 0. Got: ", args$refresh)
}
if (is.na(args$bayesian_lag_days) || args$bayesian_lag_days < 1) {
  stop("--bayesian_lag_days must be at least 1. Got: ", args$bayesian_lag_days)
}
if (is.na(args$bayesian_chains) || args$bayesian_chains < 1) {
  stop("--bayesian_chains must be at least 1. Got: ", args$bayesian_chains)
}
if (is.na(args$bayesian_parallel_chains) || args$bayesian_parallel_chains < 1) {
  stop("--bayesian_parallel_chains must be at least 1. Got: ",
       args$bayesian_parallel_chains)
}
if (is.na(args$bayesian_iter_warmup) || args$bayesian_iter_warmup < 1) {
  stop("--bayesian_iter_warmup must be at least 1. Got: ",
       args$bayesian_iter_warmup)
}
if (is.na(args$bayesian_iter_sampling) || args$bayesian_iter_sampling < 1) {
  stop("--bayesian_iter_sampling must be at least 1. Got: ",
       args$bayesian_iter_sampling)
}
if (is.na(args$bayesian_adapt_delta) ||
    args$bayesian_adapt_delta <= 0 ||
    args$bayesian_adapt_delta >= 1) {
  stop("--bayesian_adapt_delta must be greater than 0 and less than 1. Got: ",
       args$bayesian_adapt_delta)
}

cat("[DINEMITES/run] Input:", args$input, "\n")
cat("[DINEMITES/run] Model:", model_type, "\n")
cat("[DINEMITES/run] Output dir:", args$outdir, "\n")
cat("[DINEMITES/run] Simple windows: n_lags=", args$n_lags,
    ", t_lag=", ifelse(is.infinite(args$t_lag), "Inf", args$t_lag), "\n", sep = "")
cat("[DINEMITES/run] Stochastic settings: seed=", args$seed,
    ", refresh=", args$refresh, "\n", sep = "")
cat("[DINEMITES/run] Bayesian settings: lag_days=", args$bayesian_lag_days,
    ", chains=", args$bayesian_chains,
    ", parallel_chains=", args$bayesian_parallel_chains,
    ", warmup=", args$bayesian_iter_warmup,
    ", sampling=", args$bayesian_iter_sampling,
    ", adapt_delta=", args$bayesian_adapt_delta,
    ", drop_out=", args$bayesian_drop_out, "\n", sep = "")

build_time_axis_labels <- function(dataset) {
  time_map <- dataset %>%
    dplyr::distinct(.data$time, .keep_all = TRUE) %>%
    dplyr::arrange(.data$time)

  labels <- as.character(time_map$time)

  if ("date_label" %in% colnames(time_map)) {
    date_labels <- trimws(as.character(time_map$date_label))
    labels <- ifelse(!is.na(date_labels) & nchar(date_labels) > 0,
                     paste0(time_map$time, "\n", date_labels),
                     labels)
  } else if ("date_full" %in% colnames(time_map)) {
    dates <- tryCatch(as.Date(time_map$date_full),
                      error = function(e) rep(as.Date(NA), nrow(time_map)))
    labels <- ifelse(!is.na(dates),
                     paste0(time_map$time, "\n", format(dates, "%b %Y")),
                     labels)
  }

  stats::setNames(labels, as.character(time_map$time))
}

build_allele_key <- function(dataset) {
  if (!"allele" %in% colnames(dataset)) {
    return(data.frame(short_allele_id = character(),
                      locus = character(),
                      allele = character(),
                      stringsAsFactors = FALSE))
  }

  locus_values <- if ("locus" %in% colnames(dataset)) {
    dataset$locus
  } else {
    rep("Allele", nrow(dataset))
  }

  key_source <- data.frame(
    locus = as.character(locus_values),
    allele = as.character(dataset$allele),
    stringsAsFactors = FALSE
  ) %>%
    mutate(locus = ifelse(is.na(.data$locus) | trimws(.data$locus) == "",
                          "Allele", .data$locus)) %>%
    filter(!is.na(.data$allele), trimws(.data$allele) != "") %>%
    distinct(.data$locus, .data$allele) %>%
    arrange(.data$locus, .data$allele) %>%
    group_by(.data$locus) %>%
    mutate(short_allele_id = paste0(.data$locus, "-", sprintf("%02d", row_number()))) %>%
    ungroup()

  key_source %>%
    select(short_allele_id, locus, allele)
}

build_allele_axis_labeler <- function(allele_key) {
  if (is.null(allele_key) || nrow(allele_key) == 0) {
    return(function(values) as.character(values))
  }
  label_map <- stats::setNames(allele_key$short_allele_id, allele_key$allele)
  function(values) {
    full_values <- as.character(values)
    short_values <- unname(label_map[full_values])
    ifelse(is.na(short_values) | short_values == "", full_values, short_values)
  }
}

safe_subject_dir_name <- function(value) {
  safe_value <- gsub("[^A-Za-z0-9._-]+", "_", as.character(value))
  safe_value <- gsub("^_+|_+$", "", safe_value)
  ifelse(nchar(safe_value) == 0, "subject", safe_value)
}

ensure_cmdstan_path <- function() {
  current_path <- tryCatch(cmdstanr::cmdstan_path(),
                           error = function(e) "")
  if (nzchar(current_path) && dir.exists(current_path)) {
    return(current_path)
  }

  rscript_path <- Sys.which("Rscript")
  rscript_prefix <- if (nzchar(rscript_path)) {
    normalizePath(file.path(dirname(rscript_path), ".."),
                  winslash = "/", mustWork = FALSE)
  } else {
    ""
  }
  r_home_prefix <- normalizePath(file.path(R.home(), "..", ".."),
                                 winslash = "/", mustWork = FALSE)

  candidates <- unique(c(
    Sys.getenv("CMDSTAN", unset = ""),
    file.path(Sys.getenv("CONDA_PREFIX", unset = ""), "bin", "cmdstan"),
    file.path(r_home_prefix, "bin", "cmdstan"),
    file.path(rscript_prefix, "bin", "cmdstan")
  ))
  candidates <- candidates[nzchar(candidates) & dir.exists(candidates)]
  if (length(candidates) == 0) {
    stop("[DINEMITES/run] ERROR: CmdStanR is installed, but no CmdStan directory was found.")
  }

  cmdstanr::set_cmdstan_path(candidates[1])
  cat("[DINEMITES/run] CmdStan path:", candidates[1], "\n")
  invisible(candidates[1])
}

compile_dinemites_stan_model <- function(model_type, bayesian_drop_out = FALSE) {
  package_path <- system.file(package = "dinemites")
  model_files <- instantiate::stan_package_model_files(package_path)
  model_pattern <- if (model_type == "clustering") {
    "model_infection_probabilities_clusters\\.stan$"
  } else if (bayesian_drop_out) {
    "model_infection_probabilities_bayesian_drop_out\\.stan$"
  } else {
    "model_infection_probabilities_bayesian\\.stan$"
  }
  selected_models <- model_files[grepl(model_pattern, model_files)]
  if (length(selected_models) != 1) {
    stop("[DINEMITES/run] ERROR: Could not find packaged Stan model for ",
         model_type, ".")
  }
  cat("[DINEMITES/run] Preparing Stan model:", basename(selected_models), "\n")
  instantiate::stan_package_compile(models = selected_models, quiet = TRUE)
  invisible(selected_models)
}

PREVALENCE_STRIP_LABEL <- "Allele prevalence (% visits)"
READABLE_PLOT_THEME <- theme(
  text = element_text(size = 10, color = "grey10"),
  plot.title = element_text(size = 14, color = "grey10"),
  axis.title = element_text(size = 10, color = "grey10"),
  axis.text = element_text(size = 9, color = "grey10"),
  strip.text = element_text(size = 10, color = "grey10"),
  legend.text = element_text(size = 9, color = "grey10"),
  legend.title = element_text(size = 10, color = "grey10")
)

prevalence_strip <- function() {
  display_label <- sub(" \\(", "\n(", PREVALENCE_STRIP_LABEL)
  ggplot(data.frame(l = display_label, x = 1, y = 1)) +
    geom_text(aes(.data$x, .data$y, label = .data$l), angle = 270, size = 4.2,
              lineheight = 0.95, color = "grey10") +
    theme_void() +
    coord_cartesian(clip = "off")
}

apply_time_axis_to_plot <- function(plot_out, x_breaks, x_labels, x_limit_max,
                                    allele_labeler = function(values) as.character(values)) {
  patchwork_plots <- list()
  if (inherits(plot_out, "patchwork")) {
    patchwork_plots <- plot_out$plots
    if (length(patchwork_plots) == 0 && !is.null(plot_out$patches$plots)) {
      patchwork_plots <- plot_out$patches$plots
    }
  }

  x_scale <- scale_x_continuous(
    breaks = x_breaks,
    labels = x_labels,
    limits = c(min(x_breaks), x_limit_max)
  )
  x_label <- labs(x = "Day\nMonth/Year")
  x_theme <- theme(axis.text.x = element_text(size = 9, color = "grey10", lineheight = 0.95))
  allele_axis_theme <- theme(
    axis.text.y = element_text(size = 8, color = "grey10"),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(size = 10, color = "grey10")
  )
  top_axis_theme <- theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )

  if (inherits(plot_out, "patchwork") && length(patchwork_plots) >= 3) {
    top_plot <- suppressMessages(patchwork_plots[[1]] + x_scale + labs(x = NULL) +
                                   READABLE_PLOT_THEME + top_axis_theme)
    allele_plot <- suppressMessages(
      patchwork_plots[[2]] + x_scale + x_label + labs(y = "Allele ID") +
        scale_y_discrete(labels = allele_labeler) + READABLE_PLOT_THEME +
        x_theme + allele_axis_theme
    )
    prevalence_label <- prevalence_strip()

    return(patchwork::wrap_plots(
      list(top_plot, allele_plot, prevalence_label),
      design = c(patchwork::area(1, 1),
                 patchwork::area(2, 1),
                 patchwork::area(2, 2, 2, 2)),
      heights = c(1, 5),
      widths = c(18, 1.5)
    ))
  }

  if (inherits(plot_out, "patchwork") && length(patchwork_plots) == 2) {
    top_plot <- suppressMessages(patchwork_plots[[1]] + x_scale + labs(x = NULL) +
                                   READABLE_PLOT_THEME + top_axis_theme)
    allele_plot <- suppressMessages(
      patchwork_plots[[2]] + x_scale + x_label + labs(y = "Allele ID") +
        scale_y_discrete(labels = allele_labeler) + READABLE_PLOT_THEME +
        x_theme + allele_axis_theme
    )
    return(patchwork::wrap_plots(
      list(top_plot, allele_plot, prevalence_strip()),
      design = c(patchwork::area(1, 1),
                 patchwork::area(2, 1),
                 patchwork::area(2, 2, 2, 2)),
      heights = c(1, 5),
      widths = c(18, 1.5)
    ))
  }

  suppressMessages(
    plot_out + x_scale + x_label + labs(y = "Allele ID") +
      scale_y_discrete(labels = allele_labeler) + READABLE_PLOT_THEME +
      x_theme + allele_axis_theme
  )
}

# ── Check DINEMITES package availability ───────────────────────────────────
if (!requireNamespace("dinemites", quietly = TRUE)) {
  stop("[DINEMITES/run] ERROR: dinemites package not installed. ",
       "Install with: devtools::install_github('WillNickols/dinemites')")
}

library(dinemites)

# Check Stan package availability for Bayesian/clustering models.
if (model_type %in% c("bayesian", "clustering")) {
  stan_package_names <- c("instantiate", "cmdstanr", "rstan", "posterior")
  if (model_type == "clustering") {
    stan_package_names <- c(stan_package_names, "linkcomm")
  }
  missing_stan_packages <- stan_package_names[
    !vapply(stan_package_names, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_stan_packages) > 0) {
    stop("[DINEMITES/run] ERROR: ", model_type,
         " model requires missing R package(s): ",
         paste(missing_stan_packages, collapse = ", "))
  }
  invisible(ensure_cmdstan_path())
  invisible(compile_dinemites_stan_model(model_type, args$bayesian_drop_out))
}

# ── 1. Read input data ────────────────────────────────────────────────────
dataset <- read.csv(args$input, header = TRUE, sep = "\t",
                    stringsAsFactors = FALSE)

cat("[DINEMITES/run] Loaded", nrow(dataset), "rows (",
    length(unique(dataset$subject)), "subjects,",
    length(unique(dataset$allele)), "alleles )\n")

# ── 2. Fill in dataset (complete allele × subject × time grid) ─────────────
cat("[DINEMITES/run] Filling in dataset (complete grid)...\n")
dataset_filled <- fill_in_dataset(dataset)
cat("[DINEMITES/run] Filled dataset:", nrow(dataset_filled), "rows\n")

# ── 3. Run selected model ─────────────────────────────────────────────────
cat("[DINEMITES/run] Running", model_type, "model...\n")
t_start <- Sys.time()

if (model_type == "simple") {
  results <- determine_probabilities_simple(dataset_filled,
                                             n_lags = args$n_lags,
                                             t_lag  = args$t_lag)
} else if (model_type == "clustering") {
  results <- determine_probabilities_clustering(dataset_filled,
                                                refresh = args$refresh,
                                                seed = args$seed)
} else if (model_type == "bayesian") {
  cat("[DINEMITES/run] Generating Bayesian infection/lag covariates...\n")
  lag_column <- paste0("lag_", args$bayesian_lag_days)
  lag_infection_column <- paste0("lag_infection_", args$bayesian_lag_days)
  dataset_filled <- dataset_filled %>%
    add_present_infection() %>%
    add_persistent_column() %>%
    add_persistent_infection() %>%
    add_lag_column(lag_time = args$bayesian_lag_days) %>%
    add_lag_infection(lag_time = args$bayesian_lag_days)

  results <- determine_probabilities_bayesian(
    dataset_filled,
    infection_persistence_covariates = c("persistent_infection", lag_infection_column),
    infection_general_covariates = NULL,
    alleles_persistence_covariates = c("persistent", lag_column),
    chains = args$bayesian_chains,
    parallel_chains = args$bayesian_parallel_chains,
    iter_warmup = args$bayesian_iter_warmup,
    iter_sampling = args$bayesian_iter_sampling,
    refresh = args$refresh,
    adapt_delta = args$bayesian_adapt_delta,
    seed = args$seed,
    drop_out = args$bayesian_drop_out
  )
}

t_elapsed <- difftime(Sys.time(), t_start, units = "secs")
cat("[DINEMITES/run] Model completed in", round(as.numeric(t_elapsed), 1), "seconds\n")

# ── 4. Add probability column to dataset ───────────────────────────────────
# determine_probabilities_* returns a list with $probability_new (vector)
# and $fit (model object or NULL). We attach the probabilities to the dataset.
cat("[DINEMITES/run] Attaching probabilities to dataset...\n")
dataset_filled$probability_new <- results$probability_new
allele_key <- build_allele_key(dataset_filled)
if (nrow(allele_key) > 0 && "locus" %in% colnames(dataset_filled)) {
  dataset_filled <- dataset_filled %>%
    left_join(allele_key, by = c("locus", "allele"))
} else {
  dataset_filled$short_allele_id <- NA_character_
}

# ── 5. Calculate new infections and molFOI ─────────────────────────────────
cat("[DINEMITES/run] Estimating new infections...\n")
estimated_new_infections_for_plot <- estimate_new_infections(dataset_filled)
new_infections <- estimated_new_infections_for_plot
if (!"subject" %in% colnames(new_infections)) {
  new_infections$subject <- rownames(new_infections)
  new_infections <- new_infections %>%
    select(subject, everything())
}
rownames(new_infections) <- NULL

cat("[DINEMITES/run] Calculating molecular FOI...\n")
molfoi <- compute_molFOI(dataset_filled)
if (!"subject" %in% colnames(molfoi)) {
  molfoi$subject <- rownames(molfoi)
  molfoi <- molfoi %>%
    select(subject, everything())
}
rownames(molfoi) <- NULL
subjects <- unique(dataset_filled$subject)

# ── 6. Write outputs ──────────────────────────────────────────────────────
dir.create(args$outdir, showWarnings = FALSE, recursive = TRUE)
plots_dir <- file.path(args$outdir, "dinemites_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)
subjects_dir <- file.path(args$outdir, "dinemites_subjects")
dir.create(subjects_dir, showWarnings = FALSE, recursive = TRUE)

# Allele probabilities
allele_probs_path <- file.path(args$outdir, "dinemites_allele_probabilities.tsv")
write.table(dataset_filled, file = allele_probs_path, sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("[DINEMITES/run] Wrote allele probabilities:", allele_probs_path, "\n")

# Allele key for static plot row IDs
allele_key_path <- file.path(args$outdir, "dinemites_allele_key.tsv")
write.table(allele_key, file = allele_key_path, sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("[DINEMITES/run] Wrote allele key:", allele_key_path, "\n")

# New infections per subject
new_inf_path <- file.path(args$outdir, "dinemites_new_infections.tsv")
write.table(new_infections, file = new_inf_path, sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("[DINEMITES/run] Wrote new infections:", new_inf_path, "\n")

# molFOI per subject
molfoi_path <- file.path(args$outdir, "dinemites_molfoi.tsv")
write.table(molfoi, file = molfoi_path, sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("[DINEMITES/run] Wrote molFOI:", molfoi_path, "\n")

# Per-subject outputs for patient-facing review.
for (subj in subjects) {
  subject_dir <- file.path(subjects_dir, safe_subject_dir_name(subj))
  dir.create(subject_dir, showWarnings = FALSE, recursive = TRUE)

  subject_dataset_all <- dataset_filled %>%
    filter(.data$subject == subj)
  subject_observed_alleles <- subject_dataset_all %>%
    filter(!is.na(.data$allele), trimws(as.character(.data$allele)) != "",
           .data$present > 0) %>%
    distinct(.data$locus, .data$allele)

  subject_dataset <- subject_dataset_all[0, ]
  if (nrow(subject_observed_alleles) > 0 && "locus" %in% colnames(subject_dataset_all)) {
    subject_dataset <- subject_dataset_all %>%
      semi_join(subject_observed_alleles, by = c("locus", "allele"))
  }

  subject_key <- allele_key
  if (nrow(subject_observed_alleles) > 0 && nrow(allele_key) > 0) {
    subject_key <- allele_key %>%
      semi_join(subject_observed_alleles, by = c("locus", "allele"))
  } else if (nrow(allele_key) > 0) {
    subject_key <- allele_key[0, ]
  }

  write.table(subject_dataset,
              file = file.path(subject_dir, "dinemites_allele_probabilities.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(subject_key,
              file = file.path(subject_dir, "dinemites_allele_key.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(new_infections %>% filter(.data$subject == subj),
              file = file.path(subject_dir, "dinemites_new_infections.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)
  write.table(molfoi %>% filter(.data$subject == subj),
              file = file.path(subject_dir, "dinemites_molfoi.tsv"),
              sep = "\t", row.names = FALSE, quote = FALSE)

  subject_summary_lines <- c(
    "{",
    paste0('  "model_type": "', model_type, '",'),
    paste0('  "subject": "', subj, '",'),
    paste0('  "n_alleles": ', nrow(subject_observed_alleles), ','),
    paste0('  "n_timepoints": ', length(unique(subject_dataset_all$time)), ','),
    paste0('  "plot": "subject_', subj, '.png"'),
    "}"
  )
  writeLines(subject_summary_lines, file.path(subject_dir, "dinemites_summary.json"))
  cat("[DINEMITES/run] Wrote subject outputs:", subject_dir, "\n")
}

# ── 7. Generate per-subject plots ──────────────────────────────────────────
cat("[DINEMITES/run] Generating per-subject plots...\n")
empty_treatments <- data.frame(subject = character(), time = numeric())
x_breaks <- sort(unique(dataset_filled$time))
x_axis_labels <- build_time_axis_labels(dataset_filled)
x_labels <- unname(x_axis_labels[as.character(x_breaks)])
allele_labeler <- build_allele_axis_labeler(allele_key)
x_limit_max <- ifelse(max(x_breaks) == min(x_breaks), max(x_breaks) + 1, max(x_breaks) * 1.1)
n_present_alleles <- length(unique(dataset_filled$allele[dataset_filled$present > 0]))
n_loci <- if ("locus" %in% colnames(dataset_filled)) {
  length(unique(dataset_filled$locus))
} else {
  1
}
plot_width <- max(args$plot_width, 14)
plot_height <- if (args$plot_height > 0) {
  args$plot_height
} else {
  min(28, max(8, 3.5 + (0.20 * n_present_alleles) + (0.45 * n_loci)))
}

cat("[DINEMITES/run] Plot size:", plot_width, "x", round(plot_height, 1), "inches\n")

tryCatch({
  plot_list <- plot_dataset(dataset_filled,
                            treatments = empty_treatments,
                            estimated_new_infections = estimated_new_infections_for_plot,
                            output = NULL,
                            height = plot_height,
                            width = plot_width)
  for (subj in subjects) {
    plot_path <- file.path(plots_dir, paste0("subject_", subj, ".png"))
    plot_out <- plot_list[[as.character(subj)]]
    if (!is.null(plot_out)) {
      plot_out <- apply_time_axis_to_plot(plot_out, x_breaks, x_labels,
                                          x_limit_max, allele_labeler)
      ggsave(plot_path, plot = plot_out, height = plot_height,
             width = plot_width, limitsize = FALSE)
      cat("[DINEMITES/run] Plot:", plot_path, "\n")
      subject_plot_path <- file.path(subjects_dir, safe_subject_dir_name(subj),
                                     paste0("subject_", subj, ".png"))
      invisible(file.copy(plot_path, subject_plot_path, overwrite = TRUE))
    }
  }
}, error = function(e) {
  cat("[DINEMITES/run] WARNING: Could not generate DINEMITES plots:",
      conditionMessage(e), "\n")
})

# ── 8. Summary report ─────────────────────────────────────────────────────
summary_path <- file.path(args$outdir, "dinemites_summary.json")
summary_data <- list(
  model_type     = model_type,
  seed           = args$seed,
  refresh        = args$refresh,
  bayesian_lag_days = args$bayesian_lag_days,
  bayesian_chains = args$bayesian_chains,
  bayesian_parallel_chains = args$bayesian_parallel_chains,
  bayesian_iter_warmup = args$bayesian_iter_warmup,
  bayesian_iter_sampling = args$bayesian_iter_sampling,
  bayesian_adapt_delta = args$bayesian_adapt_delta,
  bayesian_drop_out = args$bayesian_drop_out,
  n_subjects     = length(subjects),
  n_alleles      = length(unique(dataset_filled$allele)),
  n_timepoints   = length(unique(dataset_filled$time)),
  runtime_secs   = round(as.numeric(t_elapsed), 1),
  subjects       = subjects,
  plots          = paste0("subject_", subjects, ".png")
)

# Write JSON manually (avoid jsonlite dependency)
json_lines <- c("{")
json_lines <- c(json_lines,
  paste0('  "model_type": "', summary_data$model_type, '",'),
  paste0('  "seed": ', summary_data$seed, ','),
  paste0('  "refresh": ', summary_data$refresh, ','),
  paste0('  "bayesian_lag_days": ', summary_data$bayesian_lag_days, ','),
  paste0('  "bayesian_chains": ', summary_data$bayesian_chains, ','),
  paste0('  "bayesian_parallel_chains": ', summary_data$bayesian_parallel_chains, ','),
  paste0('  "bayesian_iter_warmup": ', summary_data$bayesian_iter_warmup, ','),
  paste0('  "bayesian_iter_sampling": ', summary_data$bayesian_iter_sampling, ','),
  paste0('  "bayesian_adapt_delta": ', summary_data$bayesian_adapt_delta, ','),
  paste0('  "bayesian_drop_out": ', tolower(as.character(summary_data$bayesian_drop_out)), ','),
  paste0('  "n_subjects": ', summary_data$n_subjects, ','),
  paste0('  "n_alleles": ', summary_data$n_alleles, ','),
  paste0('  "n_timepoints": ', summary_data$n_timepoints, ','),
  paste0('  "runtime_secs": ', summary_data$runtime_secs, ','),
  paste0('  "subjects": [', paste0('"', summary_data$subjects, '"', collapse = ', '), '],'),
  paste0('  "plots": [', paste0('"', summary_data$plots, '"', collapse = ', '), ']')
)
json_lines <- c(json_lines, "}")
writeLines(json_lines, summary_path)

cat("[DINEMITES/run] Wrote summary:", summary_path, "\n")
cat("[DINEMITES/run] Done.\n")
