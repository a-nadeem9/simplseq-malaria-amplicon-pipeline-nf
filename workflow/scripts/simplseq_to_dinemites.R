#!/usr/bin/env Rscript
# simplseq_to_dinemites.R
# ---------------------------------------------------------------------------
# Bridge between SIMPLseq CIGAR output and DINEMITES longitudinal input format.
#
# Reads the wide seqtab_cigar.tsv and samples.csv, performs replicate
# intersection merging (keep only haplotypes present in BOTH replicates),
# converts dates (YYYY-MM → YYYY-MM-27 → days since earliest), and outputs
# a tab-delimited file with columns: allele, time, subject, locus.
# ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
})

DEFAULT_COLLECTION_YEAR <- "2022"
MONTH_ALIASES <- c(
  jan = "01", january = "01",
  feb = "02", february = "02",
  mar = "03", march = "03",
  apr = "04", april = "04",
  may = "05",
  jun = "06", june = "06",
  jul = "07", july = "07",
  aug = "08", august = "08",
  sep = "09", sept = "09", september = "09",
  oct = "10", october = "10",
  nov = "11", november = "11",
  dec = "12", december = "12"
)
MONTH_PATTERN <- paste(names(MONTH_ALIASES)[order(nchar(names(MONTH_ALIASES)), decreasing = TRUE)],
                       collapse = "|")

infer_collection_date_from_sample_id <- function(sample_id) {
  label <- tolower(as.character(sample_id))
  pattern <- paste0("(^|[^[:alnum:]])(", MONTH_PATTERN, ")([^[:alnum:]]|$)")
  match <- regexec(pattern, label, ignore.case = TRUE, perl = TRUE)
  parts <- regmatches(label, match)
  vapply(parts, function(item) {
    if (length(item) < 3) {
      return("")
    }
    month <- MONTH_ALIASES[[tolower(item[[3]])]]
    if (is.null(month) || is.na(month)) {
      return("")
    }
    paste0(DEFAULT_COLLECTION_YEAR, "-", month)
  }, character(1))
}

# ── CLI arguments ──────────────────────────────────────────────────────────
option_list <- list(
  make_option("--cigar", type = "character", help = "Path to seqtab_cigar.tsv"),
  make_option("--samples", type = "character", help = "Path to samples.csv"),
  make_option("--out", type = "character", help = "Output path for DINEMITES input TSV"),
  make_option("--min_abundance_pct", type = "double", default = 0.3,
              help = "Minimum allele abundance within each sequencing sample [default %default]"),
  make_option("--abundance_denominator", type = "character", default = "locus",
              help = "Allele abundance denominator: locus or sample [default %default]"),
  make_option("--day_of_month", type = "integer", default = 27,
              help = "Day to assume for YYYY-MM dates [default %default]")
)
args <- parse_args(OptionParser(option_list = option_list))

if (is.null(args$cigar) || is.null(args$samples) || is.null(args$out)) {
  stop("Required arguments: --cigar, --samples, --out")
}
if (is.na(args$min_abundance_pct) || args$min_abundance_pct < 0 || args$min_abundance_pct > 100) {
  stop("--min_abundance_pct must be between 0 and 100. Got: ", args$min_abundance_pct)
}
args$abundance_denominator <- tolower(trimws(as.character(args$abundance_denominator)))
if (!args$abundance_denominator %in% c("locus", "sample")) {
  stop("--abundance_denominator must be either locus or sample. Got: ",
       args$abundance_denominator)
}

cat("[DINEMITES/bridge] Reading CIGAR table:", args$cigar, "\n")
cat("[DINEMITES/bridge] Reading samples:", args$samples, "\n")
cat("[DINEMITES/bridge] Allele abundance filter:",
    args$min_abundance_pct, "% of",
    ifelse(args$abundance_denominator == "locus", "sample+locus reads", "total sample reads"),
    "\n")

# ── 1. Read inputs ─────────────────────────────────────────────────────────
cigar_wide <- read.delim(args$cigar, header = TRUE, sep = "\t",
                         check.names = FALSE, stringsAsFactors = FALSE)

samples <- read.csv(args$samples, header = TRUE, stringsAsFactors = FALSE)

# Normalise column names for samples.csv (handle minor variations)
colnames(samples) <- tolower(trimws(colnames(samples)))

required_sample_cols <- c("sample_id", "participant_id", "collection_date", "replicate")
missing_cols <- setdiff(required_sample_cols, colnames(samples))
if (length(missing_cols) > 0) {
  stop("[DINEMITES/bridge] ERROR: Missing required samples.csv columns: ",
       paste(missing_cols, collapse = ", "))
}

# ── 2. Parse column headers: "LOCUS,CIGAR" → locus + cigar ────────────────
# The first column is "sample"; remaining are "LOCUS,CIGAR" encoded haplotypes
haplotype_cols <- setdiff(colnames(cigar_wide), "sample")

parsed_cols <- data.frame(
  col_name = haplotype_cols,
  locus    = sub(",.*", "", haplotype_cols),
  cigar    = sub("^[^,]*,", "", haplotype_cols),
  stringsAsFactors = FALSE
)

cat("[DINEMITES/bridge]", nrow(parsed_cols), "haplotype columns across",
    length(unique(parsed_cols$locus)), "loci\n")

# ── 3. Pivot wide → long ──────────────────────────────────────────────────
cigar_long <- cigar_wide %>%
  pivot_longer(
    cols      = all_of(haplotype_cols),
    names_to  = "haplotype_col",
    values_to = "reads"
  ) %>%
  left_join(parsed_cols, by = c("haplotype_col" = "col_name")) %>%
  mutate(reads = suppressWarnings(as.numeric(.data$reads))) %>%
  mutate(reads = ifelse(is.na(.data$reads), 0, .data$reads)) %>%
  select(sample_id = sample, locus, cigar, reads)

# ── 4. Join metadata ──────────────────────────────────────────────────────
# Match samples.csv on sample_id
sample_type <- if ("sample_type" %in% colnames(samples)) {
  tolower(trimws(samples$sample_type))
} else {
  rep("sample", nrow(samples))
}

meta_candidates <- samples %>%
  mutate(
    .sample_type_for_filter = sample_type,
    sample_id = ifelse(is.na(.data$sample_id), "", trimws(as.character(.data$sample_id))),
    participant_id = ifelse(is.na(.data$participant_id), "", trimws(as.character(.data$participant_id))),
    collection_date = ifelse(is.na(.data$collection_date), "", trimws(as.character(.data$collection_date))),
    replicate = ifelse(is.na(.data$replicate), "", trimws(as.character(.data$replicate)))
  ) %>%
  filter(.data$.sample_type_for_filter != "negative" | is.na(.data$.sample_type_for_filter))

matched_candidates <- meta_candidates %>%
  semi_join(cigar_long %>% distinct(sample_id), by = "sample_id")

if (nrow(matched_candidates) < nrow(meta_candidates)) {
  cat("[DINEMITES/bridge] WARNING:", nrow(meta_candidates) - nrow(matched_candidates),
      "sample sheet rows are not present in the CIGAR table and will be ignored.\n")
}

if (nrow(matched_candidates) == 0) {
  stop("[DINEMITES/bridge] ERROR: No sample_id values matched between the CIGAR table and samples.csv.")
}

missing_participant <- sum(nchar(matched_candidates$participant_id) == 0)
missing_date <- sum(nchar(matched_candidates$collection_date) == 0)
if (missing_participant > 0) {
  stop("[DINEMITES/bridge] ERROR: DINEMITES requires participant_id for every matched sample. ",
       missing_participant, " matched rows are missing participant_id.")
}
if (missing_date > 0) {
  inferred_dates <- infer_collection_date_from_sample_id(matched_candidates$sample_id)
  can_infer <- nchar(matched_candidates$collection_date) == 0 & nchar(inferred_dates) > 0
  if (any(can_infer)) {
    matched_candidates$collection_date[can_infer] <- inferred_dates[can_infer]
    cat("[DINEMITES/bridge] WARNING:", sum(can_infer),
        "matched rows had month but no year in sample_id; assuming year",
        DEFAULT_COLLECTION_YEAR, "for DINEMITES only.\n")
  }
  missing_date <- sum(nchar(matched_candidates$collection_date) == 0)
  if (missing_date > 0) {
    stop("[DINEMITES/bridge] ERROR: DINEMITES requires collection_date values in YYYY-MM format. ",
         missing_date, " matched sample rows are missing collection_date and no month could be inferred from sample_id.")
  }
}

meta <- matched_candidates %>%
  select(sample_id, participant_id, collection_date, replicate)

matched_meta <- meta

sample_timepoints <- matched_meta %>%
  distinct(participant_id, collection_date) %>%
  mutate(
    date_full = as.Date(paste0(collection_date, "-",
                               sprintf("%02d", args$day_of_month)),
                        format = "%Y-%m-%d")
  )

if (any(is.na(sample_timepoints$date_full))) {
  stop("[DINEMITES/bridge] ERROR: Could not parse one or more collection_date values. ",
       "Expected YYYY-MM.")
}

sample_timepoints <- sample_timepoints %>%
  mutate(
    time = as.integer(date_full - min(date_full)),
    date_label = format(date_full, "%b %Y")
  )

cigar_meta <- cigar_long %>%
  inner_join(matched_meta, by = "sample_id")

if (nrow(cigar_meta) == 0) {
  stop("[DINEMITES/bridge] ERROR: No samples matched between CIGAR table and samples.csv. ",
       "Check that sample IDs match.")
}

cat("[DINEMITES/bridge]", nrow(cigar_meta), "rows after metadata join (",
    length(unique(cigar_meta$participant_id)), " participants)\n")

# ── 5. Determine presence (reads > 0) ─────────────────────────────────────
cigar_meta <- cigar_meta %>%
  group_by(.data$sample_id, dplyr::across(dplyr::all_of(
    if (args$abundance_denominator == "locus") "locus" else character()
  ))) %>%
  mutate(total_sample_reads = sum(.data$reads, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    min_reads_required = ifelse(
      .data$total_sample_reads > 0 & args$min_abundance_pct > 0,
      pmax(1, ceiling(.data$total_sample_reads * (args$min_abundance_pct / 100))),
      1
    ),
    allele_abundance_pct = ifelse(
      .data$total_sample_reads > 0,
      100 * .data$reads / .data$total_sample_reads,
      0
    ),
    present = as.integer(.data$reads > 0 & .data$reads >= min_reads_required)
  )

# ── 6. Replicate intersection merging ──────────────────────────────────────
# For each participant + date + locus + cigar combination:
# Keep haplotype ONLY if present in ALL replicates for that time point.
# If only one replicate exists, keep as-is.

merged <- cigar_meta %>%
  group_by(participant_id, collection_date, locus, cigar) %>%
  summarise(
    n_replicates     = n(),
    n_present        = sum(present),
    min_reads_required = max(min_reads_required, na.rm = TRUE),
    max_abundance_pct = max(allele_abundance_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # Intersection: keep only if present in ALL replicates
  filter(n_present == n_replicates & n_present > 0)

cat("[DINEMITES/bridge]", nrow(merged),
    "haplotype-timepoints after replicate intersection merging\n")

if (nrow(merged) == 0) {
  cat("[DINEMITES/bridge] WARNING: No haplotypes survived replicate intersection. ",
      "Writing sampled visits as allele=NA.\n", sep = "")
}

# ── 7. Date conversion: YYYY-MM → YYYY-MM-27 → days since earliest ────────
merged <- merged %>%
  left_join(sample_timepoints,
            by = c("participant_id", "collection_date"))

# Report date range
date_range <- range(sample_timepoints$date_full)
cat("[DINEMITES/bridge] Date range:",
    format(date_range[1]), "to", format(date_range[2]),
    "(", max(sample_timepoints$time), "days )\n")

# ── 8. Build DINEMITES input ──────────────────────────────────────────────
# Allele naming: LOCUS:CIGAR (e.g., "KELT:139C", "SERA8:.")
present_input <- merged %>%
  mutate(
    allele  = paste0(locus, ":", cigar),
    subject = participant_id
  ) %>%
  select(allele, time, subject, locus, collection_date, date_full, date_label)

empty_timepoints <- sample_timepoints %>%
  anti_join(merged %>% distinct(participant_id, collection_date),
            by = c("participant_id", "collection_date")) %>%
  transmute(
    allele = NA_character_,
    time,
    subject = participant_id,
    locus = NA_character_,
    collection_date,
    date_full,
    date_label
  )

dinemites_input <- bind_rows(present_input, empty_timepoints) %>%
  arrange(subject, time, locus, allele)

# Report summary
n_subjects  <- length(unique(dinemites_input$subject))
n_timepoints <- length(unique(dinemites_input$time))
n_alleles   <- length(unique(dinemites_input$allele[!is.na(dinemites_input$allele)]))
n_loci      <- length(unique(dinemites_input$locus[!is.na(dinemites_input$locus)]))

cat("[DINEMITES/bridge] Output summary:\n")
cat("  Subjects:    ", n_subjects, "\n")
cat("  Time points: ", n_timepoints, "\n")
cat("  Unique alleles:", n_alleles, "\n")
cat("  Loci:        ", n_loci, "\n")
cat("  Total rows:  ", nrow(dinemites_input), "\n")

# ── 9. Write output ───────────────────────────────────────────────────────
dir.create(dirname(args$out), showWarnings = FALSE, recursive = TRUE)
write.table(dinemites_input, file = args$out, sep = "\t",
            row.names = FALSE, quote = FALSE)

cat("[DINEMITES/bridge] Wrote", args$out, "\n")
cat("[DINEMITES/bridge] Done.\n")
