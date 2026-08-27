#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(facets))

parse_args <- function(args) {
  out <- list()
  idx <- 1
  while (idx <= length(args)) {
    key <- args[[idx]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
    key <- substring(key, 3)
    if (idx == length(args) || startsWith(args[[idx + 1]], "--")) {
      out[[key]] <- TRUE
      idx <- idx + 1
    } else {
      out[[key]] <- args[[idx + 1]]
      idx <- idx + 2
    }
  }
  out
}

require_arg <- function(args, name) {
  value <- args[[name]]
  if (is.null(value) || identical(value, TRUE)) {
    stop(sprintf("Missing required argument --%s", name))
  }
  value
}

first_existing <- function(df, candidates) {
  for (candidate in candidates) {
    if (candidate %in% colnames(df)) {
      return(candidate)
    }
  }
  stop(sprintf("None of the expected columns were found: %s", paste(candidates, collapse = ", ")))
}

normalize_gbuild <- function(build) {
  if (build == "GRCh38") {
    return("hg38")
  }
  if (build == "GRCh37") {
    return("hg19")
  }
  build
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

counts_file <- require_arg(args, "counts-file")
sample_id <- require_arg(args, "sample-id")
genome_build <- normalize_gbuild(require_arg(args, "genome-build"))
output_purity <- require_arg(args, "output-purity")
output_allele_specific <- require_arg(args, "output-allele-specific")
output_rds <- require_arg(args, "output-rds")

cval <- as.numeric(if (!is.null(args[["cval"]])) args[["cval"]] else 150)
ndepth <- as.numeric(if (!is.null(args[["ndepth"]])) args[["ndepth"]] else 35)
ndepthmax <- as.numeric(if (!is.null(args[["ndepthmax"]])) args[["ndepthmax"]] else 1000)

rcmat <- readSnpMatrix(filename = counts_file)
xx <- preProcSample(rcmat, ndepth = ndepth, ndepthmax = ndepthmax, gbuild = genome_build)
oo <- procSample(xx, cval = cval)
fit <- emcncf(oo)

cncf <- fit$cncf
chrom_col <- first_existing(cncf, c("chrom", "Chrom"))
start_col <- first_existing(cncf, c("start", "loc.start"))
end_col <- first_existing(cncf, c("end", "loc.end"))
log2_col <- first_existing(cncf, c("cnlr.median", "segmean"))
tcn_col <- first_existing(cncf, c("tcn.em", "tcn"))
lcn_col <- first_existing(cncf, c("lcn.em", "lcn"))
cf_col <- first_existing(cncf, c("cf.em", "cf"))

major_cn <- cncf[[tcn_col]] - cncf[[lcn_col]]
minor_cn <- cncf[[lcn_col]]
loh <- ifelse(!is.na(minor_cn) & minor_cn == 0, "true", "false")

purity_row <- data.frame(
  patient_id = sample_id,
  purity = fit$purity,
  ploidy = fit$ploidy,
  diplogr = fit$dipLogR,
  stringsAsFactors = FALSE
)
write.table(purity_row, file = output_purity, sep = "\t", row.names = FALSE, quote = FALSE)

allele_specific <- data.frame(
  patient_id = sample_id,
  chromosome = cncf[[chrom_col]],
  start = cncf[[start_col]],
  end = cncf[[end_col]],
  log2 = cncf[[log2_col]],
  total_cn = cncf[[tcn_col]],
  major_cn = major_cn,
  minor_cn = minor_cn,
  segment_cf = cncf[[cf_col]],
  loh = loh,
  stringsAsFactors = FALSE
)
write.table(allele_specific, file = output_allele_specific, sep = "\t", row.names = FALSE, quote = FALSE)

saveRDS(fit, file = output_rds)
