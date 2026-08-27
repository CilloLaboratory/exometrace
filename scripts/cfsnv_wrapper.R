#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(cfSNV))

parse_args <- function(argv) {
  if (length(argv) == 0) {
    stop("missing subcommand")
  }
  command <- argv[[1]]
  rest <- argv[-1]
  args <- list()
  i <- 1
  while (i <= length(rest)) {
    token <- rest[[i]]
    if (!startsWith(token, "--")) {
      stop(sprintf("unexpected positional argument '%s'", token))
    }
    key <- gsub("-", "_", substring(token, 3))
    if (i == length(rest) || startsWith(rest[[i + 1]], "--")) {
      args[[key]] <- TRUE
      i <- i + 1
    } else {
      args[[key]] <- rest[[i + 1]]
      i <- i + 2
    }
  }
  list(command = command, args = args)
}

require_arg <- function(args, key) {
  value <- args[[key]]
  if (is.null(value) || identical(value, TRUE) || identical(value, "")) {
    stop(sprintf("missing required argument --%s", gsub("_", "-", key)))
  }
  value
}

tool_path <- function(..., env = NULL) {
  if (!is.null(env)) {
    env_value <- Sys.getenv(env, "")
    if (nzchar(env_value)) {
      return(env_value)
    }
  }
  for (candidate in c(...)) {
    path <- Sys.which(candidate)
    if (nzchar(path)) {
      return(path)
    }
  }
  stop(sprintf("required tool not found: %s", paste(c(...), collapse = ", ")))
}

read_blocked_positions <- function(path) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    return(character())
  }
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  rows <- lines[!startsWith(lines, "#") & nzchar(lines)]
  if (!length(rows)) {
    return(character())
  }
  fields <- strsplit(rows, "\t", fixed = TRUE)
  refs <- vapply(fields, function(x) if (length(x) >= 4) x[[4]] else "", character(1))
  alts <- vapply(fields, function(x) if (length(x) >= 5) x[[5]] else "", character(1))
  chroms <- vapply(fields, `[[`, character(1), 1)
  poss <- vapply(fields, `[[`, character(1), 2)
  paste(chroms, poss, refs, alts, sep = ":")
}

filter_blocked_positions <- function(variant_list, blocked_path) {
  blocked <- read_blocked_positions(blocked_path)
  if (!length(blocked) || !nrow(variant_list)) {
    return(variant_list)
  }
  keys <- paste(variant_list$CHROM, variant_list$POS, variant_list$REF, variant_list$ALT, sep = ":")
  variant_list[!(keys %in% blocked), , drop = FALSE]
}

write_vcf <- function(variant_list, output_path, tumor_sample, normal_sample) {
  tmp_vcf <- sub("\\.gz$", "", output_path)
  if (identical(tmp_vcf, output_path)) {
    stop("output path must end with .gz")
  }

  info_lines <- c(
    "##fileformat=VCFv4.2",
    "##source=cfSNV_R_wrapper",
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth placeholder\">",
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths placeholder\">"
  )
  header <- paste(
    c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", tumor_sample, normal_sample),
    collapse = "\t"
  )

  con <- file(tmp_vcf, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(info_lines, header), con)
  if (nrow(variant_list)) {
    rows <- apply(variant_list, 1, function(row) {
      paste(c(row[["CHROM"]], row[["POS"]], row[["ID"]], row[["REF"]], row[["ALT"]],
              row[["QUAL"]], row[["FILTER"]], row[["INFO"]], "GT:DP:AD", "./.:0:0,0", "./.:0:0,0"),
            collapse = "\t")
    })
    writeLines(rows, con)
  }
  close(con)
  bgzip <- tool_path("bgzip")
  tabix <- tool_path("tabix")
  system2(bgzip, c("-f", tmp_vcf))
  system2(tabix, c("-f", "-p", "vcf", output_path))
}

run_stdprep <- function(args) {
  output_dir <- normalizePath(require_arg(args, "output_dir"), mustWork = FALSE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cfSNV::getbam_align(
    fastq1 = require_arg(args, "fastq1"),
    fastq2 = require_arg(args, "fastq2"),
    reference = require_arg(args, "reference"),
    SNP.database = require_arg(args, "snp_database"),
    samtools.dir = tool_path("samtools"),
    picard.dir = tool_path(env = "CFSNV_PICARD_JAR"),
    bedtools.dir = tool_path("bedtools"),
    GATK.dir = tool_path(env = "CFSNV_GATK_JAR"),
    bwa.dir = tool_path("bwa"),
    sample.id = require_arg(args, "sample_id"),
    output.dir = output_dir,
    java.dir = tool_path("java")
  )
}

run_cfdnaprep <- function(args) {
  output_dir <- normalizePath(require_arg(args, "output_dir"), mustWork = FALSE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  cfSNV::getbam_align_after_merge(
    fastq1 = require_arg(args, "fastq1"),
    fastq2 = require_arg(args, "fastq2"),
    reference = require_arg(args, "reference"),
    SNP.database = require_arg(args, "snp_database"),
    samtools.dir = tool_path("samtools"),
    picard.dir = tool_path(env = "CFSNV_PICARD_JAR"),
    bedtools.dir = tool_path("bedtools"),
    GATK.dir = tool_path(env = "CFSNV_GATK_JAR"),
    bwa.dir = tool_path("bwa"),
    flash.dir = tool_path("flash2", "flash"),
    sample.id = require_arg(args, "sample_id"),
    output.dir = output_dir,
    java.dir = tool_path("java")
  )
}

run_detectmuts <- function(args) {
  results <- cfSNV::variant_calling(
    plasma.unmerged = require_arg(args, "tumor_bam"),
    normal = require_arg(args, "normal_bam"),
    plasma.merged.extendedFrags = require_arg(args, "extended_bam"),
    plasma.merge.notCombined = require_arg(args, "not_combined_bam"),
    target.bed = require_arg(args, "targets"),
    reference = require_arg(args, "reference"),
    SNP.database = require_arg(args, "snp_database"),
    samtools.dir = tool_path("samtools"),
    picard.dir = tool_path(env = "CFSNV_PICARD_JAR"),
    bedtools.dir = tool_path("bedtools"),
    sample.id = require_arg(args, "sample_id"),
    MIN_HOLD_SUPPORT_COUNT = as.integer(require_arg(args, "min_hold_support")),
    MIN_PASS_SUPPORT_COUNT = as.integer(require_arg(args, "min_pass_support")),
    java.dir = tool_path("java"),
    python.dir = tool_path("python3", "python")
  )
  variant_list <- filter_blocked_positions(results$variant.list, args$blocked_positions)
  write_vcf(
    variant_list = variant_list,
    output_path = require_arg(args, "output"),
    tumor_sample = require_arg(args, "tumor_sample"),
    normal_sample = require_arg(args, "normal_sample")
  )
}

main <- function() {
  parsed <- parse_args(commandArgs(trailingOnly = TRUE))
  command <- parsed$command
  args <- parsed$args
  switch(
    command,
    STDprep = run_stdprep(args),
    cfDNAprep = run_cfdnaprep(args),
    DetectMuts = run_detectmuts(args),
    stop(sprintf("unsupported subcommand '%s'", command))
  )
}

main()
