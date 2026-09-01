#!/usr/bin/env Rscript

ensure_writable_cfsnv_library <- function() {
  lib_root <- Sys.getenv("CFSNV_R_LIB_ROOT", "")
  if (!nzchar(lib_root)) {
    return(invisible(NULL))
  }

  dir.create(lib_root, showWarnings = FALSE, recursive = TRUE)
  source_pkg <- find.package("cfSNV", quiet = TRUE)
  if (!length(source_pkg)) {
    stop("cfSNV package is not installed")
  }

  source_pkg <- source_pkg[[1]]
  target_pkg <- file.path(lib_root, basename(source_pkg))
  if (!dir.exists(target_pkg)) {
    ok <- file.copy(source_pkg, lib_root, recursive = TRUE)
    if (!ok) {
      stop(sprintf("failed to copy cfSNV package from %s to %s", source_pkg, lib_root))
    }
  }

  .libPaths(unique(c(lib_root, .libPaths())))
  suppressPackageStartupMessages(library("cfSNV", lib.loc = lib_root, character.only = TRUE))
  invisible(target_pkg)
}

configure_cfsnv_tmpdir <- function(package_root, tmpdir) {
  if (is.null(package_root) || !nzchar(package_root)) {
    return(invisible(NULL))
  }

  if (.Platform$OS.type == "windows") {
    stop("cfSNV wrapper requires Unix-style symlink support for tmpdir configuration")
  }

  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
  extdata_dir <- file.path(package_root, "extdata")
  dir.create(extdata_dir, showWarnings = FALSE, recursive = TRUE)
  package_tmpdir <- file.path(extdata_dir, "tmp")
  if (file.exists(package_tmpdir) || dir.exists(package_tmpdir) || nzchar(Sys.readlink(package_tmpdir))) {
    unlink(package_tmpdir, recursive = TRUE, force = TRUE)
  }

  linked <- tryCatch(
    isTRUE(file.symlink(tmpdir, package_tmpdir)),
    warning = function(...) FALSE,
    error = function(...) FALSE
  )
  resolved_link <- Sys.readlink(package_tmpdir)
  if (!linked || !nzchar(resolved_link) || normalizePath(resolved_link, mustWork = FALSE) != normalizePath(tmpdir, mustWork = FALSE)) {
    stop(sprintf("failed to bind cfSNV extdata tmpdir %s to %s", package_tmpdir, tmpdir))
  }

  invisible(package_tmpdir)
}

cfsnv_package_root <- ensure_writable_cfsnv_library()

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

tool_path <- function(..., env = NULL, default_paths = character()) {
  if (!is.null(env)) {
    env_value <- Sys.getenv(env, "")
    if (nzchar(env_value)) {
      return(env_value)
    }
  }
  for (candidate in default_paths) {
    if (nzchar(candidate) && file.exists(candidate)) {
      return(candidate)
    }
  }
  for (candidate in c(...)) {
    path <- Sys.which(candidate)
    if (nzchar(path)) {
      return(path)
    }
  }
  known_paths <- c(c(...), default_paths)
  stop(sprintf("required tool not found: %s", paste(known_paths[nzchar(known_paths)], collapse = ", ")))
}

require_existing_path <- function(path, label) {
  normalized <- normalizePath(path, mustWork = FALSE)
  if (!file.exists(normalized)) {
    stop(sprintf("%s does not exist: %s", label, normalized))
  }
  normalized
}

resolve_java_path <- function() {
  explicit <- Sys.getenv("CFSNV_JAVA", "")
  if (nzchar(explicit)) {
    return(explicit)
  }

  java_home <- Sys.getenv("JAVA_HOME", "")
  if (nzchar(java_home)) {
    candidate <- file.path(java_home, "bin", "java")
    if (file.exists(candidate)) {
      return(candidate)
    }
  }

  tool_path("java", default_paths = c("/opt/conda/bin/java", "/usr/bin/java"))
}

resolve_picard_path <- function() {
  tool_path(
    env = "CFSNV_PICARD_JAR",
    default_paths = c("/usr/local/share/cfsnv-tools/picard.jar")
  )
}

resolve_gatk3_path <- function() {
  tool_path(
    env = "CFSNV_GATK_JAR",
    default_paths = c("/usr/local/share/cfsnv-tools/GenomeAnalysisTK.jar", "/opt/gatk3/GenomeAnalysisTK.jar")
  )
}

stage_cfsnv_input_file <- function(source_path, target_path) {
  dir.create(dirname(target_path), showWarnings = FALSE, recursive = TRUE)
  if (file.exists(target_path) || nzchar(Sys.readlink(target_path))) {
    unlink(target_path, recursive = FALSE, force = TRUE)
  }

  linked <- tryCatch(
    isTRUE(file.symlink(source_path, target_path)),
    warning = function(...) FALSE,
    error = function(...) FALSE
  )
  if (!linked) {
    copied <- file.copy(source_path, target_path, overwrite = TRUE)
    if (!isTRUE(copied)) {
      stop(sprintf("failed to stage cfSNV input %s at %s", source_path, target_path))
    }
  }
}

stage_cfsnv_detectmuts_inputs <- function(tmpdir, sample_id, tumor_bam, normal_bam, extended_bam, not_combined_bam) {
  deeplearn_dir <- file.path(tmpdir, "deeplearn")
  tumor_bam <- require_existing_path(tumor_bam, "tumor BAM")
  normal_bam <- require_existing_path(normal_bam, "normal BAM")
  extended_bam <- require_existing_path(extended_bam, "extended-fragment BAM")
  not_combined_bam <- require_existing_path(not_combined_bam, "non-overlapping BAM")

  aliases <- list(
    tumor_bam = c(
      sprintf("%s.recal.bam", sample_id),
      sprintf("%s.paired-reads.bam", sample_id)
    ),
    normal_bam = c(
      sprintf("%s.normal.recal.bam", sample_id),
      sprintf("%s.normal-blood.bam", sample_id)
    ),
    extended_bam = c(
      sprintf("%s.extendedFrags.recal.bam", sample_id),
      sprintf("%s.paired-reads.extendedFrags.recal.bam", sample_id)
    ),
    not_combined_bam = c(
      sprintf("%s.notCombined.recal.bam", sample_id),
      sprintf("%s.paired-reads.notCombined.recal.bam", sample_id)
    )
  )
  inputs <- list(
    tumor_bam = tumor_bam,
    normal_bam = normal_bam,
    extended_bam = extended_bam,
    not_combined_bam = not_combined_bam
  )

  for (label in names(aliases)) {
    source_path <- inputs[[label]]
    for (alias_name in aliases[[label]]) {
      stage_cfsnv_input_file(source_path, file.path(deeplearn_dir, alias_name))
    }

    source_index <- paste0(source_path, ".bai")
    if (file.exists(source_index)) {
      for (alias_name in aliases[[label]]) {
        stage_cfsnv_input_file(source_index, file.path(deeplearn_dir, paste0(alias_name, ".bai")))
      }
    }
  }

  invisible(deeplearn_dir)
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
  tmpdir <- Sys.getenv("TMPDIR", unset = file.path(output_dir, "tmp"))
  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
  configure_cfsnv_tmpdir(cfsnv_package_root, tmpdir)
  cfSNV::getbam_align(
    fastq1 = require_arg(args, "fastq1"),
    fastq2 = require_arg(args, "fastq2"),
    reference = require_arg(args, "reference"),
    SNP.database = require_arg(args, "snp_database"),
    samtools.dir = tool_path("samtools"),
    picard.dir = resolve_picard_path(),
    bedtools.dir = tool_path("bedtools"),
    GATK.dir = resolve_gatk3_path(),
    bwa.dir = tool_path("bwa"),
    sample.id = require_arg(args, "sample_id"),
    output.dir = output_dir,
    java.dir = resolve_java_path()
  )
}

run_cfdnaprep <- function(args) {
  output_dir <- normalizePath(require_arg(args, "output_dir"), mustWork = FALSE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  tmpdir <- Sys.getenv("TMPDIR", unset = file.path(output_dir, "tmp"))
  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
  configure_cfsnv_tmpdir(cfsnv_package_root, tmpdir)
  cfSNV::getbam_align_after_merge(
    fastq1 = require_arg(args, "fastq1"),
    fastq2 = require_arg(args, "fastq2"),
    reference = require_arg(args, "reference"),
    SNP.database = require_arg(args, "snp_database"),
    samtools.dir = tool_path("samtools"),
    picard.dir = resolve_picard_path(),
    bedtools.dir = tool_path("bedtools"),
    GATK.dir = resolve_gatk3_path(),
    bwa.dir = tool_path("bwa"),
    flash.dir = tool_path("flash2", "flash"),
    sample.id = require_arg(args, "sample_id"),
    output.dir = output_dir,
    java.dir = resolve_java_path()
  )
}

run_detectmuts <- function(args) {
  output_dir <- dirname(normalizePath(require_arg(args, "output"), mustWork = FALSE))
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  tmpdir <- Sys.getenv("TMPDIR", unset = file.path(output_dir, "tmp"))
  dir.create(tmpdir, showWarnings = FALSE, recursive = TRUE)
  configure_cfsnv_tmpdir(cfsnv_package_root, tmpdir)
  sample_id <- require_arg(args, "sample_id")
  tumor_bam <- require_existing_path(require_arg(args, "tumor_bam"), "tumor BAM")
  normal_bam <- require_existing_path(require_arg(args, "normal_bam"), "normal BAM")
  extended_bam <- require_existing_path(require_arg(args, "extended_bam"), "extended-fragment BAM")
  not_combined_bam <- require_existing_path(require_arg(args, "not_combined_bam"), "non-overlapping BAM")
  stage_cfsnv_detectmuts_inputs(
    tmpdir = tmpdir,
    sample_id = sample_id,
    tumor_bam = tumor_bam,
    normal_bam = normal_bam,
    extended_bam = extended_bam,
    not_combined_bam = not_combined_bam
  )
  results <- cfSNV::variant_calling(
    plasma.unmerged = tumor_bam,
    normal = normal_bam,
    plasma.merged.extendedFrags = extended_bam,
    plasma.merge.notCombined = not_combined_bam,
    target.bed = require_arg(args, "targets"),
    reference = require_arg(args, "reference"),
    SNP.database = require_arg(args, "snp_database"),
    samtools.dir = tool_path("samtools"),
    picard.dir = resolve_picard_path(),
    bedtools.dir = tool_path("bedtools"),
    sample.id = sample_id,
    MIN_HOLD_SUPPORT_COUNT = as.integer(require_arg(args, "min_hold_support")),
    MIN_PASS_SUPPORT_COUNT = as.integer(require_arg(args, "min_pass_support")),
    java.dir = resolve_java_path(),
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
