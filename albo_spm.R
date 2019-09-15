#!/usr/bin/env Rscript

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# ALBOPICTUS SEX PREDICTION MODEL
#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

VERSION <- "0.99.1"

suppressPackageStartupMessages({
  suppressWarnings({
    library(optparse, quietly=TRUE)
    library(e1071, quietly=TRUE)
    library(ps, quietly=TRUE)
  })
})

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# PARSE COMMANDLINE ARGUMENTS

option_list <- list(
  make_option("--r1",
              metavar="R1_FASTQ",
              help = "R1 FASTQ file [REQUIRED]",
              type = "character"),
  make_option("--r2",
              metavar="R2_FASTQ",
              help = "R2 FASTQ file [REQUIRED]",
              type = "character"),
  make_option("--output",
              metavar="OUTPUT_FILENAME",
              help = "Output filename [REQUIRED]",
              type = "character"),
  make_option("--sample_name",
              help = "Sample name",
              type = "character"),
  make_option("--ctrl_depth",
              help = "Output control depth of coverage to output",
              action = "store_true",
              default = FALSE),
  make_option("--save_rdata_file",
              metavar="RDATA_FILENAME",
              help = "RData filename to save objects for debugging",
              type = "character"),
  make_option("--keep_tmp_files",
              help = "Keep tmp files for debugging",
              action = "store_true",
              default = FALSE),
  make_option("--threads",
              metavar="N",
              help = "Number of threads to use [default: 2]",
              type = "integer",
              default = 2),
  make_option("--version",
              help = "Print version and exit",
              action = "store_true",
              default = FALSE)
)

parser <- OptionParser(
  usage = paste("%prog [OPTIONS] --r1 [R1_FASTQ] --r2 [R2_FASTQ] --output [OUTPUT_FILENAME]\n",
                "Predict sex of Aedes albopictus from ddRAD sequencing ",
                "(nlaIII and mluCI restriction enzymes).",
                sep=""),
  option_list = option_list
)

# Parse command line arguments
INVALID_MESSAGE <- "Invalid command line arguments. Use --help for help."
tryCatch({
  suppressWarnings(
    arguments <- parse_args(parser, positional_arguments=0)
  )},
  error = function(e) {
    message(INVALID_MESSAGE)
    quit(status=2)
  }
)

opts <- arguments$options

# Print version
if (opts$version) {
  cat(basename(commandArgs()[4]), as.character(VERSION), "\n")
  quit(save = "no")
}

# Require r1 r2 output filename to be in opts
MISSING_OPTS_MESSAGE <- paste("Invalid command. --r1, --r2, and --output are",
                              "required. Use --help for help.")
if (! all(c("r1", "r2", "output") %in% names(opts))) {
  message(MISSING_OPTS_MESSAGE)
  quit(status=2)
}

# Check files exist
program_dir <- dirname(strsplit(commandArgs()[4], "=")[[1]][2])
model_rdata <- paste0(program_dir, "/model/model.RData")
feature_bed_path <- paste0(program_dir, "/index/depth_features.bed")
bt2_index <- paste0(program_dir, "/index/features")
MISSING_FILES_MESSAGE <- paste0("Files are missing from program directory at ",
                                program_dir, "/")
if (! (file.exists(model_rdata) &
       file.exists(feature_bed_path) &
       file.exists(paste0(program_dir, "/index/features.fa")) &
       file.exists(paste0(program_dir, "/index/features.rev.2.bt2")))) {
  message(MISSING_FILES_MESSAGE)
  quit(status=2)
}

# Get PID and tmp suffix
ps <- ps_handle()
pid <- strsplit(capture.output(ps), "[=,]")[[1]][2]
tmp_suffix <- paste0(".", pid, ".tmp")

message("Running: ", paste(commandArgs(), collapse=" "))
load(model_rdata)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# BOWTIE2 ALIGNMENT

output_filename <- opts$output
bam_output <- paste0(output_filename, tmp_suffix, ".bam")
bt2_threads <- max(opts$threads - 1, 1)
bt2_cmd <- sprintf("bowtie2 -p %d -x %s -1 %s -2 %s 2> /dev/null | samtools view -bS -f 3 - | samtools sort -o %s -",
                   bt2_threads, bt2_index, opts$r1, opts$r2, bam_output)

message("Running command:\n  ", bt2_cmd)
system(bt2_cmd)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# BEDTOOLS COVERAGE

bedtools_output <- paste0(output_filename, tmp_suffix)
bedtools_cmd <- sprintf("bedtools coverage -a %s -b %s > %s",
                        feature_bed_path, bam_output, bedtools_output)

message("Running command:\n  ", bedtools_cmd)
system(bedtools_cmd)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# LOAD AND NORMALISE COUNTS

# Load feature counts file
coverage_colnames <- c("scaffold_start_end", "start", "end", "coverage", 
                       "covered_bases", "total_bases", "percent_coverage")
feature_counts <- read.table(bedtools_output, header=FALSE,
                             stringsAsFactors=FALSE,
                             col.names=coverage_colnames)

# Get region names
tmp <- do.call(rbind, strsplit(feature_counts$scaffold_start_end, split="[:-]"))
feature_counts$region_names <- paste(tmp[,1], tmp[,2], tmp[,3], sep="_")
stopifnot(feature_counts$region_names == rad_df$region_name)
ctrl_regions <- rad_df[rad_df$type == "control", "region_name"]
feature_regions <- rad_df[rad_df$type == "feature", "region_name"]

# Estimate fragment size
rad_df$counts <- feature_counts[,"coverage"]
rad_df$log_counts <- log2(feature_counts[,"coverage"] + 1)
tmp <- rad_df[rad_df$region_name %in% ctrl_regions,]
fit <- lowess(x=tmp$fragment_length, y=tmp$log_counts, f=0.1)
fragment_size_estimate <- fit$x[which.max(fit$y)]
message("Fragment size estimate: ", fragment_size_estimate)

# Print warning if fragment size is less than 270 bp
if (fragment_size_estimate < 270) {
  message(paste0("Warning: Fragment size estimate is low. Prediction best ",
                 "works for fragment size ~300 bp."))
}

# Calculate normalisation factor based on fragment size estimate and the
# geometric mean within a window around the fragment size estimate
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}
window_size <- 25
ctrl_counts <- tmp[(tmp$fragment_length >= fragment_size_estimate - window_size) &
                     (tmp$fragment_length <= fragment_size_estimate + window_size),
                   "counts"]
ctrl_depth <- gm_mean(ctrl_counts)
message("Sample control depth: ", round(ctrl_depth, 2))

# Quit if ctrl depth is less than 1 and print warning if ctrl depth is 
# less than 10
rm_cmd <- paste("rm", bam_output, bedtools_output)
if (ctrl_depth < 1) {
  if (! opts$keep_tmp_files) {
    system(rm_cmd)
  }
  stop("Depth of sample is too low to predict sex classification. Exiting.")
} else if (ctrl_depth < 10) {
  message("Warning: Depth of sample is low. Classification may in inaccurate.")
}

# Normalise counts
feature_region_counts <- 
  feature_counts[feature_counts$region_names %in% feature_regions,]
norm_counts <- feature_region_counts[,"coverage"] / ctrl_depth * 1e3
log_norm_counts <- log2(norm_counts + 0.5)

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# SVM PREDICTION

# Wrangle data
x <- matrix(log_norm_counts, nrow=1)
colnames(x) <- feature_region_counts$region_names

# Predict sex
svm_pred <- predict(mod, x, probability=TRUE)
pred <- as.character(svm_pred)
pred <- ifelse(pred == "F", "Female", "Male")
prob <- attr(svm_pred, "probabilities")
message(sprintf("Sample was classified as %s with %0.1f%% probability.", 
                pred, max(prob) * 100))

# Output classification
if (! is.null(opts$sample_name)) {
  output <- data.frame(sample=opts$sample_name, class=pred, 
                       probability=round(max(prob), 4),
                       prob_female=round(prob[,"F"], 4), 
                       prob_male=round(prob[,"M"], 4))
} else {
  output <- data.frame(class=pred, 
                       probability=round(max(prob), 4),
                       prob_female=round(prob[,"F"], 4), 
                       prob_male=round(prob[,"M"], 4))
}
if (opts$ctrl_depth) {
  output$fragment_size_estimate <- fragment_size_estimate
  output$ctrl_depth <- round(ctrl_depth, 2)
}
write.table(output, file=output_filename, row.names=FALSE, sep="\t", 
            quote=FALSE)

# Save RData
if (! is.null(opts$save_rdata_file)) {
  save(feature_counts, ctrl_depth, log_norm_counts, output,
       file=opts$save_rdata_file)
}

#=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
# REMOVE TMP FILES

if (! opts$keep_tmp_files) {
  system(rm_cmd)
}
