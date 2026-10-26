# 10_prepare_methylation_data.R
# Download / prepare raw DNA methylation data (Illumina 450K) from GDC.
# 10_prepare_methylation_data.R
# Prepare DNA methylation data (Illumina 450K) for the PRIMARY dataset (TCGA-BRCA).
#
# DESIGN NOTE: GDC's beta-value files are read directly rather than via
# TCGAbiolinks::GDCprepare(). GDCprepare() builds a full SummarizedExperiment
# (~485,000 probes x 895 samples, several GB in memory) and requires the
# SeSAMe annotation stack. Reading the files directly is faster, has a far
# smaller memory and disk footprint, and lets us retain only the most
# variable probes — which is required for modelling regardless, given
# 485,000 probes against ~1,000 patients.
#
# Files are tab-delimited, one per aliquot, with no header:
#   column 1 = CpG probe ID, column 2 = beta value
#
# Output: data/raw/TCGA-BRCA/TCGA-BRCA_methylation_betas.rds
#         (matrix: top N most variable probes x samples)

library(data.table)
library(yaml)

config <- yaml::read_yaml("config/config.yaml")
primary_id <- config$datasets$primary$name
raw_dir <- config$paths$data_raw

# Number of most-variable probes to retain (configurable)
N_PROBES <- 20000

meth_dir <- file.path(raw_dir, primary_id, "DNA_Methylation")
beta_files <- list.files(meth_dir, pattern = "level3betas\\.txt$",
                         recursive = TRUE, full.names = TRUE)

message("Found ", length(beta_files), " beta-value files.")
if (length(beta_files) == 0) stop("No methylation files found under ", meth_dir)

# --- Map file -> TCGA barcode via the GDC sample sheet in the query results ---
# The UUID is the file's parent folder name; we resolve barcodes afterwards
# using the GDC query results (script 12 handles ID harmonisation in full).
file_uuids <- basename(dirname(beta_files))

# --- Pass 1: read first file to establish probe order ---
message("Reading probe index from first file...")
first <- data.table::fread(beta_files[1], header = FALSE,
                           col.names = c("probe", "beta"),
                           colClasses = c("character", "numeric"))
probe_ids <- first$probe
n_probes <- length(probe_ids)
message("Probes per sample: ", n_probes)

# --- Pass 2: read all files into a matrix ---
message("Reading all samples (this takes a few minutes)...")
beta_mat <- matrix(NA_real_, nrow = n_probes, ncol = length(beta_files),
                   dimnames = list(probe_ids, file_uuids))

pb <- txtProgressBar(min = 0, max = length(beta_files), style = 3)
for (i in seq_along(beta_files)) {
  d <- data.table::fread(beta_files[i], header = FALSE,
                         col.names = c("probe", "beta"),
                         colClasses = c("character", "numeric"))
  # Guard against probe-order differences between files
  if (!identical(d$probe, probe_ids)) {
    d <- d[match(probe_ids, d$probe), ]
  }
  beta_mat[, i] <- d$beta
  setTxtProgressBar(pb, i)
}
close(pb)

message("Full matrix dimensions: ", paste(dim(beta_mat), collapse = " x "))

# --- Reduce to the most variable probes ---
message("Selecting top ", N_PROBES, " most variable probes...")
probe_sd <- apply(beta_mat, 1, sd, na.rm = TRUE)
keep <- order(probe_sd, decreasing = TRUE, na.last = NA)[seq_len(min(N_PROBES, sum(!is.na(probe_sd))))]
beta_top <- beta_mat[keep, , drop = FALSE]

message("Reduced matrix dimensions: ", paste(dim(beta_top), collapse = " x "))

# --- Save ---
out_path <- file.path(raw_dir, primary_id, paste0(primary_id, "_methylation_betas.rds"))
saveRDS(beta_top, out_path)

# Also save the full probe SD vector — script 11 uses it for annotation summaries
saveRDS(probe_sd, file.path(raw_dir, primary_id, paste0(primary_id, "_probe_sd.rds")))

message("Saved: ", out_path)
message("Size on disk: ", round(file.size(out_path)/1024^2, 1), " MB")

