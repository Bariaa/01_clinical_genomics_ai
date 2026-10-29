# 14_generate_epigenomic_figures.R
# Publication-quality figures for the methylation layer.
#
# Outputs:
#   outputs/figures/methylation_pca.png
#   outputs/figures/methylation_heatmap.png
#   outputs/figures/methylation_beta_distribution.png
#   outputs/figures/methylation_by_clinical_group.png

library(ggplot2)
library(dplyr)
library(yaml)
library(pheatmap)

config <- yaml::read_yaml("config/config.yaml")
proc_dir <- config$paths$data_processed

integrated <- readRDS(file.path(proc_dir, "multiomics_integrated.rds"))
meth_cols <- grep("^cg", names(integrated), value = TRUE)
beta_mat <- as.matrix(integrated[, meth_cols])
rownames(beta_mat) <- integrated$patient_id

message("Figure input: ", paste(dim(beta_mat), collapse = " x "))
dir.create("outputs/figures", recursive = TRUE, showWarnings = FALSE)

W <- config$plots$figure_width
H <- config$plots$figure_height
DPI <- config$plots$dpi

# --- 1. PCA ---
message("Computing PCA...")
pca <- prcomp(beta_mat, center = TRUE, scale. = FALSE)
var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  patient_id = rownames(beta_mat),
  vital_status = integrated$vital_status,
  TP53_mutated = factor(integrated$TP53_mutated, levels = c(0, 1),
                        labels = c("TP53 wild-type", "TP53 mutated"))
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, colour = TP53_mutated)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_colour_manual(values = c("#1C7293", "#F96167")) +
  labs(
    title = "PCA of DNA Methylation Profiles — TCGA-BRCA",
    subtitle = paste0(nrow(beta_mat), " patients, ", length(meth_cols), " variable CpG probes"),
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)"),
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

ggsave("outputs/figures/methylation_pca.png", p_pca, width = W, height = H, dpi = DPI)
message("✅ methylation_pca.png")

# --- 2. Heatmap: top 500 most variable probes ---
message("Building heatmap...")
probe_sd <- apply(beta_mat, 2, sd)
top500 <- names(sort(probe_sd, decreasing = TRUE))[1:500]
hm_mat <- t(beta_mat[, top500])

annot_col <- data.frame(
  `Vital status` = integrated$vital_status,
  TP53 = factor(integrated$TP53_mutated, levels = c(0, 1),
                labels = c("wild-type", "mutated")),
  check.names = FALSE,
  row.names = integrated$patient_id
)

annot_colours <- list(
  TP53 = c(`wild-type` = "#CADCFC", mutated = "#F96167"),
  `Vital status` = c(alive = "#1C7293", dead = "#21295C")
)

png("outputs/figures/methylation_heatmap.png",
    width = W, height = H, units = "in", res = DPI)
pheatmap(
  hm_mat,
  show_rownames = FALSE,
  show_colnames = FALSE,
  annotation_col = annot_col,
  annotation_colors = annot_colours,
  color = colorRampPalette(c("#065A82", "#F4F7FA", "#F96167"))(100),
  main = "DNA Methylation — Top 500 Most Variable CpG Probes",
  fontsize = 10,
  treeheight_row = 25,
  treeheight_col = 25
)
dev.off()
message("✅ methylation_heatmap.png")

# --- 3. Beta value distribution ---
message("Building beta distribution plot...")
set.seed(config$project$seed)
sample_probes <- sample(meth_cols, min(2000, length(meth_cols)))
beta_long <- data.frame(beta = as.vector(beta_mat[, sample_probes]))

p_dist <- ggplot(beta_long, aes(x = beta)) +
  geom_histogram(bins = 60, fill = "#065A82", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = c(0.3, 0.7), linetype = "dashed", colour = "#F96167") +
  labs(
    title = "Distribution of DNA Methylation Beta Values — TCGA-BRCA",
    subtitle = "Dashed lines mark conventional hypo- (<0.3) and hypermethylation (>0.7) thresholds",
    x = "Beta value (proportion methylated)",
    y = "Frequency"
  ) +
  theme_minimal(base_size = 12)

ggsave("outputs/figures/methylation_beta_distribution.png", p_dist,
       width = W, height = H, dpi = DPI)
message("✅ methylation_beta_distribution.png")

# --- 4. Mean methylation by clinical group ---
message("Building clinical group comparison...")
group_df <- data.frame(
  mean_methylation = rowMeans(beta_mat),
  vital_status = integrated$vital_status,
  TP53 = factor(integrated$TP53_mutated, levels = c(0, 1),
                labels = c("TP53 wild-type", "TP53 mutated"))
) %>% filter(!is.na(vital_status))

p_group <- ggplot(group_df, aes(x = TP53, y = mean_methylation, fill = TP53)) +
  geom_boxplot(outlier.alpha = 0.4, width = 0.55) +
  scale_fill_manual(values = c("#1C7293", "#F96167")) +
  labs(
    title = "Global Mean Methylation by TP53 Mutation Status",
    subtitle = paste0(nrow(group_df), " patients with all three omics layers"),
    x = NULL,
    y = "Mean beta value across variable probes"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

ggsave("outputs/figures/methylation_by_clinical_group.png", p_group,
       width = 7, height = H, dpi = DPI)
message("✅ methylation_by_clinical_group.png")

message("\nAll four figures saved to outputs/figures/")

# --- 5. PCA coloured by pathologic stage ---
message("Building stage-coloured PCA...")

clinical_full <- read.csv("data/processed/clinical_cleaned.csv", stringsAsFactors = FALSE)
clinical_full$patient_id <- toupper(trimws(clinical_full$patient_id))

stage_lookup <- clinical_full[, c("patient_id", "ajcc_pathologic_stage")]
stage_lookup$stage_group <- dplyr::case_when(
  grepl("^Stage IV",  stage_lookup$ajcc_pathologic_stage) ~ "Stage IV",
  grepl("^Stage III", stage_lookup$ajcc_pathologic_stage) ~ "Stage III",
  grepl("^Stage II",  stage_lookup$ajcc_pathologic_stage) ~ "Stage II",
  grepl("^Stage I",   stage_lookup$ajcc_pathologic_stage) ~ "Stage I",
  TRUE ~ NA_character_
)

pca_stage <- pca_df %>%
  dplyr::left_join(stage_lookup[, c("patient_id", "stage_group")], by = "patient_id") %>%
  dplyr::filter(!is.na(stage_group))

message("Patients with usable stage: ", nrow(pca_stage), " / ", nrow(pca_df))

p_pca_stage <- ggplot(pca_stage, aes(x = PC1, y = PC2, colour = stage_group)) +
  geom_point(alpha = 0.75, size = 2) +
  scale_colour_manual(values = c("Stage I"   = "#CADCFC",
                                 "Stage II"  = "#1C7293",
                                 "Stage III" = "#065A82",
                                 "Stage IV"  = "#F96167")) +
  labs(
    title = "PCA of DNA Methylation Profiles by Pathologic Stage",
    subtitle = paste0(nrow(pca_stage), " patients with AJCC stage recorded"),
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)"),
    colour = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

ggsave("outputs/figures/methylation_pca_by_stage.png", p_pca_stage,
       width = W, height = H, dpi = DPI)
message("✅ methylation_pca_by_stage.png")
