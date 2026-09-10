# Extended Data Figure 2 | Normalization scale factors and stage-specific detection

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "readxl", "patchwork", "ggpubr", "openxlsx", "showtext"), "SomaDataIO")
source("utils/figure_style.R")
source("utils/feature_mapping.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
showtext_auto(); showtext_opts(dpi = 600)
label_style <- list(size = 12, face = "bold")
paths <- c(metadata = "data/study_metadata.xlsx", feature_metadata = "data/feature_metadata.tsv.gz", profiles = "data/protein_profiles_long.tsv.gz",
           detection = "results/detection_status.tsv.gz")
adat_paths <- c(
    SOM_P1_B1 = upstream_path("raw", "SOM_ACM_241017", "normalized-adat",
                              "FDDX.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"),
    SOM_P1_B2 = upstream_path("raw", "SOM_ACM_251113", "LYQ.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"),
    SOM_P2_B1 = upstream_path("raw", "SOM_ISG_250919", "Fudan_18Sep2025_Fudan_19Sep2025_Step7_FinalNormStep_MedNormExt-R.adat"),
    SOM_P2_B2 = upstream_path("raw", "SOM_ISH_260331", "NovaSeqX-Fudan260331_dragen-protein-quant_Step3_SampleNorm.adat")
)
upstream_metadata <- upstream_path("metadata", "metadata_v3.0-ms_master.xlsx")
missing_inputs <- c(paths[!file.exists(paths)], adat_paths[!file.exists(adat_paths)], upstream_metadata[!file.exists(upstream_metadata)])
if (length(missing_inputs)) stop("Extended Data Figure 2 is missing input files:\n", paste(missing_inputs, collapse = "\n"))

# 1. Inputs and SomaScan scale factors ----
meta_batch <- read_xlsx(paths["metadata"], sheet = "batch")
meta_sample <- read_xlsx(paths["metadata"], sheet = "sample")
raw_meta_sample <- read_xlsx(upstream_metadata, sheet = "sample")
feature_metadata <- fread(paths["feature_metadata"]) %>% as_tibble()
long_df <- fread(paths["profiles"]) %>% as_tibble() %>% filter_batch_analysis_features(feature_metadata, strict_platforms = character())

extract_scale_factors <- function(file_path, batch_name) {
    adat <- suppressWarnings(SomaDataIO::read_adat(file_path))
    meta_sub <- raw_meta_sample %>% filter(Batch == batch_name)
    match_ext <- if ("ExtIdentifier" %in% colnames(adat)) sum(adat$ExtIdentifier %in% meta_sub$RawName) else 0
    match_smp <- if ("SampleID" %in% colnames(adat)) sum(adat$SampleID %in% meta_sub$RawName) else 0
    id_col <- if (match_ext > match_smp) "ExtIdentifier" else "SampleID"
    sf_hyb <- grep("^HybNorm_1_ScaleFactor$|^HybControlNormScale$", colnames(adat), value = TRUE, ignore.case = TRUE)[1]
    sf_term_cols <- grep("^NormScale_[0-9]|^MedNormExt_.*_ScaleFactor", colnames(adat), value = TRUE)
    if (!length(sf_term_cols)) sf_term_cols <- grep("^MedNormInt_.*_ScaleFactor", colnames(adat), value = TRUE)
    if (is.na(sf_hyb) || !length(sf_term_cols)) stop(batch_name, " ADAT lacks required normalization scale-factor columns.")

    adat %>%
        select(id_val = all_of(id_col), HybNorm_SF = all_of(sf_hyb), all_of(sf_term_cols)) %>%
        inner_join(meta_sub, by = c("id_val" = "RawName")) %>%
        filter(Include == TRUE | Sample %in% c("BLK", "CAL", "QC")) %>%
        pivot_longer(all_of(sf_term_cols), names_to = "Terminal_SF_Type", values_to = "Terminal_SF") %>%
        mutate(
            Batch = batch_name,
            Dilution_Bin = case_when(
                str_detect(Terminal_SF_Type, "NormScale_20|MedNorm.*_0\\.2_") ~ "High conc",
                str_detect(Terminal_SF_Type, "NormScale_0\\.5|MedNorm.*_0\\.005_") ~ "Med conc",
                str_detect(Terminal_SF_Type, "NormScale_0\\.005|MedNorm.*_5e-05_") ~ "Low conc",
                TRUE ~ "Unknown"
            )
        )
}

som_scale_factor <- bind_rows(lapply(names(adat_paths), function(batch_name) extract_scale_factors(adat_paths[[batch_name]], batch_name))) %>%
    mutate(Sample = factor(Sample, levels = c("M", "Y", "P", "X", "F", "N", "CAL", "QC", "BLK")))
som_sf_hyb <- som_scale_factor %>% distinct(Batch, Sample, ColName, HybNorm_SF)

# 2. Panel a: readout-normalization scale factors ----
p_left <- ggplot(som_sf_hyb %>% filter(Batch %in% c("SOM_P1_B1", "SOM_P1_B2")), aes(Sample, HybNorm_SF, fill = Sample)) +
    geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.6, color = "black", linewidth = 0.2) +
    geom_jitter(alpha = 1, size = 1, width = 0.15, shape = 21, color = "white", stroke = 0.3) +
    facet_wrap(~Batch) +
    scale_fill_manual(values = sample_color) +
    labs(y = "Readout scale factor", x = NULL) +
    plasmix_theme +
    theme(legend.position = "none", panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_right <- ggplot(som_sf_hyb %>% filter(Batch %in% c("SOM_P2_B1", "SOM_P2_B2")), aes(Sample, HybNorm_SF, fill = Sample)) +
    geom_boxplot(outlier.shape = NA, alpha = 1, width = 0.6, color = "black", linewidth = 0.2) +
    geom_jitter(alpha = 1, size = 1, width = 0.15, shape = 21, color = "white", stroke = 0.3) +
    facet_wrap(~Batch) +
    scale_fill_manual(values = sample_color) +
    labs(y = NULL, x = NULL) +
    plasmix_theme +
    theme(legend.position = "none", panel.grid.major = element_blank(), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))

p_scale_factors <- p_left + p_right + plot_layout(widths = c(1, 1)) & theme(plot.margin = margin(2.5, 3, 2.5, 3))

# 3. Stage-specific feature detection ----
affinity_batches <- c("OLK_P2_B1", "OLK_P2_B2", "SOM_P1_B1", "SOM_P1_B2", "SOM_P2_B1", "SOM_P2_B2")
meta_batch_af <- meta_batch %>% filter(Batch %in% affinity_batches)
meta_sample_af <- meta_sample %>% filter(Batch %in% affinity_batches, Sample %in% c("M", "Y", "P", "X", "F", "N", "BLK"))

calculate_strict_detection_status <- function(df_data, threshold_df, sample_metadata, stage_name) {
    group_sizes <- sample_metadata %>% filter(Sample != "BLK") %>% count(Batch, Sample, name = "ExpectedN")
    feature_grid <- df_data %>% distinct(Batch, Platform, UniqueID)

    group_status <- df_data %>%
        filter(Sample != "BLK") %>%
        inner_join(threshold_df, by = c("Platform", "Batch", "UniqueID")) %>%
        mutate(Detected = is.finite(Log2_Value) & Log2_Value > LoD_Threshold) %>%
        group_by(Batch, Platform, UniqueID, Sample) %>%
        summarize(DetectedN = sum(Detected), .groups = "drop") %>%
        right_join(feature_grid %>% crossing(Sample = c("M", "Y", "P", "X", "F", "N")),
                   by = c("Batch", "Platform", "UniqueID", "Sample")) %>%
        left_join(group_sizes, by = c("Batch", "Sample")) %>%
        mutate(DetectedN = replace_na(DetectedN, 0L), GroupPass = is.finite(ExpectedN) & DetectedN / ExpectedN > 0.5) %>%
        select(Batch, Platform, UniqueID, Sample, GroupPass) %>%
        pivot_wider(names_from = Sample, values_from = GroupPass, values_fill = FALSE)

    for (sample_name in c("M", "Y", "P", "X", "F", "N")) {
        if (!sample_name %in% names(group_status)) group_status[[sample_name]] <- FALSE
    }
    group_status %>%
        mutate(PassedGroup = rowSums(across(all_of(c("M", "Y", "P", "X", "F", "N")))) ,
               IsAboveLoD = PassedGroup > 0, JumpNote = NA_character_, Stage = stage_name) %>%
        select(Batch, Platform, UniqueID, M, Y, P, X, F, N, PassedGroup, JumpNote, IsAboveLoD, Stage)
}

df_raw <- long_df %>%
    filter(Batch %in% affinity_batches, Platform %in% c("OLK", "SOM"), ProcessLevel == "Raw") %>%
    mutate(Linear_Value = 2^Value, Log2_Value = Value)
df_readout <- long_df %>%
    filter(Batch %in% affinity_batches, Platform %in% c("OLK", "SOM"), ProcessLevel == "HybNorm") %>%
    mutate(Linear_Value = 2^Value, Log2_Value = Value)

raw_thresholds <- df_raw %>% filter(Sample == "BLK") %>%
    group_by(Platform, Batch, UniqueID) %>%
    summarize(BLK_Median = median(Linear_Value, na.rm = TRUE), BLK_MAD = mad(Linear_Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(LoD_Threshold = log2(pmax(BLK_Median + 3 * BLK_MAD, 1e-9)))
readout_thresholds <- df_readout %>% filter(Sample == "BLK") %>%
    group_by(Platform, Batch, UniqueID) %>%
    summarize(BLK_Median = median(Linear_Value, na.rm = TRUE), BLK_MAD = mad(Linear_Value, na.rm = TRUE), .groups = "drop") %>%
    mutate(LoD_Threshold = log2(pmax(BLK_Median + 3 * BLK_MAD, 1e-9)))

status_raw <- calculate_strict_detection_status(df_raw, raw_thresholds, meta_sample_af, "Raw readout")
status_readout <- calculate_strict_detection_status(df_readout, readout_thresholds, meta_sample_af, "Readout norm")
status_plate <- fread(paths["detection"]) %>% as_tibble() %>%
    semi_join(distinct(long_df, Platform, Batch, UniqueID), by = c("Platform", "Batch", "UniqueID")) %>%
    filter(Batch %in% affinity_batches) %>%
    mutate(Stage = "Plate norm") %>%
    select(any_of(c("Batch", "Platform", "UniqueID", "M", "Y", "P", "X", "F", "N", "PassedGroup", "JumpNote", "IsAboveLoD", "Stage")))
if (!"IsAboveLoD" %in% names(status_plate)) stop("detection_status.tsv.gz lacks IsAboveLoD.")

feature_detection_status <- bind_rows(status_raw, status_readout, status_plate)
detection_summary <- feature_detection_status %>%
    group_by(Platform, Batch, Stage) %>%
    summarize(Total = n_distinct(UniqueID), Detected = sum(IsAboveLoD, na.rm = TRUE), `Detection rate (%)` = 100 * Detected / Total,
              .groups = "drop") %>%
    mutate(Stage = factor(Stage, levels = c("Raw readout", "Readout norm", "Plate norm")))

stage_colors <- c("Raw readout" = "#c6d5e9", "Readout norm" = "#358DB9FF", "Plate norm" = "#26185F")
p_detection <- ggplot(detection_summary, aes(Batch, `Detection rate (%)`, fill = Stage)) +
    geom_col(position = position_dodge(width = 0.6), width = 0.6) +
    scale_fill_manual(values = stage_colors) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0)) +
    labs(y = "Overall detection rate (%)", x = NULL, fill = "Stage") +
    plasmix_theme +
    theme(axis.text.x = element_text(angle = 20, hjust = 1), panel.grid.major.x = element_blank(), legend.position = "top",
          legend.margin = margin(0, 5, 0, -10), legend.title = element_text(size = 7.5, face = "bold", vjust = 0.5),
          legend.text = element_text(size = 7.5, vjust = 0.5, margin = margin(0, 0, 0, 3)))

# 4. Assemble and export ----
figure_ed2 <- ggarrange(p_scale_factors, p_detection, nrow = 1, widths = c(1.5, 1), labels = c("a", "b"),
                                label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2, font.label = label_style)
ggsave("figures/ed2_normalization_detection.pdf", figure_ed2, width = 10, height = 2.8)
ggsave("figures/ed2_normalization_detection.png", figure_ed2, width = 10, height = 2.8, dpi = 600, bg = "white")

source_data <- list(
    ED2a_scale_factors = som_sf_hyb %>% arrange(Batch, Sample, ColName),
    ED2b_detection_summary = detection_summary %>% mutate(Stage = as.character(Stage)) %>% arrange(Platform, Batch, Stage),
    ED2b_feature_status = feature_detection_status %>% arrange(Platform, Batch, Stage, UniqueID),
    ED2b_raw_thresholds = raw_thresholds %>% arrange(Platform, Batch, UniqueID),
    ED2b_readout_thresholds = readout_thresholds %>% arrange(Platform, Batch, UniqueID)
)
write.xlsx(source_data, "tables/SourceData_EDFigure2.xlsx", overwrite = TRUE, keepNA = TRUE, na.string = "NA")
message("Extended Data Figure 2 and its source data were exported.")
