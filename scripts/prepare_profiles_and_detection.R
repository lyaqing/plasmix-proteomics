# Prepare public protein profiles and detection status

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "readxl", "SomaDataIO", "openxlsx"))
set.seed(2026)
options(stringsAsFactors = FALSE)
# 1. Metadata and reference mapping ----
raw_ref <- fread(upstream_path("references/uniprotkb_AND_model_organism_9606_2026_04_01.tsv.gz"), data.table = FALSE)
uniprot_ref <- raw_ref %>%
  filter(Reviewed == "reviewed") %>%
  select(
      Protein.ID = Entry, Protein.Name = `Entry Name`, Protein.Fullname = `Protein names`,
      Gene.Symbol = `Gene Names (primary)`, Gene.Synonyms = `Gene Names (synonym)`
  ) %>%
  mutate(Gene.Symbol = as.character(Gene.Symbol))
meta_batch <- read_xlsx(upstream_path("metadata/metadata_v3.0-ms_master.xlsx"), sheet = "batch")
meta_sample <- read_xlsx(upstream_path("metadata/metadata_v3.0-ms_master.xlsx"), sheet = "sample")
meta_variance <- read_xlsx(upstream_path("metadata/metadata_v3.0-ms_master.xlsx"), sheet = "variance")

long_data_list <- list()

# 2. Shared profile-import helpers ----
clean_and_sort_str <- function(str_vec) {
  sapply(str_vec, function(x) {
    if (is.na(x) || x == "") return(NA)
    parts <- unlist(strsplit(as.character(x), "[:;_|,]"))
    parts <- trimws(parts)
    parts <- parts[parts != ""]
    if (length(parts) == 0) return(NA)
    paste(sort(unique(parts)), collapse = "|")
  })
}

standardize_long <- function(df, batch_id, proc_level) {
  meta_sub <- meta_sample %>% filter(Batch == batch_id)
  if(!"LOD" %in% colnames(df)) {
    df$LOD <- NA_real_
  }
  df %>%
    inner_join(meta_sub %>% select(RawName, ColName, Platform, Sample, Include), by = "RawName") %>%
    mutate(Batch = batch_id, ProcessLevel = proc_level) %>%
    select(Batch, Platform, Sample, ColName, RawName, AssayID, TargetName, UniProtID, Value, LOD, ProcessLevel, Include)
}

# 3. SomaScan and Illumina Protein Prep ----
read_soma_steps <- function(base_dir, file_prefix, step_map, batch_id, anno_file, meta_sample) {
  if (batch_id %in% c("SOM_P2_B1", "SOM_P2_B2")) {
    anno_df <- read_excel(anno_file, range = "B7:I10464")
    tgt_col <- grep("Target", colnames(anno_df), value = TRUE)[1]
    anno_df <- anno_df %>%
      rename(SeqId = SeqID, UniProtID = `UniProt ID`, TargetName = !!sym(tgt_col))
  } else {
    anno_df <- read_annotations(anno_file)
    tgt_col <- grep("Target", colnames(anno_df), value = TRUE)[1]
    anno_df <- anno_df %>% rename(TargetName = !!sym(tgt_col), UniProtID = UniProt)
  }
  anno <- anno_df %>%
    mutate(AssayID = paste0("seq.", gsub("-", ".", SeqId))) %>%
    filter(Organism == "Human") %>%
    mutate(UniProtID = clean_and_sort_str(UniProtID)) %>%
    select(AssayID, TargetName, UniProtID)
  meta_sub <- meta_sample %>% filter(Batch == batch_id) %>% select(Batch, RawName, Sample)
  soma_list <- list()
  for (unified_name in names(step_map)) {
    suffix <- step_map[[unified_name]]
    fname <- paste0(file_prefix, suffix)
    full_path <- file.path(base_dir, fname)
    if (file.exists(full_path)) {
      message(sprintf("Processing [%s]: %s -> %s", batch_id, unified_name, fname))
      adat <- suppressMessages(read_adat(full_path))
      buffer_dat <- adat %>% filter(SampleType %in% c("Buffer", "Blank"))
      if(nrow(buffer_dat) > 0) {
        lod_vals <- buffer_dat %>%
          select(starts_with("seq")) %>%
          summarize(across(everything(), ~ median(.x, na.rm = TRUE) + 3 * mad(.x, na.rm = TRUE)), .groups = "drop") %>%
          pivot_longer(everything(), names_to = "AssayID", values_to = "LOD_Linear") %>%
          mutate(LOD = log2(LOD_Linear)) %>%
          select(-LOD_Linear)
      } else {
        message(sprintf("Warning: No Buffer samples found in %s. LOD set to NA.", fname))
        lod_vals <- tibble(AssayID = colnames(adat)[grep("^seq", colnames(adat))], LOD = NA_real_)
      }
      id_col <- ifelse("ExtIdentifier" %in% colnames(adat), "ExtIdentifier", "SampleID")
      df_long_raw <- adat %>%
        select(all_of(id_col), starts_with("seq")) %>%
        pivot_longer(cols = starts_with("seq"), names_to = "AssayID", values_to = "Intensity") %>%
        mutate(Value = log2(Intensity)) %>%
        rename(RawName = all_of(id_col)) %>%
        inner_join(anno, by = "AssayID") %>%
        left_join(lod_vals, by = "AssayID")
      df_long_original <- df_long_raw
      soma_list[[unified_name]] <- standardize_long(df_long_original, batch_id, unified_name)
      }
  }
  result <- bind_rows(soma_list)
  return(result)
}

step_map_11k <- c(
    "Raw" = ".adat", "HybNorm" = ".hybNorm.adat", "MedNormInt" = ".hybNorm.medNormInt.adat",
    "PlateScale" = ".hybNorm.medNormInt.plateScale.adat",
    "Calibrate" = ".hybNorm.medNormInt.plateScale.calibrate.adat",
    "ANML-SMP" = ".hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.adat"
)
long_data_list[["SOM_P1_B1"]] <- read_soma_steps(
    base_dir = upstream_path("raw/SOM_ACM_241017/normalized-adat"), file_prefix = "FDDX",
    step_map = step_map_11k, batch_id = "SOM_P1_B1",
    anno_file = upstream_path("references/SomaScan_11K_Annotated_Content.xlsx"), meta_sample = meta_sample
)

long_data_list[["SOM_P1_B2"]] <- read_soma_steps(
    base_dir = upstream_path("raw/SOM_ACM_251113"), file_prefix = "LYQ", step_map = step_map_11k,
    batch_id = "SOM_P1_B2", anno_file = upstream_path("references/SomaScan_11K_Annotated_Content.xlsx"),
    meta_sample = meta_sample
)
step_map_ipp <- c(
    "Raw" = "Step0_Raw-R.adat", "HybNorm" = "Step1_HybNorm-R.adat", "MedNormInt" = "Step2_MedNormInt-R.adat",
    "PlateScale" = "Step5_CrossPlatformPlateScale-R.adat", "Calibrate" = "Step6_CrossPlatformCalibrate-R.adat",
    "MedNormExt" = "Step7_FinalNormStep_MedNormExt-R.adat"
)
long_data_list[["SOM_P2_B1"]] <- read_soma_steps(
    base_dir = upstream_path("raw/SOM_ISG_250919"), file_prefix = "Fudan_18Sep2025_Fudan_19Sep2025_",
    step_map = step_map_ipp, batch_id = "SOM_P2_B1",
    anno_file = upstream_path("references/Customer Facing - Protein Prep SOMAmer Annotations SEPT2025.xlsx"),
    meta_sample = meta_sample
)
step_map_ipp_v2 <- c(
    'Raw' = 'Step0_Raw.adat', 'HybNorm' = 'Step1_ReadoutNorm.adat', 'Calibrate' = 'Step2_PlateNorm.adat',
    'MedNormExt' = 'Step3_SampleNorm.adat'
)
long_data_list[["SOM_P2_B2"]] <- read_soma_steps(
    base_dir = upstream_path("raw/SOM_ISH_260331"), file_prefix = "NovaSeqX-Fudan260331_dragen-protein-quant_",
    step_map = step_map_ipp_v2, batch_id = "SOM_P2_B2",
    anno_file = upstream_path("references/Customer Facing - Protein Prep SOMAmer Annotations SEPT2025.xlsx"),
    meta_sample = meta_sample
)

# 4. Olink ----
process_old_olink <- function(filepath, bat_id) {
  raw <- read.csv(filepath)
  if(bat_id == "OLK_P1_B2") {
    raw[raw$Index == 73, ]$SampleID <- "SC_1"
    raw[raw$Index == 85, ]$SampleID <- "SC_2"
  }
  raw %>%
    mutate(
        Assay = case_when(Assay == "NTproBNP" & UniProt == "NTproBNP" ~ "NT-proBNP", TRUE ~ Assay),
        UniProt = case_when(Assay == "NTproBNP" & UniProt == "NTproBNP" ~ "NT-proBNP", TRUE ~ UniProt)
    ) %>%
    mutate(UniProtID = clean_and_sort_str(UniProt), TargetName = clean_and_sort_str(Assay)) %>%
    select(SampleID, OlinkID, TargetName, UniProtID, NPX, LOD) %>%
    rename(RawName = SampleID, AssayID = OlinkID, Value = NPX) %>%
    standardize_long(bat_id, "NPX")
}
long_data_list[["OLK_P1_B1"]] <- process_old_olink(upstream_path("raw/OLK_SQT_220303/Explore384_Cardiometabolic.csv"), "OLK_P1_B1")
long_data_list[["OLK_P1_B2"]] <- process_old_olink(upstream_path("raw/OLK_SQT_220830/Fudan_15Plasma_NPX.csv"), "OLK_P1_B2")

process_olink_ht <- function(filepath, bat_id) {
  raw_olk <- read_xlsx(filepath)
  df_ctrl <- raw_olk %>%
    filter(AssayType == "ext_ctrl") %>%
    select(SampleID, Block, Count) %>%
    rename(IC_Count = Count) %>%
    mutate(ScaleFactor = IC_Count / median(IC_Count, na.rm = TRUE))
  df_olk_base <- raw_olk %>%
    filter(AssayType == "assay") %>%
    mutate(UniProtID = clean_and_sort_str(UniProt), TargetName = clean_and_sort_str(Assay))
  df_count_raw <- df_olk_base %>%
    select(SampleID, OlinkID, TargetName, UniProtID, Count) %>%
    rename(RawName = SampleID, AssayID = OlinkID, Value = Count) %>%
    mutate(Value = log2(Value), LOD = NA_real_) %>%
    standardize_long(bat_id, "Raw")
  df_hybnorm <- df_olk_base %>%
    left_join(df_ctrl, by = c("SampleID", "Block")) %>%
    mutate(Value = log2(Count / ScaleFactor), LOD = NA_real_) %>%
    select(SampleID, OlinkID, TargetName, UniProtID, Value) %>%
    rename(RawName = SampleID, AssayID = OlinkID) %>%
    standardize_long(bat_id, "HybNorm")
  df_npx <- df_olk_base %>%
    select(SampleID, OlinkID, TargetName, UniProtID, NPX, LOD) %>%
    mutate(LOD = as.numeric(LOD)) %>%
    rename(RawName = SampleID, AssayID = OlinkID, Value = NPX) %>%
    standardize_long(bat_id, "NPX")
  bind_rows(df_count_raw, df_hybnorm, df_npx)
}

long_data_list[["OLK_P2_B1"]] <- process_olink_ht(upstream_path("raw/OLK_SNT_250116/JZ202412260940_OLINK_NPX_withLOD.xlsx"), "OLK_P2_B1")
long_data_list[["OLK_P2_B2"]] <- process_olink_ht(upstream_path("raw/OLK_SNT_260112/JZ202511210920_OLINK_NPX_withLOD.xlsx"), "OLK_P2_B2")

# 5. NULISA ----
calc_nulisa_lod <- function(df, nc_cols, data_type, ic_scale_df = NULL) {
  if(length(nc_cols) == 0) return(NULL)
  nc_mat <- df %>% select(targetName, all_of(nc_cols)) %>% column_to_rownames("targetName")
  if (data_type == "HybNorm") {
    if (is.null(ic_scale_df)) stop("Need IC scale dataframe for HybNorm LOD")
    sf_nc <- ic_scale_df %>% filter(RawName %in% nc_cols)
    sf_vec <- sf_nc$ScaleFactor[match(nc_cols, sf_nc$RawName)]
    nc_mat_adj <- nc_mat + 0.1
    nc_mat_norm <- sweep(nc_mat_adj, 2, sf_vec, "/")
    lod_vec_linear <- apply(nc_mat_norm, 1, function(x) mean(x, na.rm=TRUE) + 3*sd(x, na.rm=TRUE))
    lod_vec <- log2(lod_vec_linear)
  } else if (data_type == "Raw") {
    nc_mat_adj <- nc_mat + 0.1
    lod_vec_linear <- apply(nc_mat_adj, 1, function(x) mean(x, na.rm=TRUE) + 3*sd(x, na.rm=TRUE))
    lod_vec <- log2(lod_vec_linear)
  } else {
    lod_vec <- apply(nc_mat, 1, function(x) mean(x, na.rm=TRUE) + 3*sd(x, na.rm=TRUE))
  }
  tibble(TargetName = names(lod_vec), LOD = lod_vec)
}

process_nulisa_step <- function(df_wide, val_col, level_name, lod_df, anno_df, ic_df=NULL) {
  df_long <- df_wide %>%
    filter(targetName != "mCherry") %>%
    pivot_longer(cols = -c(targetName), names_to = "RawName", values_to = "RawVal") %>%
    inner_join(anno_df, by = "targetName") %>%
    mutate(AssayID = NA_character_, TargetName = targetName)
  if (level_name == "HybNorm") {
    df_long <- df_long %>%
      left_join(ic_df, by = "RawName") %>%
      mutate(Value = log2((RawVal + 0.1) / ScaleFactor)) %>%
      select(-IC_Count, -ScaleFactor)
  } else if (level_name == "Raw") {
    df_long <- df_long %>% mutate(Value = log2(RawVal + 0.1))
  } else {
    df_long <- df_long %>% mutate(Value = RawVal)
  }
  if (!is.null(lod_df)) {
    df_long <- df_long %>% left_join(lod_df, by = "TargetName")
  } else {
    df_long <- df_long %>% mutate(LOD = NA_real_)
  }
  df_long %>%
    select(RawName, AssayID, TargetName, UniProtID, Value, LOD) %>%
    standardize_long("NLS_P1_B1", level_name)
}

raw_nls_npq <- read_xlsx(upstream_path("raw/NLS_CIM_250603/20250603_SLM data analysis_NULISA_Analysis_Software_NPQ_Values_Annotation_LYQ.xlsx"), sheet = 1)
raw_nls_cnt <- read_xlsx(upstream_path("raw/NLS_CIM_250603/20250603_SLM data analysis_NULISA_Analysis_Software_NPQ_Values_Annotation_LYQ.xlsx"), sheet = 2)
anno_nls <- read_xlsx(upstream_path("raw/NLS_CIM_250603/20250603_SLM data analysis_NULISA_Analysis_Software_NPQ_Values_Annotation_LYQ.xlsx"), sheet = 4) %>%
  select(targetName, UniProt) %>% rename(UniProtID = UniProt)

ic_row_raw <- raw_nls_cnt %>% filter(targetName == "mCherry")
ic_vals_long <- ic_row_raw %>%
  pivot_longer(-targetName, names_to = "RawName", values_to = "IC_Count") %>%
  select(RawName, IC_Count)

ic_median <- median(ic_vals_long$IC_Count, na.rm = TRUE)
ic_vals_long <- ic_vals_long %>%
  mutate(ScaleFactor = IC_Count / ic_median)

nc_cols <- meta_sample %>% filter(Batch == "NLS_P1_B1", Sample == "BLK") %>% pull(RawName)
lod_raw <- calc_nulisa_lod(raw_nls_cnt, nc_cols, "Raw")
lod_hyb <- calc_nulisa_lod(raw_nls_cnt, nc_cols, "HybNorm", ic_scale_df = ic_vals_long)
lod_npq <- calc_nulisa_lod(raw_nls_npq, nc_cols, "NPQ")

df_nls_raw <- process_nulisa_step(raw_nls_cnt, "Count", "Raw", lod_raw, anno_nls)
df_nls_hyb <- process_nulisa_step(raw_nls_cnt, "Count", "HybNorm", lod_hyb, anno_nls, ic_vals_long)
df_nls_npq <- process_nulisa_step(raw_nls_npq, "NPQ", "NPQ", lod_npq, anno_nls)

long_data_list[["NLS_P1_B1"]] <- bind_rows(df_nls_raw, df_nls_hyb, df_nls_npq)

# 6. AAgAtlas antibody arrays ----
process_aag <- function(file_path, sheet_idx, aag_anno, bat_id, proc_level) {
    raw <- read_excel(file_path, sheet = sheet_idx) %>%
        select(-any_of(c("Column NO.", "Row NO."))) %>%
        inner_join(aag_anno, by = "Symbol") %>%
        pivot_longer(
            cols = -c(Symbol, TargetName, UniProtID), names_to = "RawName", values_to = "Intensity"
        ) %>%
        mutate(AssayID = NA_character_, Value = log2(Intensity)) %>%
        select(RawName, AssayID, TargetName, UniProtID, Value) %>%
        standardize_long(bat_id, proc_level)
}
aag_anno <- read_xlsx(upstream_path("raw/AAG_PEM_250103/AAgAtlas分子注释_YL.xlsx"), "Patch") %>% select(-Pos, -Rename, -`Vendor UniProtID`, -Note) %>% unique()
df_aag_raw_iga <- process_aag(upstream_path("raw/AAG_PEM_250103/ruimin/SNR数据_INSM.xlsx"), "635原始数据", aag_anno, "AAG_P1_B1", "Raw")
df_aag_snr_iga <- process_aag(upstream_path("raw/AAG_PEM_250103/ruimin/SNR数据_INSM.xlsx"), "635SNR数据", aag_anno, "AAG_P1_B1", "SNR")
long_data_list[["AAG_P1_B1"]] <- bind_rows(df_aag_raw_iga, df_aag_snr_iga)
df_aag_raw_igg <- process_aag(upstream_path("raw/AAG_PEM_250103/ruimin/SNR数据_INSM.xlsx"), "532原始数据", aag_anno, "AAG_P2_B1", "Raw")
df_aag_snr_igg <- process_aag(upstream_path("raw/AAG_PEM_250103/ruimin/SNR数据_INSM.xlsx"), "532SNR数据", aag_anno, "AAG_P2_B1", "SNR")
long_data_list[["AAG_P2_B1"]] <- bind_rows(df_aag_raw_igg, df_aag_snr_igg)

# 7. DIA mass spectrometry ----
process_dia <- function(file_path, bat_id) {
  if(grepl(".xlsx$", file_path)) raw <- read_xlsx(file_path) else raw <- fread(file_path)
  pg_col <- "Protein.Group"
  gene_col <- "Genes"
  meta_sub <- meta_sample %>% filter(Batch == bat_id)
  cols_to_keep <- c(pg_col, gene_col, intersect(colnames(raw), meta_sub$RawName))
  raw %>%
    select(all_of(cols_to_keep)) %>%
    rename(PG_Raw = !!sym(pg_col), Genes_Raw = !!sym(gene_col)) %>%
    mutate(
        UniProtID = clean_and_sort_str(PG_Raw), TargetName = clean_and_sort_str(Genes_Raw),
        AssayID = NA_character_
    ) %>%
    pivot_longer(cols = -c(PG_Raw, Genes_Raw, UniProtID, TargetName, AssayID), names_to = "RawName", values_to = "Intensity") %>%
    filter(!is.na(Intensity)) %>%
    mutate(Value = log2(Intensity)) %>%
    standardize_long(bat_id, "Intensity")
}

long_data_list[["DIA_P1_B1"]] <- process_dia(upstream_path("raw/DIA_AMS_241120/DIANN_FDU32.txt"), "DIA_P1_B1")
long_data_list[["DIA_P2_B1"]] <- process_dia(upstream_path("raw/DIA_TMO_241118/report.pg_matrix.xlsx"), "DIA_P2_B1")
long_data_list[["DIA_P3_B1"]] <- process_dia(upstream_path("raw/DIA_WLU_241129/WAN20241129gaohh_ShilmGroup_24minDIA_DIANN181_report.pg_matrix.tsv"), "DIA_P3_B1")
long_data_list[["DIA_P4_B1"]] <- process_dia(upstream_path("raw/DIA_WLU_241129/WAN20241129gaohh_ShilmGroup_24minDIA_OmniProt_DIANN181_report.pg_matrix.tsv"), "DIA_P4_B1")
long_data_list[["DIA_P5_B1"]] <- process_dia(upstream_path("raw/DIA_IPM_250305/20250305/report.pg_matrix.tsv"), "DIA_P5_B1")
long_data_list[["DIA_P5_B2"]] <- process_dia(upstream_path("raw/DIA_IPM_250423/FDU_Plasma_20250423/report.pg_matrix.tsv"), "DIA_P5_B2")

# 8. Integrate features, processing tiers and batch statistics ----
message(">>> Starting Part 6: Final Output Generation...")
target_levels <- c("ANML-SMP", "MedNormExt", "NPX", "NPQ", "SNR", "Intensity")
temp_long <- rbindlist(long_data_list, fill = TRUE)
temp_long <- temp_long[!is.na(Value)]
length(unique(temp_long$Batch))
length(unique(temp_long$ColName))

raw_stats <- temp_long[ProcessLevel %in% target_levels & Include == TRUE,
                      .(n_rows = uniqueN(paste(Platform, AssayID, TargetName, UniProtID))),
                      by = .(Batch, Platform)] %>%
    rename(Raw_Assay_Count = n_rows)
message(">>> Track 1 Stats calculated exclusively on Include == TRUE samples.")

temp_long[, Distinction_Key := case_when(
    Platform == "SOM" ~ AssayID, Platform %in% c("NLS", "AAG") ~ TargetName,
    Platform %in% c("DIA", "OLK") ~ fcoalesce(UniProtID, TargetName), TRUE ~ TargetName
)]
temp_long[is.na(Distinction_Key) | Distinction_Key == ""]

feat_meta_pre <- unique(temp_long[ProcessLevel %in% target_levels,
                                  .(Platform, AssayID, TargetName, UniProtID, Distinction_Key)])

uniprot_info <- as.data.table(uniprot_ref)[, .(UniProtID = Protein.ID, Protein_Full_Name = Protein.Fullname)]
uniprot_info <- unique(uniprot_info, by = "UniProtID")
feat_meta <- feat_meta_pre %>%
    left_join(uniprot_info, by = "UniProtID") %>%
    group_by(Platform, UniProtID) %>%
    mutate(
        Is_Tau_Primary = (UniProtID == "P10636" & TargetName %in% c("MAPT", "tTau")),
        Is_Tau_Variant = (UniProtID == "P10636" & grepl("pTau", TargetName))
    ) %>%
    arrange(desc(Is_Tau_Primary), Is_Tau_Variant, Distinction_Key) %>%
    mutate(
        Rank = row_number(),
        Suffix = Distinction_Key,
        UniqueID = case_when(
            Platform %in% c("DIA", "OLK") ~ UniProtID,
            is.na(UniProtID) ~ Suffix,
            Rank == 1 ~ UniProtID,
            TRUE ~ paste0(UniProtID, "_", Suffix)
        )
    ) %>%
    ungroup() %>%
    mutate(Is_Protein_Group = grepl("[|;]", UniProtID), Is_Unknown = is.na(UniProtID) | UniProtID == "") %>%
    as.data.table()

map_dt <- unique(feat_meta[, .(Platform, Distinction_Key, UniqueID, UniProtID)])
setkey(temp_long, Platform, Distinction_Key)
setkey(map_dt, Platform, Distinction_Key)

message(">>> Merging Metadata...")
merged_long <- merge(temp_long, map_dt, by = c("Platform", "Distinction_Key", "UniProtID"), all.x = TRUE)
long_df <- merged_long[
    ,
    .(
        Value = mean(Value, na.rm = TRUE), LOD = mean(LOD, na.rm = TRUE),
        RawName = paste(unique(RawName), collapse = ";"), AssayID = paste(unique(AssayID), collapse = ";"),
        TargetName = paste(unique(TargetName), collapse = ";"), Include = any(Include == TRUE, na.rm = TRUE)
    ),
    by = .(Platform, Batch, Sample, ColName, ProcessLevel, UniqueID, UniProtID)
]

tier_dict <- c(
    "Raw" = "Readout", "HybNorm" = "Baseline", "SNR" = "Baseline", "Intensity" = "Baseline",
    "Calibrate" = "Calibrated", "NPX" = "Calibrated", "NPQ" = "Calibrated", "ANML-SMP" = "Reshaped",
    "MedNormExt"= "Reshaped", "PlateScale"= "Intermediate", "MedNormInt"= "Intermediate"
)
long_df <- long_df %>% mutate(DataTier = tier_dict[ProcessLevel])

before_rows <- nrow(merged_long)
after_rows <- nrow(long_df)
duplicate_groups <- merged_long[, .N, by = .(Platform, Batch, Sample, ColName, ProcessLevel, UniqueID, UniProtID)][N > 1]
total_dup_rows <- sum(duplicate_groups$N)
n_groups <- nrow(duplicate_groups)
expected_reduction <- total_dup_rows - n_groups
message(sprintf("Rows before aggregation: %s, rows after aggregation: %s, observed reduction: %s", format(before_rows, big.mark=","), format(after_rows, big.mark=","), format(before_rows - after_rows, big.mark=",")))
message(sprintf("duplicate groups: %s, duplicate-group rows: %s, expected reduction: %s", format(n_groups, big.mark=","), format(total_dup_rows, big.mark=","), format(expected_reduction, big.mark=",")))
validation_ok <- before_rows - after_rows == expected_reduction
message(sprintf("Aggregation validation: %s (%s = %s - %s)", ifelse(validation_ok, "passed", "failed"),
                format(before_rows - after_rows, big.mark = ","), format(total_dup_rows, big.mark = ","),
                format(n_groups, big.mark = ",")))
if (!validation_ok) stop("Unexpected row loss during duplicate aggregation.")
dup_platforms <- duplicate_groups[, .(Groups = .N, Dup_Rows = sum(N)), by = Platform]
message("Duplicate sources:\n", paste(capture.output(print(dup_platforms)), collapse="\n"))
if(any(grepl(";", long_df$TargetName))) {warning("TargetName  contains semicolons unexpectedly")} else {message("✓ TargetName  contains no semicolons, as expected")}

message(">>> Calculating Detailed Batch Statistics...")
batch_features <- long_df[Include == TRUE, .(UniqueID = unique(UniqueID)), by = .(Batch, Platform)]
feat_type_map <- feat_meta[, .(UniqueID, Platform, Is_Protein_Group, Is_Unknown, UniProtID)] %>%
    unique() %>%
    mutate(
        Mapping_Type = case_when(
            Is_Unknown ~ "Unknown",
            Is_Protein_Group ~ "Group",
            Platform == "SOM" & grepl("_seq", UniqueID) ~ "Multi",
            Platform %in% c("NLS", "AAG") & grepl("_", UniqueID) ~ "Multi",
            TRUE ~ "Single"
        )
    ) %>% as.data.table()

calc_detailed_stats <- function(chunk, platform_name) {
    dt <- merge(chunk, feat_type_map[Platform == platform_name], by = "UniqueID", all.x = TRUE)
    dt[, n_assays_per_tgt := .N, by = UniProtID]
    list(
        Analysis_Assays = uniqueN(dt$UniqueID), Assay_1to1 = sum(dt$Mapping_Type == "Single"),
        Assay_Multi = sum(dt$Mapping_Type == "Multi"), Assay_Group = sum(dt$Mapping_Type == "Group"),
        Assay_Unknown = sum(dt$Mapping_Type == "Unknown"),
        Tgt_Total = uniqueN(dt$UniProtID[!is.na(dt$UniProtID)]),
        Tgt_Single = uniqueN(dt$UniProtID[!is.na(dt$UniProtID) & dt$n_assays_per_tgt == 1]),
        Tgt_Multi = uniqueN(dt$UniProtID[!is.na(dt$UniProtID) & dt$n_assays_per_tgt > 1]),
        Tgt_Group = uniqueN(dt$UniProtID[dt$Mapping_Type == "Group"])
    )
}

analysis_stats_list <- list()
for(bat in unique(batch_features$Batch)) {
    res <- calc_detailed_stats(batch_features[Batch == bat], batch_features[Batch == bat, unique(Platform)])
    res$Batch <- bat
    analysis_stats_list[[bat]] <- res
}
analysis_stats <- rbindlist(analysis_stats_list)

batch_stats_final <- raw_stats %>%
    left_join(analysis_stats) %>%
    select(
        Batch, `Total features` = Raw_Assay_Count, `Unique-target features` = Assay_1to1,
        `Co-targeting features` = Assay_Multi, `Protein-group features` = Assay_Group,
        `Unmapped features` = Assay_Unknown, `Total target proteins` = Tgt_Total,
        `Single-feature proteins` = Tgt_Single, `Multi-feature proteins` = Tgt_Multi,
        `Protein groups` = Tgt_Group
    ) %>%
    arrange(Batch)

# Recorded notebook output (reference only; regenerate by running this script):
# Imported batches: 17; imported sample columns: 684.
# Rows before aggregation: 15,604,193; rows after aggregation: 15,601,013; observed reduction: 3,180.
# Duplicate groups: 2,756; duplicate-group rows: 5,936; expected reduction: 3,180.
# Aggregation validation passed; duplicate source: AAG.

# 9. Export the article-specific release profile and metadata ----
# The release profile retains the selected M/Y/P/X/F/N measurements and the BLK/CAL/QC controls used elsewhere in the manuscript. Internal raw names, tube identifiers, barcodes, laboratory codes and unrelated samples are not exported.

valid_features <- unique(feat_meta$UniqueID[!feat_meta$Is_Protein_Group & !feat_meta$Is_Unknown])
message("Analysis features: ", length(valid_features))

dia_p04062 <- feat_meta %>%
    filter(Platform == "DIA", grepl("P04062", UniProtID)) %>%
    pull(UniqueID) %>%
    unique()
if (length(dia_p04062) == 1 && identical(dia_p04062, "P04062")) {
    message("DIA check passed: GBA/GBA1 was merged to P04062.")
} else if (length(dia_p04062) == 0) {
    warning("DIA check failed: P04062 was not found.")
} else {
    warning("DIA check failed: multiple P04062 forms were found: ", paste(dia_p04062, collapse = ", "))
}

if (any(grepl("_seq", feat_meta$UniqueID))) {
    message("SomaScan isoform check passed: sequence-specific identifiers were preserved.")
} else {
    warning("SomaScan isoform check failed: no sequence-specific identifier was found.")
}

tau_check <- feat_meta %>%
    filter(grepl("P10636", UniqueID)) %>%
    pull(UniqueID) %>%
    unique()
message("NULISA tau forms: ", length(tau_check))
print(tau_check)

release_core_samples <- c("M", "Y", "P", "X", "F", "N")
release_control_samples <- c("BLK", "CAL", "QC")

# Public profiles retain valid features, included reference samples and release controls.
release_long_df <- as_tibble(long_df) %>%
    filter(
        UniqueID %in% valid_features,
        (Include %in% TRUE & Sample %in% release_core_samples) | Sample %in% release_control_samples
    ) %>%
    select(-any_of(c("RawName", "Include"))) %>%
    arrange(Platform, Batch, Sample, ColName, ProcessLevel, UniqueID)

release_colnames <- unique(release_long_df$ColName)
release_batches <- unique(release_long_df$Batch)

release_sample_metadata <- meta_sample %>%
    filter(ColName %in% release_colnames) %>%
    select(Batch, Platform, ColName, Sample, Replicate) %>%
    distinct() %>%
    arrange(Batch, Sample, Replicate, ColName)

release_batch_metadata <- meta_batch %>%
    filter(Batch %in% release_batches) %>%
    select(Batch, Platform, Principle, Protocol, Protocol_A, Protocol_B, Protocol_Detail) %>%
    distinct() %>%
    arrange(Batch)

release_variance_metadata <- meta_variance %>%
    filter(Batch %in% release_batches) %>%
    arrange(Batch)

missing_metadata_columns <- setdiff(release_colnames, release_sample_metadata$ColName)
if (length(missing_metadata_columns) > 0) {
    stop(
        "Release profile columns missing from the public sample metadata: ",
        paste(missing_metadata_columns, collapse = ", ")
    )
}

fwrite(as.data.table(release_long_df), "data/protein_profiles_long.tsv.gz", sep = "\t", na = "NA")
fwrite(as.data.table(feat_meta), "data/feature_metadata.tsv.gz", sep = "\t", na = "NA")

write.xlsx(
    list(batch = release_batch_metadata, sample = release_sample_metadata, variance = release_variance_metadata),
    "data/study_metadata.xlsx", overwrite = TRUE, keepNA = TRUE, na.string = "NA"
)

cat("Release profile rows:", format(nrow(release_long_df), big.mark = ","), "\n")
cat("Release sample measurements:", length(release_colnames), "\n")
cat("Release batches:", length(release_batches), "\n")

# Recorded notebook output (reference only; regenerate by running this script):
# Analysis features: 13,836
# DIA check passed: GBA/GBA1 was merged to P04062.
# SomaScan isoform check passed: sequence-specific identifiers were preserved.
# NULISA tau forms: 6 (P10636, pTau181, pTau205, pTau212, pTau217 and pTau231).
# Release profile rows: 9,030,107; release sample measurements: 488; release batches: 17.

# 10. Platform-specific detection status ----
# Detection is deliberately platform-specific and is not interpreted as a common sensitivity scale. The anchor tiers are:
# - DIA and AAgAtlas: `Baseline`
# - Olink, NULISA and SomaScan: `Calibrated`
# A sample group passes when more than half of its included technical measurements exceed the platform-specific criterion. `IsDetected` is the primary status field. `IsAboveLoD` is retained as a backward-compatible alias for existing downstream code.

valid_features <- unique(feat_meta$UniqueID[!feat_meta$Is_Protein_Group & !feat_meta$Is_Unknown])
long_df_filter <- long_df %>%
    filter(UniqueID %in% valid_features, Include == TRUE, DataTier %in% c("Baseline", "Calibrated", "Reshaped"))

# A group passes when more than half of its expected replicates are detected.
calc_detection_status <- function(long_data, meta_sample) {
    analysis_samples <- c("M", "Y", "P", "X", "F", "N")
    long_data <- long_data %>% filter(Sample %in% analysis_samples)
    meta_core_sample <- meta_sample %>% filter(Sample %in% analysis_samples)

    group_sizes <- meta_core_sample %>%
        filter(Include == TRUE) %>%
        group_by(Batch, Sample) %>%
        summarize(expected_n = n(), .groups = "drop")

    batch_expected_total <- group_sizes %>%
        group_by(Batch) %>%
        summarize(ExpectedN = sum(expected_n), .groups = "drop")

    df_det_check <- long_data %>%
        mutate(
            IsDetectedMeasurement = case_when(
                Platform == "AAG" ~ Value >= log2(3),
                Platform == "DIA" ~ !is.na(Value) & Value > 0,
                Platform %in% c("OLK", "NLS", "SOM") & !is.na(LOD) ~ Value > LOD,
                TRUE ~ FALSE
            )
        )

    global_stats <- df_det_check %>%
        group_by(DataTier, ProcessLevel, Batch, Platform, UniqueID) %>%
        summarize(DetectedN = sum(IsDetectedMeasurement, na.rm = TRUE), .groups = "drop") %>%
        left_join(batch_expected_total, by = "Batch")

    status_summary <- df_det_check %>%
        group_by(DataTier, ProcessLevel, Batch, Platform, UniqueID, Sample) %>%
        summarize(DetectedN = sum(IsDetectedMeasurement, na.rm = TRUE), .groups = "drop") %>%
        left_join(group_sizes, by = c("Batch", "Sample")) %>%
        mutate(is_group_pass = (DetectedN / expected_n) > 0.5) %>%
        pivot_wider(
            id_cols = c(DataTier, ProcessLevel, Batch, Platform, UniqueID), names_from = Sample,
            values_from = is_group_pass, values_fill = FALSE
        )

    for (sample_name in c("M", "Y", "P", "X", "F", "N")) {
        if (!sample_name %in% colnames(status_summary)) {
            status_summary[[sample_name]] <- FALSE
        }
    }

    status_summary %>%
        mutate(
            PassedGroup = as.integer(M) + as.integer(F) + as.integer(Y) +
                as.integer(P) + as.integer(X) + as.integer(N),
            JumpNote = case_when(
                PassedGroup >= 1 & !(M | F | N) ~ "ArtifactualJump",
                PassedGroup == 1 & (M | F | N) ~ "SingleSourceOnly",
                PassedGroup > 0 ~ "RobustDetection",
                TRUE ~ "NotDetected"
            )
        ) %>%
        left_join(global_stats, by = c("DataTier", "ProcessLevel", "Batch", "Platform", "UniqueID")) %>%
        select(
            DataTier, ProcessLevel, Batch, Platform, UniqueID, M, F, Y, P, X, N, ExpectedN, DetectedN,
            PassedGroup, JumpNote
        )
}

detection_status_multi <- calc_detection_status(long_df_filter, meta_sample) %>%
    mutate(IsDetected = PassedGroup > 0, IsAboveLoD = IsDetected)

print(table(detection_status_multi[, c("DataTier", "JumpNote", "Batch")]))
print(table(detection_status_multi[, c("DataTier", "IsDetected", "Batch")]))

# Anchor tiers follow the original LoD analysis: Baseline for DIA/AAgAtlas and Calibrated for affinity platforms.
detection_status <- detection_status_multi %>%
    filter(
        (Platform %in% c("DIA", "AAG") & DataTier == "Baseline") |
            (Platform %in% c("OLK", "SOM", "NLS") & DataTier == "Calibrated")
    ) %>%
    rename(AnchorTier = DataTier)

message("Anchor-tier extraction check:")
print(table(Platform = detection_status$Platform, AnchorTier = detection_status$AnchorTier))

expected_unique_ids <- length(unique(long_df_filter$UniqueID))
actual_unique_ids <- length(unique(detection_status$UniqueID))
if (expected_unique_ids == actual_unique_ids) {
    message(
        sprintf(
            "Validation passed: all %d features received an anchor-tier detection status.",
            actual_unique_ids
        )
    )
} else {
    warning(
        sprintf(
            "Feature-count mismatch: expected %d, observed %d.",
            expected_unique_ids,
            actual_unique_ids
        )
    )
}

fwrite(as.data.table(detection_status), "results/detection_status.tsv.gz", sep = "\t", na = "NA")

# Recorded notebook output (reference only; regenerate by running this script):
# Anchor-tier extraction: AAG Baseline 2,006; DIA Baseline 8,303; NLS Calibrated 127;
# OLK Calibrated 11,540; SOM Calibrated 42,064.
# Validation passed: all 13,158 features received an anchor-tier detection status.

# 11. Detection-summary source data ----
df_counts <- detection_status %>%
    mutate(Passed_Count = as.integer(M) + as.integer(F) + as.integer(Y) +
        as.integer(P) + as.integer(X) + as.integer(N)) %>%
    arrange(Batch)

analyzed_counts <- df_counts %>%
    count(Batch, name = "Analyzed features")

threshold_wide <- expand.grid(Batch = unique(df_counts$Batch), Threshold = 1:6) %>%
    rowwise() %>%
    mutate(
        Feature_Count = sum(df_counts$Batch == Batch & df_counts$Passed_Count >= Threshold),
        Feature_Percent = 100 * Feature_Count / sum(df_counts$Batch == Batch),
        Display = sprintf(
            "%s (%.1f%%)",
            format(Feature_Count, big.mark = ",", scientific = FALSE, trim = TRUE),
            Feature_Percent
        )
    ) %>%
    ungroup() %>%
    select(Batch, Threshold, Display) %>%
    pivot_wider(names_from = Threshold, values_from = Display, names_prefix = "k = ") %>%
    rename(
        `Detected in ≥ 1 group` = `k = 1`,
        `Detected in ≥ 2 groups` = `k = 2`,
        `Detected in ≥ 3 groups` = `k = 3`,
        `Detected in ≥ 4 groups` = `k = 4`,
        `Detected in ≥ 5 groups` = `k = 5`,
        `Detected in 6 groups` = `k = 6`
    )

st1_table <- batch_stats_final %>%
    left_join(analyzed_counts, by = "Batch") %>%
    left_join(threshold_wide, by = "Batch") %>%
    select(
        Batch, `Total features`, `Unique-target features`, `Co-targeting features`,
        `Protein-group features`, `Unmapped features`, `Total target proteins`,
        `Single-feature proteins`, `Multi-feature proteins`, `Protein groups`,
        `Analyzed features`, `Detected in ≥ 1 group`, `Detected in ≥ 2 groups`,
        `Detected in ≥ 3 groups`, `Detected in ≥ 4 groups`,
        `Detected in ≥ 5 groups`, `Detected in 6 groups`
    ) %>%
    arrange(Batch)

stopifnot(
    !anyDuplicated(st1_table$Batch),
    setequal(st1_table$Batch, batch_stats_final$Batch),
    !anyNA(st1_table)
)

write.xlsx(list(Analytical_feature_statistics = st1_table), "tables/SourceData_AnalyticalFeatureStatistics.xlsx", overwrite = TRUE, keepNA = TRUE, na.string = "NA")

# 12. Final validation summary ----
profile_outputs <- c(
    "data/protein_profiles_long.tsv.gz", "data/feature_metadata.tsv.gz", "data/study_metadata.xlsx",
    "results/detection_status.tsv.gz", "tables/SourceData_AnalyticalFeatureStatistics.xlsx"
)
stopifnot(all(file.exists(profile_outputs)))

cat("Internal batches processed:", n_distinct(long_df$Batch), "\n")
cat("Release batches:", n_distinct(release_long_df$Batch), "\n")
cat("Analysis features:", n_distinct(valid_features), "\n")
cat("Release profile rows:", format(nrow(release_long_df), big.mark = ","), "\n")
cat("Detection-status rows:", format(nrow(detection_status), big.mark = ","), "\n")
cat("All expected outputs are present.\n")

# Recorded notebook output (reference only; regenerate by running this script):
# Internal batches processed: 17
# Release batches: 17
# Analysis features: 13,836
# Release profile rows: 9,030,107
# Detection-status rows: 64,040
# All expected outputs are present.
