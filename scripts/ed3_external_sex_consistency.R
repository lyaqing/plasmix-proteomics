# Extended Data Figure 3 | External-cohort sex-effect consistency

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "dplyr", "tidyr", "tibble", "stringr", "readxl", "readr", "ggplot2", "ggrepel", "patchwork", "grid", "ggh4x", "ggpp", "ggpubr", "openxlsx", "showtext"))
source("utils/figure_style.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
label_style <- list(size = 12, face = "bold")
calc_fdr <- function(p_vals) {
    p_adj <- rep(NA, length(p_vals))
    valid_idx <- !is.na(p_vals)
    if(sum(valid_idx) > 0) {
        p_adj[valid_idx] <- p.adjust(p_vals[valid_idx], method = "BH")
    }
    return(p_adj)
}
canonicalize_uniprot <- function(u_vec) {
    sapply(u_vec, function(x) {
        if (is.na(x) || x == "") return(NA_character_)
        x_clean <- gsub("-[0-9]+", "", x)
        parts <- strsplit(x_clean, "[|_]")[[1]]
        paste(sort(unique(trimws(parts))), collapse = ";")
    }, USE.NAMES = FALSE)
}

uniprot_ref <- read_tsv(upstream_path("references", "uniprotkb_AND_model_organism_9606_2026_04_01.tsv.gz"), show_col_types = FALSE) %>%
    filter(Reviewed == "reviewed") %>%
    transmute(Protein.ID = Entry, Protein.Name = gsub("_HUMAN$", "", `Entry Name`), Gene.Symbol = `Gene Names (primary)`) %>%
    distinct(Gene.Symbol, .keep_all = TRUE)

map_3072 <- read_excel(upstream_path("references", "olink-explore-3072-ht-assay-list.xlsx"), sheet = "3072") %>%
    transmute(Gene_3072 = `Gene name`, UniProt = canonicalize_uniprot(`UniProt ID`)) %>%
    filter(!is.na(UniProt), UniProt != "") %>% distinct()
map_ht <- read_excel(upstream_path("references", "olink-explore-3072-ht-assay-list.xlsx"), sheet = "HT") %>%
    transmute(Gene_HT = `Gene name`, UniProt = canonicalize_uniprot(`UniProt ID`)) %>%
    filter(!is.na(UniProt), UniProt != "") %>% distinct()
mapping_ledger <- full_join(map_3072, map_ht, by = "UniProt")

# 1. External-cohort inputs and standardization ----
# Cohort 1: Nature UKBiobank (Olink 3K)
nature_ukb <- read_excel(upstream_path("cohort", "nature_ukb_iceland", "Nature_ST14_ST15.xlsx"), sheet = "ST14_olink_age_sex_corr", skip = 3) %>%
    mutate(
        uniprot = case_when(
            uniprot == "NTproBNP" & (is.na(gene_name) | gene_name == "NA") ~ "NT-proBNP",
            uniprot == "A0A0B4J2D5" & gene_name == "C21orf33" ~ "P0DPI2",
            TRUE ~ uniprot
        ),
        gene_name = case_when(
            uniprot == "NT-proBNP" & (is.na(gene_name) | gene_name == "NA") ~ "NT-proBNP",
            uniprot == "P0DPI2" & gene_name == "C21orf33" ~ "GATD3",
            TRUE ~ gene_name
        )
    ) %>%
    mutate(across(everything(), ~ifelse(. == "NA", NA, .))) %>%
    select(UniProt_Raw = uniprot, Gene.Symbol = gene_name, beta_raw = beta_sex_bi, pval_raw = pval_sex_bi) %>%
    filter(!is.na(beta_raw) & !is.na(pval_raw)) %>%
    mutate(UniProt = canonicalize_uniprot(UniProt_Raw)) %>%
    filter(UniProt %in% mapping_ledger$UniProt) %>%
    mutate(
        beta_raw = gsub(",", ".", beta_raw),
        pval_raw = gsub(",", ".", pval_raw),
        Platform = "Olink_3K",
        Cohort = "UKBiobank",
        Source = "Nature_2023",
        Effect_Size = -1 * as.numeric(beta_raw),
        Pval = as.numeric(pval_raw),
        FDR = calc_fdr(Pval),
        Classification = ifelse(FDR < 0.05, "Significant", "Non-significant")
    ) %>%
    select(UniProt, Gene.Symbol, Effect_Size, Pval, FDR, Platform, Source, Cohort, Classification)

# Cohort 2: Nature Iceland (SomaScan 5K)
nature_ice <- read_excel(upstream_path("cohort", "nature_ukb_iceland", "Nature_ST14_ST15.xlsx"), sheet = "ST15_somascan_age_sex_corr", skip = 3) %>%
    select(Gene.Symbol = gene_name, beta_raw = `beta sex no normalization`, pval_raw = `pval sex no normalization`) %>%
    left_join(uniprot_ref, by = "Gene.Symbol") %>%
    filter(!is.na(Protein.ID), !is.na(beta_raw), !is.na(pval_raw)) %>%
    rename(UniProt = Protein.ID) %>%
    mutate(
        Platform = "SomaScan_5K",
        Cohort = "Iceland",
        Source = "Nature_2023",
        Effect_Size = -1 * as.numeric(beta_raw),
        Pval = as.numeric(pval_raw),
        FDR = calc_fdr(Pval),
        Classification = ifelse(FDR < 0.05, "Significant", "Non-significant")
    ) %>%
    select(UniProt, Gene.Symbol, Effect_Size, Pval, FDR, Platform, Source, Cohort, Classification)

ht_dict_clean <- mapping_ledger %>% select(UniProt, Gene.Symbol = Gene_HT) %>% filter(!is.na(Gene.Symbol) & Gene.Symbol != "")
process_science <- function(file_path, cohort_name) {
    read_excel(file_path, sheet = 2) %>%
        as.data.frame() %>%
        filter(term == "SexM") %>%
        select(Gene.Symbol = Protein, estimate, p.value) %>%
        left_join(ht_dict_clean, by = "Gene.Symbol") %>%
        mutate(
            Platform = "Olink_HT",
            Cohort = cohort_name,
            Source = "Science_2025",
            Effect_Size = as.numeric(estimate),
            Pval = as.numeric(p.value),
            FDR = calc_fdr(Pval),
            Classification = ifelse(FDR < 0.05, "Significant", "Non-significant")
        ) %>%
        select(UniProt, Gene.Symbol, Effect_Size, Pval, FDR, Platform, Source, Cohort, Classification)
}
science_wellness <- process_science(upstream_path("cohort", "science_bamse_wellness", "science.adx2678_data_s2.xlsx"), "Wellness")
science_bamse <- process_science(upstream_path("cohort", "science_bamse_wellness", "science.adx2678_data_s3.xlsx"), "BAMSE")

print(dim(nature_ukb));print(length(unique(nature_ukb$UniProt)))
print(dim(nature_ice));print(length(unique(nature_ice$UniProt)))
print(dim(science_wellness));print(length(unique(science_wellness$UniProt)))
print(dim(science_bamse));print(length(unique(science_bamse$UniProt)))
master_long <- bind_rows(nature_ukb, nature_ice, science_wellness, science_bamse)

duplicates_in_raw <- master_long %>% group_by(UniProt, Cohort, Platform, Source) %>% filter(n() > 1) %>% arrange(Cohort, UniProt) %>% ungroup()
n_dup_groups <- duplicates_in_raw %>% distinct(UniProt, Cohort) %>% nrow()
cat(sprintf("\n[Diagnostic] master_long contains %d groups (%d rows) with probe redundancy or mapping overlap.\n", n_dup_groups, nrow(duplicates_in_raw)))

# 2. Tier classification ----

fc_threshold <- 0.05
master_long_tagged <- master_long %>%
    mutate(
        Is_Strong_Sig = (Classification == "Significant" & abs(Effect_Size) >= fc_threshold),
        Direction_Vote = case_when(
            Is_Strong_Sig & Effect_Size > 0 ~ 1,
            Is_Strong_Sig & Effect_Size < 0 ~ -1,
            TRUE ~ 0
        )
    )

tier_summary <- master_long_tagged %>%
    group_by(UniProt) %>%
    summarize(
        n_strong_sig = n_distinct(Cohort[Is_Strong_Sig == TRUE], na.rm = TRUE),
        is_conflict  = any(Direction_Vote == 1, na.rm = TRUE) & any(Direction_Vote == -1, na.rm = TRUE),
        has_olink    = any(Is_Strong_Sig & Platform %in% c("Olink_3K", "Olink_HT"), na.rm = TRUE),
        has_soma     = any(Is_Strong_Sig & Platform == "SomaScan_5K", na.rm = TRUE),
        is_cross_platform = (has_olink & has_soma),
        .groups = "drop"
    ) %>%
    mutate(
        Tier = case_when(
            is_conflict ~ "Conflict",
            n_strong_sig >= 3 & is_cross_platform ~ "Tier1",
            (n_strong_sig >= 3 & !is_cross_platform) | (n_strong_sig == 2 & is_cross_platform) ~ "Tier2",
            (n_strong_sig == 2 & !is_cross_platform) | (n_strong_sig == 1) ~ "Tier3",
            TRUE ~ "NotSig"
        )
    )
sex_ext_df <- master_long %>%
    left_join(tier_summary %>% select(UniProt, Tier), by = "UniProt") %>%
    left_join(uniprot_ref %>% select(Protein.ID, Ref_Gene = Gene.Symbol, Protein.Name), by = c("UniProt" = "Protein.ID")) %>%
    mutate(Gene.Symbol = coalesce(na_if(Gene.Symbol, ""), Ref_Gene, Protein.Name)) %>%
    select(UniProt, Gene.Symbol, Effect_Size, Pval, FDR, Platform, Source, Cohort, Classification, Tier, Protein.Name)

tier_cross_check <- sex_ext_df %>%
    filter(!is.na(UniProt), !is.na(Tier)) %>%
    distinct(UniProt, Tier) %>%
    count(UniProt, name = "n_tiers") %>%
    filter(n_tiers > 1)
print(tier_cross_check)

print(table(unique(sex_ext_df[, c('UniProt', 'Tier')])$Tier))
# Conflict   NotSig    Tier1    Tier2    Tier3
#       27     4118      100      275     1689

fwrite(sex_ext_df, "results/external_sex_effects_assay_level.tsv.gz", sep = "\t", na = "NA")

# 3. Consistency of two large-scale studies (Tier Overview) ----
# Collapse duplicated assay records within each cohort/platform to one UniProt-level estimate.
# Do not aggregate UKBiobank with Iceland or Wellness with BAMSE at the source level.
cohort_level_df <- sex_ext_df %>%
    filter(!is.na(UniProt), !is.na(Effect_Size)) %>%
    group_by(UniProt, Source, Cohort, Platform) %>%
    summarize(
        Gene.Symbol = if (all(is.na(Gene.Symbol) | Gene.Symbol == "")) NA_character_ else as.character(na.omit(Gene.Symbol[Gene.Symbol != ""])[1]),
        Beta = median(Effect_Size, na.rm = TRUE),
        Pval_min = suppressWarnings(min(Pval, na.rm = TRUE)),
        FDR_min = suppressWarnings(min(FDR, na.rm = TRUE)),
        Strong_Class = ifelse(
            any(FDR < 0.05 & abs(Effect_Size) >= fc_threshold, na.rm = TRUE),
            "Significant",
            "Non-significant"
        ),
        .groups = "drop"
    ) %>%
    left_join(tier_summary %>% select(UniProt, Tier), by = "UniProt")

stopifnot(
    cohort_level_df %>%
        count(UniProt, Source, Cohort, Platform) %>%
        filter(n > 1) %>%
        nrow() == 0
)

external_sex_effects_long <- cohort_level_df %>%
    transmute(UniProt, Gene.Symbol, Effect_Size = Beta, Pval = Pval_min, FDR = FDR_min, Platform, Source, Cohort, Classification = ifelse(FDR_min < 0.05, "Significant", "Non-significant"), Tier)
external_sex_tier_summary <- tier_summary %>%
    left_join(sex_ext_df %>% group_by(UniProt) %>% summarize(Gene.Symbol = if (all(is.na(Gene.Symbol) | Gene.Symbol == "")) NA_character_ else first(na.omit(Gene.Symbol[Gene.Symbol != ""])), .groups = "drop"), by = "UniProt")
fwrite(external_sex_effects_long, "results/external_sex_effects_long.tsv.gz", sep = "\t", na = "NA")
fwrite(external_sex_tier_summary, "results/external_sex_tier_summary.tsv.gz", sep = "\t", na = "NA")

dup_direction_check <- sex_ext_df %>%
    filter(!is.na(UniProt), !is.na(Effect_Size)) %>%
    group_by(UniProt, Source, Cohort, Platform) %>%
    summarize(
        n_records = n(),
        n_pos = sum(Effect_Size > 0, na.rm = TRUE),
        n_neg = sum(Effect_Size < 0, na.rm = TRUE),
        has_direction_conflict = n_pos > 0 & n_neg > 0,
        .groups = "drop"
    ) %>%
    filter(n_records > 1 | has_direction_conflict)

cat("\n=== [Diagnostic 0] Within-cohort probe duplication and directional conflict ===\n")
print(dup_direction_check %>% filter(has_direction_conflict == TRUE))

source_summary_df <- cohort_level_df %>%
    group_by(UniProt, Source) %>%
    summarize(
        Gene.Symbol = ifelse(
            all(is.na(Gene.Symbol) | Gene.Symbol == ""),
            NA_character_,
            first(na.omit(Gene.Symbol))
        ),
        Tier = first(Tier),
        Beta_source = mean(Beta, na.rm = TRUE),
        n_cohorts_source = n_distinct(Cohort),
        .groups = "drop"
    )

plot_data_a <- source_summary_df %>%
    select(UniProt, Gene.Symbol, Tier, Source, Beta_source, n_cohorts_source) %>%
    pivot_wider(
        names_from = Source,
        values_from = c(Beta_source, n_cohorts_source),
        names_sep = "__"
    ) %>%
    filter(
        !is.na(Beta_source__Nature_2023),
        !is.na(Beta_source__Science_2025)
    ) %>%
    mutate(
        Tier = factor(Tier, levels = c("NotSig", "Conflict", "Tier3", "Tier2", "Tier1"))
    ) %>%
    arrange(Tier)

data_bg_a <- plot_data_a %>% filter(Tier == "NotSig")
data_fg_a <- plot_data_a %>% filter(Tier != "NotSig")

all_values_a <- c(plot_data_a$Beta_source__Nature_2023, plot_data_a$Beta_source__Science_2025)
min_val_a <- min(all_values_a, na.rm = TRUE)
max_val_a <- max(all_values_a, na.rm = TRUE)
padding_a <- (max_val_a - min_val_a) * 0.05
common_limits_a <- c(min_val_a - padding_a, max_val_a + padding_a)
range_len_a <- max_val_a - min_val_a + 2 * padding_a

stats_df_a <- plot_data_a %>%
    group_by(Tier) %>%
    summarize(
        n = n(),
        r_p = ifelse(n >= 3, cor(Beta_source__Nature_2023, Beta_source__Science_2025, method = "pearson", use = "complete.obs"), NA_real_),
        r_s = ifelse(n >= 3, cor(Beta_source__Nature_2023, Beta_source__Science_2025, method = "spearman", use = "complete.obs"), NA_real_),
        .groups = "drop"
    ) %>%
    arrange(desc(Tier)) %>%
    mutate(
        parse_label = case_when(
            Tier == "NotSig" | is.na(r_p) | n < 3 ~ sprintf("'%s' ~ '(' * italic(n) == %d * ')'", Tier, n),
            TRUE ~ sprintf("'%s' ~ '(' * italic(n) == %d * ',' ~ italic(r) == %.2f * ',' ~ italic(rho) == %.2f * ')'", Tier, n, r_p, r_s)
        ),
        x_pos = min_val_a + range_len_a * 0.3,
        y_pos = (min_val_a - padding_a) + (n() - row_number() + 0.75) * (range_len_a * 0.05)
    )

label_data_p1 <- plot_data_a %>% arrange(desc((abs(Beta_source__Nature_2023) + abs(Beta_source__Science_2025)) / 2)) %>% slice_head(n = 10)
p1 <- ggplot(plot_data_a, aes(x = Beta_source__Nature_2023, y = Beta_source__Science_2025)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey40") +
    geom_point(data = data_bg_a, aes(color = "NotSig"), alpha = 0.3, size = 0.5) +
    geom_point(data = data_fg_a, aes(color = Tier), alpha = 0.7, size = 1.2) +
    geom_text_repel(data = label_data_p1, aes(label = Gene.Symbol), nudge_x = -0.1, size = 2.8, color = "black", max.overlaps = 20, box.padding = 0.5, seed = 2026) +
    geom_point(data = stats_df_a, aes(x = x_pos - range_len_a * 0.01, y = y_pos + range_len_a * 0.001, color = Tier), size = 1.5, alpha = 1, show.legend = FALSE) +
    geom_text(data = stats_df_a, aes(x = x_pos + range_len_a * 0.02, y = y_pos, label = parse_label, color = Tier), parse = TRUE, hjust = 0, size = 2.8, show.legend = FALSE) +
    scale_color_manual(values = c("Tier1"="#E64B35", "Tier2"="#4DBBD5", "Tier3"="#00A087", "Conflict"="#F39B7F", "NotSig"="grey80")) +
    labs(x = expression(beta[sex] ~ " (UKBiobank & Iceland)"), y = expression(beta[sex] ~ " (Wellness & BAMSE)")) +
    plasmix_theme +
    scale_x_continuous(n.breaks = 6, expand = expansion(mult = c(0.0, 0.05))) +
    scale_y_continuous(n.breaks = 6, expand = expansion(mult = c(0.0, 0.05))) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35), axis.line = element_blank(), panel.grid.major = element_blank(), legend.position = "none") +
    coord_cartesian(xlim = common_limits_a, ylim = common_limits_a)

# 4. Intra-study replicability (UniProt-level cohort pair) ----
darken_hex <- function(hex_colors, factor = 0.7) {
    sapply(hex_colors, function(h) {
        v <- col2rgb(h)
        rgb(floor(pmin(255, pmax(0, v[1] * factor))), floor(pmin(255, pmax(0, v[2] * factor))), floor(pmin(255, pmax(0, v[3] * factor))), maxColorValue = 255)
    })
}

draw_consistency_plot <- function(cohort_level_df, source_x, source_y, x_label, y_label) {
    df_x <- cohort_level_df %>%
        filter(Cohort == source_x) %>%
        select(UniProt, Gene.Symbol, Beta_X = Beta, FDR_X = FDR_min) %>%
        mutate(Class_X = ifelse(FDR_X < 0.05, "Significant", "Non-significant"))
    df_y <- cohort_level_df %>%
        filter(Cohort == source_y) %>%
        select(UniProt, Beta_Y = Beta, FDR_Y = FDR_min) %>%
        mutate(Class_Y = ifelse(FDR_Y < 0.05, "Significant", "Non-significant"))
    plot_df <- inner_join(df_x, df_y, by = "UniProt")
    stopifnot(plot_df %>% count(UniProt) %>% filter(n > 1) %>% nrow() == 0)
    spec_x <- paste(source_x, "specific")
    spec_y <- paste(source_y, "specific")
    plot_df <- plot_df %>%
        mutate(Consistency = case_when(
            Class_X == "Significant" & Class_Y == "Significant" & sign(Beta_X) == sign(Beta_Y) ~ "Consistent",
            Class_X == "Significant" & Class_Y == "Significant" & sign(Beta_X) != sign(Beta_Y) ~ "Conflict",
            Class_X == "Significant" & Class_Y != "Significant" ~ spec_x,
            Class_Y == "Significant" & Class_X != "Significant" ~ spec_y,
            TRUE ~ "Neither significant"
        )) %>%
        mutate(Consistency = factor(Consistency, levels = c("Neither significant", "Conflict", spec_x, spec_y, "Consistent"))) %>%
        arrange(Consistency)
    my_colors <- c("Consistent" = "#2ca02c", "Conflict" = "#d62728", "Neither significant" = "grey85")
    my_colors[spec_x] <- "#ff7f0e"
    my_colors[spec_y] <- "#1f77b4"
    text_colors <- darken_hex(my_colors, factor = 0.9)
    all_vals <- c(plot_df$Beta_X, plot_df$Beta_Y)
    min_v <- min(all_vals, na.rm = TRUE)
    max_v <- max(all_vals, na.rm = TRUE)
    pad <- (max_v - min_v) * 0.05
    limits <- c(min_v - pad, max_v + pad)
    real_rng <- limits[2] - limits[1]
    stats_df <- plot_df %>%
        group_by(Consistency) %>%
        summarize(
            n = n_distinct(UniProt),
            r_p = ifelse(n >= 3, cor(Beta_X, Beta_Y, method = "pearson", use = "complete.obs"), NA_real_),
            r_s = ifelse(n >= 3, cor(Beta_X, Beta_Y, method = "spearman", use = "complete.obs"), NA_real_),
            .groups = "drop"
        ) %>%
        mutate(Consistency = factor(Consistency, levels = c("Consistent", spec_x, spec_y, "Conflict", "Neither significant"))) %>%
        arrange(Consistency) %>%
        mutate(
            pt_color = my_colors[as.character(Consistency)],
            txt_color = text_colors[as.character(Consistency)],
            parse_label = case_when(
                Consistency == "Neither significant" | is.na(r_p) | n < 3 ~ sprintf("'%s' ~ '(' * italic(n) == %d * ')'", Consistency, n),
                TRUE ~ sprintf("'%s' ~ '(' * italic(n) == %d * ',' ~ italic(r) == %.2f * ',' ~ italic(rho) == %.2f * ')'", Consistency, n, r_p, r_s)
            ),
            x_pos = limits[1] + real_rng * 0.03,
            y_pos = limits[2] - real_rng * 0.03 - (row_number() - 1) * (real_rng * 0.05)
        )
    label_data <- plot_df %>% filter(Consistency == "Consistent") %>% arrange(desc(abs(Beta_X + Beta_Y))) %>% slice_head(n = 10)
    p <- ggplot(plot_df, aes(x = Beta_X, y = Beta_Y)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
        geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "grey40") +
        geom_point(aes(color = Consistency), alpha = 0.6, size = 1.2, show.legend = FALSE) +
        geom_text_repel(data = label_data, aes(label = Gene.Symbol), size = 2.8, color = "black", max.overlaps = 20, box.padding = 0.5, seed = 42) +
        geom_point(data = stats_df, aes(x = x_pos, y = y_pos + real_rng * 0.001, color = I(pt_color)), size = 1.5, alpha = 1) +
        geom_text(data = stats_df, aes(x = x_pos + real_rng * 0.025, y = y_pos, label = parse_label, color = I(txt_color)), parse = TRUE, hjust = 0, size = 2.8) +
        scale_color_manual(values = my_colors) +
        scale_x_continuous(limits = limits, expand = c(0.02, 0.02), breaks = c(-2, 0, 2, 4, 6, 8)) +
        scale_y_continuous(limits = limits, expand = c(0.02, 0.02), breaks = c(-2, 0, 2, 4, 6, 8)) +
        labs(x = x_label, y = y_label) +
        plasmix_theme +
        theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35), axis.line = element_blank(), panel.grid.major = element_blank(), legend.position = "none")
    return(list(plot = p, plot_df = plot_df))
}

res_p2 <- draw_consistency_plot(cohort_level_df, "UKBiobank", "Iceland", expression(beta[sex] ~ " (UK Biobank)"), expression(beta[sex] ~ " (Iceland)"))
p2 <- res_p2$plot
res_p3 <- draw_consistency_plot(cohort_level_df, "Wellness", "BAMSE", expression(beta[sex] ~ " (Wellness)"), expression(beta[sex] ~ " (BAMSE)"))
p3 <- res_p3$plot

selected_batches <- c("DIA_P1_B1", "DIA_P2_B1", "DIA_P3_B1", "DIA_P4_B1", "DIA_P5_B1", "DIA_P5_B2", "OLK_P2_B1", "OLK_P2_B2", "SOM_P1_B1", "SOM_P1_B2", "SOM_P2_B1", "SOM_P2_B2")
meta_batch_ht <- read_xlsx("data/study_metadata.xlsx", sheet = "batch") %>% filter(Batch %in% selected_batches)
if (n_distinct(meta_batch_ht$Batch) != length(selected_batches)) stop("study_metadata.xlsx does not contain all 12 high-throughput batches")
batch_palette <- if (exists("batch_color")) batch_color else setNames(scales::hue_pal()(length(selected_batches)), selected_batches)
missing_batch_colors <- setdiff(selected_batches, names(batch_palette))
if (length(missing_batch_colors)) batch_palette[missing_batch_colors] <- scales::hue_pal()(length(missing_batch_colors))

dea_df_multi <- fread("results/dea_df_multi.tsv.gz")
plasmix_mf_dep <- dea_df_multi %>%
    filter(Pair == "M/F", Batch %in% selected_batches,
           (Platform == "DIA" & DataTier == "Baseline") |
           (Platform %in% c("OLK", "SOM") & DataTier == "Calibrated")) %>%
    mutate(UniProt = tstrsplit(UniqueID, "_", fixed = TRUE)[[1]])

all_tier1_proteins <- sex_ext_df %>% filter(Tier == "Tier1") %>% pull(UniProt) %>% unique()
n_total_tier1 <- length(all_tier1_proteins)

tier1_consensus <- sex_ext_df %>%
    filter(Tier == "Tier1") %>%
    rename(ChortPlatform = Platform) %>%
    group_by(UniProt, ChortPlatform, Gene.Symbol) %>%
    summarize(Consensus_LogFC = mean(Effect_Size, na.rm = TRUE), .groups = "drop")

plot_data_bio <- plasmix_mf_dep %>% inner_join(tier1_consensus, by = "UniProt")

cohort_labeller <- c(
    "Olink_3K"    = "UK Biobank\n(Olink 3072)",
    "Olink_HT"    = "Wellness & BAMSE\n(Olink HT)",
    "SomaScan_5K" = "Iceland\n(SomaScan 5K)"
)

plot_data_bio <- plot_data_bio %>%
    mutate(
        ChortPlatform_Label = cohort_labeller[ChortPlatform],
        Platform = factor(Platform, levels = c("SOM", "OLK", "DIA"))
    ) %>%
    filter(!is.na(Platform), !is.na(ChortPlatform_Label))

batch_stats_calc <- plot_data_bio %>%
    group_by(ChortPlatform_Label, Platform, Batch) %>%
    summarize(
        n_prot = n_distinct(UniProt),
        r_p = ifelse(n_prot >= 3, cor(Consensus_LogFC, logFC, method = "pearson", use = "complete.obs"), NA_real_),
        r_s = ifelse(n_prot >= 3, cor(Consensus_LogFC, logFC, method = "spearman", use = "complete.obs"), NA_real_),
        .groups = "drop"
    )

stats_calc <- plot_data_bio %>%
    group_by(ChortPlatform_Label, Platform) %>%
    summarize(n_prot = n_distinct(UniProt), .groups = "drop") %>%
    left_join(
        batch_stats_calc %>%
            group_by(ChortPlatform_Label, Platform) %>%
            summarize(
                r_p_med = median(r_p, na.rm = TRUE),
                r_p_min = min(r_p, na.rm = TRUE),
                r_p_max = max(r_p, na.rm = TRUE),
                r_s_med = median(r_s, na.rm = TRUE),
                r_s_min = min(r_s, na.rm = TRUE),
                r_s_max = max(r_s, na.rm = TRUE),
                .groups = "drop"
            ),
        by = c("ChortPlatform_Label", "Platform")
    )

bio_stats_npc <- stats_calc %>% group_by(ChortPlatform_Label, Platform) %>% reframe(
    item = c("coverage", "pearson", "spearman"), npc_y = c(0.96, 0.88, 0.8),
    label_expr = c(
        sprintf('plain("Coverage:")~%d/"%d"~plain("(%.0f%%)")', n_prot, n_total_tier1, n_prot / n_total_tier1 * 100),
        sprintf('plain("Batch")~italic(r)[med]~"="~%.2f~"("~%.2f*"–"*%.2f~")"', r_p_med, r_p_min, r_p_max),
        sprintf('plain("Batch")~italic(rho)[med]~"="~%.2f~"("~%.2f*"–"*%.2f~")"', r_s_med, r_s_min, r_s_max)
    )
)

integer_breaks <- function(x) {
    breaks <- pretty(x)
    breaks[breaks %% 1 == 0]
}

p5 <- ggplot(plot_data_bio, aes(x = logFC, y = Consensus_LogFC)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "black") +
    geom_point(aes(color = Batch), alpha = 0.8, size = 1) +
    facet_grid2(ChortPlatform_Label ~ Platform, scales = "free", independent = "all") +
    scale_x_continuous(breaks = integer_breaks) +
    scale_y_continuous(breaks = integer_breaks) +
    geom_text_npc(data = bio_stats_npc, aes(npcx = 0.05, npcy = npc_y, label = label_expr), parse = TRUE, hjust = 0, size = 2.8, inherit.aes = FALSE) +
    labs(y = expression("External" ~ beta[sex]), x = expression("Plasmix observed" ~ log[2](M/F)), color = "Batch") +
    scale_color_manual(values = batch_palette) +
    plasmix_theme +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
        axis.line = element_blank(),panel.grid.major = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "none",
        panel.grid = element_blank()
    ) +
    facetted_pos_scales(
        x = lapply(split(plot_data_bio, list(plot_data_bio$Platform, plot_data_bio$ChortPlatform_Label)),
            function(df) {
                r_val <- range(c(df$Consensus_LogFC, df$logFC), na.rm = TRUE)
                scale_x_continuous(limits = r_val, breaks = integer_breaks)
            }),
        y = lapply(split(plot_data_bio, list(plot_data_bio$Platform, plot_data_bio$ChortPlatform_Label)),
            function(df) {
                r_val <- range(c(df$Consensus_LogFC, df$logFC), na.rm = TRUE)
                scale_y_continuous(limits = r_val, breaks = integer_breaks)
            })
    )
print(p5)

tier1_ids <- sex_ext_df %>% filter(Tier == "Tier1") %>% pull(UniProt) %>% unique()
prot_stats_tier1 <- plasmix_mf_dep %>%
    filter(UniProt %in% tier1_ids) %>%
    group_by(UniProt) %>%
    summarize(
        n_platforms = n_distinct(Platform),
        max_abs_logFC = max(abs(logFC), na.rm = TRUE),
        .groups = "drop"
    ) %>%
    filter(n_platforms >= 2)

cat("\n=== [Diagnostic 1] Tier 1 maximum absolute effect quantiles ===\n")
print(quantile(prot_stats_tier1$max_abs_logFC, probs = c(0.1, 0.25, 0.5, 0.75, 0.9, 1), na.rm = TRUE))
cat("\n=== [Diagnostic 2] Targets retained at common effect-size thresholds ===\n")
fc_thresholds <- c(1.1, 1.2, 1.5, 2.0)
print(tibble(
    Threshold_Label = paste0(">=", (fc_thresholds - 1) * 100, "% Change"),
    n_Proteins = sapply(log2(fc_thresholds), function(x) sum(prot_stats_tier1$max_abs_logFC >= x, na.rm = TRUE))
))

target_prots_id <- prot_stats_tier1 %>%
    arrange(desc(max_abs_logFC)) %>%
    slice_head(n = 50) %>%
    pull(UniProt)

consensus_stats <- sex_ext_df %>%
    filter(Tier == "Tier1", UniProt %in% target_prots_id) %>%
    group_by(UniProt, Gene.Symbol) %>%
    summarize(Consensus_LogFC = mean(Effect_Size, na.rm = TRUE), .groups = "drop") %>%
    arrange(Consensus_LogFC)

target_symbols <- consensus_stats$Gene.Symbol
gene_map <- consensus_stats %>% select(UniProt, Gene.Symbol = Gene.Symbol)

analysis_df <- plasmix_mf_dep %>%
    filter(!is.na(Batch), UniProt %in% target_prots_id) %>%
    left_join(gene_map, by = "UniProt")

y_levels <- unique(meta_batch_ht$Batch)
sexm_prots_id <- unique(consensus_stats[consensus_stats$Consensus_LogFC>0,]$UniProt)

plot_heatmap <- analysis_df %>%
    filter(Batch %in% y_levels) %>%
    mutate(Type = ifelse(UniProt %in% sexm_prots_id, "Male high", "Female high")) %>%
    complete(Batch = factor(y_levels, levels = rev(y_levels)), nesting(UniProt, Gene.Symbol, Type)) %>%
    mutate(Batch = factor(Batch, levels = rev(y_levels))) %>%
    mutate(Gene.Symbol = factor(Gene.Symbol, levels = target_symbols))

val_range <- range(plot_heatmap$logFC, na.rm = TRUE)
min_val <- val_range[1]
max_val <- val_range[2]

my_breaks <- c(min_val, -2.0, -1, 0, 1, 2.5, 5, max_val)
my_colors <- c("#193b84", "#6F83B2", "#A8B6CC", "white", "#F3B2A6", "#db5b4a", "#d43a25", "#a91603")
p_ht <- ggplot(plot_heatmap, aes(y = Gene.Symbol, x = Batch, fill = logFC)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_gradientn(
        colors = my_colors,
        values = scales::rescale(my_breaks, from = c(min_val, max_val)),
        limits = c(min_val, max_val),
        breaks = c(-5, -2.5, 0, 2.5, 5, 7.5),
        name = expression(bold(log[2]*FC)),
        na.value = "grey40",
        guide = guide_colorbar(barwidth = unit(4, "cm"), barheight = unit(0.25, "cm"), title.position = "left", title.vjust = 1)) +
    scale_y_discrete(expand = c(0, 0), drop = FALSE) +
    labs(x = NULL, y = NULL) +
    plasmix_theme +
    scale_x_discrete(expand = c(0,0)) +
    theme(
        plot.title = element_text(size = 11, hjust = 0.5, vjust = -1),
        axis.text.x = element_text(size = 7.5, angle = 30, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 7.5),
        plot.margin = margin(5, 5, 5, 5),
        legend.position = "bottom",
        legend.text = element_text(margin = margin(t = 2, r = 0, b = 0, l = 0)),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-8, 0, 2, 0)
    )

showtext_auto()
showtext_opts(dpi = 600)
row1 <- ggarrange(p2, p3, p1, ncol = 3, align = "hv", labels = c("a", "b", "c"),
                  font.label = label_style, label.x = 0, label.y = 1)
row2 <- ggarrange(p5, p_ht, nrow = 1, labels = c("d", "e"), widths = c(2, 1),
                  font.label = label_style, label.x = 0, label.y = 1)
final_assembly <- ggarrange(row1, row2, nrow = 2, heights = c(1, 2.1))
ggsave("figures/ed3_external_sex_consistency.pdf", final_assembly, width = 10, height = 10)
ggsave("figures/ed3_external_sex_consistency.png", final_assembly, width = 10, height = 10, dpi = 600, bg = "white")

# External-cohort validation source data
format_pvalue <- function(x) {
  x <- as.numeric(x)
  out <- character(length(x))
  for (i in seq_along(x)) {
    if (is.na(x[i])) { out[i] <- NA_character_ }
    else if (x[i] < 0.001) { out[i] <- formatC(x[i], format = "e", digits = 2) }
    else { out[i] <- formatC(x[i], format = "f", digits = 3) }
  }
  return(out)
}
st8_wide <- sex_ext_df %>%
    mutate(Tier = factor(Tier, levels = c("Tier1", "Tier2", "Tier3", "Conflict", "NotSig"))) %>%
    select(
        Identifier = UniProt,
        `Gene symbol` = Gene.Symbol,
        Tier = Tier,
        Cohort = Cohort,
        `Effect size` = Effect_Size,
        `P-value` = Pval,
        FDR = FDR
    ) %>%
    pivot_wider(
        names_from = Cohort,
        values_from = c(`Effect size`, `P-value`, FDR),
        names_glue = "{Cohort} {.value}",
        values_fn = list(
            `Effect size` = mean,
            `P-value` = min,
            FDR = min
        )
    ) %>%
    mutate(across(contains("P-value") | contains("FDR"), as.numeric)) %>%
    mutate(across(contains("Effect size"), as.numeric)) %>%
    mutate(across(contains("P-value") | contains("FDR"), format_pvalue)) %>%
    mutate(across(contains("Effect size"), ~ round(., 3))) %>%
    arrange(Tier)
output_list <- list(
    "External_cohort_effects" = st8_wide,
    "Tier_summary" = external_sex_tier_summary,
    "UKB_Iceland" = res_p2$plot_df,
    "Wellness_BAMSE" = res_p3$plot_df,
    "Source_comparison" = plot_data_a,
    "Plasmix_Tier1" = plot_data_bio,
    "Batch_correlations" = batch_stats_calc,
    "Correlation_summary" = stats_calc,
    "Tier1_heatmap" = plot_heatmap
)
write.xlsx(output_list, "tables/SourceData_EDFigure3.xlsx", overwrite = TRUE, keepNA = TRUE, na.string = "NA")
message("Extended Data Figure 3 completed; shared results: results/external_sex_effects_long.tsv.gz and results/external_sex_tier_summary.tsv.gz")
