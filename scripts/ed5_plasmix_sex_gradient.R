# Extended Data Figure 5 | Plasmix and UK Biobank sex gradients

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "readxl", "openxlsx", "cowplot", "patchwork", "ggrepel"), c("clusterProfiler", "org.Hs.eg.db"))
source("utils/figure_style.R")
source("utils/differential_analysis.R")
plasmix_theme <- if (is.function(theme_plasmix)) theme_plasmix() else theme_plasmix
label_style <- list(size = 12, face = "bold")
options(stringsAsFactors = FALSE)
set.seed(2026)

# 1. Inputs ----
required_inputs <- c("results/dea_df_multi.tsv.gz", "results/external_sex_effects_long.tsv.gz", "data/feature_metadata.tsv.gz")
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing Extended Data Figure 5 inputs: ", paste(missing_inputs, collapse = ", "))
dea_df_multi <- fread("results/dea_df_multi.tsv.gz")
sex_ext_df <- fread("results/external_sex_effects_long.tsv.gz")
feat_meta <- fread("data/feature_metadata.tsv.gz") %>% filter(!Is_Protein_Group, !Is_Unknown)
olink_ht_ids <- dea_df_multi %>% filter(Batch %in% c("OLK_P2_B1", "OLK_P2_B2"), ProcessLevel == "NPX") %>% pull(UniqueID) %>% unique()
plasmix_ht_npx <- feat_meta %>% filter(Platform == "OLK", UniqueID %in% olink_ht_ids)

pid_code <- fread(upstream_path("cohort", "profile_ukb", "coding143.tsv")) %>%
    mutate(TargetName = str_trim(str_split(meaning, ";", simplify = TRUE)[, 1]))
pr_data <- fread(upstream_path("cohort", "profile_ukb", "pr_data.txt")) %>% column_to_rownames("eid") %>% t() %>% as.data.frame()
ph_data <- fread(upstream_path("cohort", "profile_ukb", "ph_data.txt")) %>%
    mutate(ethnicity = factor(ethnicity, levels = c("White", "Black", "Asian", "Chinese", "Mixed", "Other")))

# 2. Olink 3072-to-HT mapping ----
canonicalize_uniprot <- function(u_vec) {
    sapply(u_vec, function(x) {
        if (is.na(x) || x == "") return(NA_character_)
        x_clean <- gsub("-[0-9]+", "", x)
        parts <- strsplit(x_clean, "[|_]")[[1]]
        paste(sort(unique(trimws(parts))), collapse = ";")
    }, USE.NAMES = FALSE)
}

raw_3072_df <- read_excel(upstream_path("references", "olink-explore-3072-ht-assay-list.xlsx"), sheet = "3072")
map_3072 <- raw_3072_df %>%
    select(Gene_3072 = `Gene name`, UniProt_3072_Raw = `UniProt ID`) %>%
    distinct() %>%
    mutate(UniProt = canonicalize_uniprot(UniProt_3072_Raw)) %>%
    group_by(UniProt) %>%
    summarize(
        Gene_3072 = paste(unique(Gene_3072), collapse = " ; "),
        UniProt_3072_Raw = paste(unique(UniProt_3072_Raw), collapse = " ; "),
        .groups = "drop"
    )

map_ht <- read_excel(upstream_path("references", "olink-explore-3072-ht-assay-list.xlsx"), sheet = "HT") %>%
    select(Gene_HT = `Gene name`, UniProt_HT_Raw = `UniProt ID`) %>%
    distinct() %>%
    mutate(UniProt = canonicalize_uniprot(UniProt_HT_Raw)) %>%
    group_by(UniProt) %>%
    summarize(
        Gene_HT = paste(unique(Gene_HT), collapse = " ; "),
        UniProt_HT_Raw = paste(unique(UniProt_HT_Raw), collapse = " ; "),
        .groups = "drop"
    )

map_plasmix <- plasmix_ht_npx %>%
    distinct(Gene_Plasmix = TargetName, UniProt_Plasmix_Raw = UniProtID) %>%
    mutate(UniProt = canonicalize_uniprot(UniProt_Plasmix_Raw)) %>%
    group_by(UniProt) %>%
    summarize(
        Gene_Plasmix = paste(unique(Gene_Plasmix), collapse = " ; "),
        UniProt_Plasmix_Raw = paste(unique(UniProt_Plasmix_Raw), collapse = " ; "),
        .groups = "drop"
    )

mapping_ledger <- map_3072 %>%
    full_join(map_ht, by = "UniProt") %>%
    full_join(map_plasmix, by = "UniProt") %>%
    select(
        UniProt,
        Gene_3072, UniProt_3072_Raw,
        Gene_HT, UniProt_HT_Raw,
        Gene_Plasmix, UniProt_Plasmix_Raw
    ) %>%
    arrange(UniProt)

fwrite(mapping_ledger, "results/olink_mapping_ledger.tsv.gz", sep = "\t", na = "NA")
cat(sprintf("Olink mapping ledger: %d normalized UniProt anchors.\n", nrow(mapping_ledger)))

# 3. UK Biobank-to-Plasmix target matching ----
ukb_dict <- pid_code %>% mutate(TargetName_UKB = str_trim(str_split(meaning, ";", simplify = TRUE)[, 1]))
plasmix_targets <- plasmix_ht_npx %>% distinct(TargetName, OlinkID = AssayID)

# Use the mapping ledger as the cross-generation translation layer.
ultimate_match <- ukb_dict %>%
    # Attach the 3072 annotation and retain the normalized UniProt anchor.
    inner_join(
        mapping_ledger %>%
            select(UniProt, Gene_3072, Gene_Plasmix) %>%
            mutate(UniProt_Clean = canonicalize_uniprot(UniProt)) %>%
            separate_rows(Gene_3072, sep = " ; ") %>%
            filter(!is.na(Gene_3072) & Gene_3072 != ""),
        by = c("TargetName_UKB" = "Gene_3072")
    ) %>%
    # Connect the translated Plasmix gene name to the observed HT assay.
    inner_join(plasmix_targets, by = c("Gene_Plasmix" = "TargetName")) %>%
    # Deduplicate at the normalized UniProt anchor.
    # This prevents multi-protein annotations from creating one-to-many target expansion.
    distinct(coding, meaning, TargetName_UKB, UniProt = UniProt_Clean, Gene_Plasmix, OlinkID)

cat(sprintf("Matched targets anchored by normalized UniProt identifiers: %d\n", nrow(ultimate_match)))

# 4. UK Biobank and Plasmix M/F effects ----
# Build the UK Biobank sample metadata required by the limma wrapper.
meta_ukb <- ph_data %>%
    mutate(
        ColName = as.character(eid),
        Sample = ifelse(sex_id == 1, "M", "F")
    ) %>%
    select(ColName, Sample) %>%
    filter(ColName %in% colnames(pr_data)) # Keep only participants represented in the NPX matrix.

# Estimate the UK Biobank male-versus-female effect from participant-level NPX values.
# Rows are Olink coding identifiers, columns are participant IDs, and values are already on the NPX log2 scale.
ukb_dea_raw <- dea_limma_flexible(expr_mat = pr_data, meta_mat = meta_ukb, contrast_pair = "M/F", min_samples_per_group = 3)

# Standardize the UK Biobank result key.
ukb_res <- ukb_dea_raw %>%
    select(coding = UniqueID, UKB_logFC = logFC) %>%
    mutate(coding = as.integer(coding))

# Extract Plasmix Olink HT NPX results from the two selected batches.
# Parse the UniProt accession from the Plasmix feature identifier before matching.
plasmix_res <- dea_df_multi %>%
    filter(Batch %in% c("OLK_P2_B1", "OLK_P2_B2"), ProcessLevel == "NPX", Pair == "M/F") %>%
    mutate(Raw_UniProt = sapply(strsplit(as.character(UniqueID), "_", fixed = TRUE), `[`, 1),
           UniProt = canonicalize_uniprot(Raw_UniProt)) %>%
    group_by(UniProt) %>%
    # Require detection in both batches and complete agreement in M/F effect direction.
    filter(n_distinct(Batch) == 2, all(logFC > 0, na.rm = TRUE) | all(logFC < 0, na.rm = TRUE)) %>%
    summarize(Plasmix_logFC = median(logFC, na.rm = TRUE), .groups = "drop")

# Assemble the matched plotting and enrichment table.
base_df <- ultimate_match %>%
    inner_join(ukb_res, by = "coding") %>%
    inner_join(plasmix_res, by = "UniProt") %>%
    left_join(sex_ext_df %>% distinct(UniProt, Tier), by = "UniProt") %>%
    mutate(
        Tier = ifelse(is.na(Tier), "NotSig", as.character(Tier)),
        # Operational excess-shift rule: UK Biobank change <= 5% and Plasmix change >= twofold.
        is_excess_shift = (abs(UKB_logFC) <= log2(1.05)) & (abs(Plasmix_logFC) >= log2(2))
    )
cat(sprintf("Matched Olink assay pairs retained for analysis: %d\n", nrow(base_df)))
cat(sprintf("Excess-shift proteins: %d\n", sum(base_df$is_excess_shift)))

# 5. Panel a: UK Biobank versus Plasmix effects ----
# Calculate global and Tier 1 Pearson and Spearman correlations.
n_overall <- nrow(base_df)
cor_p_overall <- cor.test(base_df$UKB_logFC, base_df$Plasmix_logFC, method = "pearson")
cor_s_overall <- cor.test(base_df$UKB_logFC, base_df$Plasmix_logFC, method = "spearman")

tier1_df <- base_df %>% filter(Tier == "Tier1")
n_tier1 <- nrow(tier1_df)
cor_p_tier1 <- cor.test(tier1_df$UKB_logFC, tier1_df$Plasmix_logFC, method = "pearson")
cor_s_tier1 <- cor.test(tier1_df$UKB_logFC, tier1_df$Plasmix_logFC, method = "spearman")

# Build separate parsed annotation lines for the global and Tier 1 summaries.
cor_label_overall <- sprintf(
    "plain('Global')~'(' * italic(n) == %d * '):'~italic(r) == '%.2f'*','~italic(rho) == '%.2f'",
    n_overall, cor_p_overall$estimate, cor_s_overall$estimate
)
cor_label_tier1 <- sprintf(
    "plain('Tier 1')~'(' * italic(n) == %d * '):'~italic(r) == '%.2f'*','~italic(rho) == '%.2f'",
    n_tier1, cor_p_tier1$estimate, cor_s_tier1$estimate
)

# Order plotting categories so highlighted proteins are drawn above the background.
plot_data_a <- base_df %>%
    mutate(Plot_Category = case_when(is_excess_shift ~ "Excess-Shift", Tier == "Tier1" ~ "Tier1", TRUE ~ "Other")) %>%
    arrange(match(Plot_Category, c("Other", "Tier1", "Excess-Shift")))

# Label the three Tier 1 proteins with the largest absolute Plasmix effects.
label_data_a <- plot_data_a %>%
    filter(Plot_Category %in% c("Tier1")) %>%
    group_by(Plot_Category) %>% arrange(desc(abs(Plasmix_logFC))) %>% slice_head(n = 3) %>% ungroup()

lim <- max(abs(range(plot_data_a$Plasmix_logFC, plot_data_a$UKB_logFC)))

# Preserve the established panel-a plotting layers.
p_a <- ggplot(plot_data_a, aes(x = UKB_logFC, y = Plasmix_logFC)) +
    # Reference lines mark the operational excess-shift thresholds and the identity relationship.
    geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = c(-log2(1.05), log2(1.05)), linetype = "dotted", color = "grey60", linewidth = 0.3) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", alpha = 0.5, linewidth = 0.3) +
    # Fit the global linear relationship without altering the highlighted subsets.
    geom_smooth(method = "lm", color = "#D85170", se = TRUE, alpha = 0.5, fullrange = TRUE, linewidth = 0.3) +
    # Draw background, Tier 1 and excess-shift proteins as separate layers.
    geom_point(data = filter(plot_data_a, Plot_Category == "Other"), color = "grey80", alpha = 0.5, size = 0.8) +
    geom_point(data = filter(plot_data_a, Plot_Category == "Tier1"), color = "#D85170", alpha = 0.8, size = 1) +
    geom_point(data = filter(plot_data_a, Plot_Category == "Excess-Shift"), shape = 4, color = "#4d8cd5", size = 0.8) +
    # Add labels only for the selected Tier 1 proteins.
    geom_text_repel(
        data = label_data_a,
        aes(label = Gene_Plasmix, color = Plot_Category),
        size = 2.5, fontface = "bold",
        max.overlaps = Inf,
        show.legend = FALSE,
        seed = 42,
        nudge_x = 1.5, nudge_y = -0.5   # Preserve the established label displacement.
    ) +
    scale_color_manual(values = c("Tier1" = "#D85170", "Excess-Shift" = "#4d8cd5")) +
    # Use matched log2(M/F) labels on both axes.
    labs(x = expression(UKBiobank~log[2](M/F)), y = expression(Plasmix~log[2](M/F))) +
    # Add the two correlation summaries as independently aligned annotations.
    annotate("text", x = lim * 0.05, y = -Inf, label = cor_label_overall, parse = TRUE, hjust = 0, vjust = -2.2, size = 2.5, color = "black") +
    annotate("text", x = lim * 0.05, y = -Inf, label = cor_label_tier1, parse = TRUE, hjust = 0, vjust = -1, size = 2.5, color = "black") +
    plasmix_theme +
    scale_x_continuous(n.breaks = 6) + scale_y_continuous(n.breaks = 4) +
    theme(panel.grid.major = element_blank(), panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35), axis.line = element_blank())

# 6. Panel b: UK Biobank resampling ----
# The observed UK Biobank effect uses the full available male and female participant groups.
# Direct resampling repeatedly estimates mean male NPX minus mean female NPX from 54 male and 51 female participants.
# Pseudo-pool resampling averages linearized NPX values before returning to the log2 scale.
set.seed(2026)
n_iter <- 5000
# Identify the available male and female participant columns.
male_eids <- meta_ukb %>% filter(Sample == "M") %>% pull(ColName)
female_eids <- meta_ukb %>% filter(Sample == "F") %>% pull(ColName)
# Use exactly the same matched target set for observed and resampled effect-range widths.
valid_codings <- as.character(base_df$coding)
# Slice the participant matrices once before iteration to reduce repeated indexing overhead.
pr_m_mat <- pr_data[valid_codings, as.character(male_eids)]
pr_f_mat <- pr_data[valid_codings, as.character(female_eids)]
resample_stats <- vector("list", n_iter)
for (i in 1:n_iter) {
    # Sample participant columns without replacement.
    s_m <- sample(ncol(pr_m_mat), 54, replace = FALSE)
    s_f <- sample(ncol(pr_f_mat), 51, replace = FALSE)
    mat_m <- pr_m_mat[, s_m]
    mat_f <- pr_f_mat[, s_f]
    # A. Direct logFC
    fc_direct <- rowMeans(mat_m, na.rm = TRUE) - rowMeans(mat_f, na.rm = TRUE)
    # Pseudo-pool effect: linearize NPX, average, and transform back to log2.
    pool_m <- log2(rowMeans(2^mat_m, na.rm = TRUE))
    pool_f <- log2(rowMeans(2^mat_f, na.rm = TRUE))
    fc_pool <- pool_m - pool_f
    resample_stats[[i]] <- tibble(
        Iteration = i,
        Direct_Width = quantile(fc_direct, 0.975, na.rm = TRUE) - quantile(fc_direct, 0.025, na.rm = TRUE),
        Pool_Width = quantile(fc_pool, 0.975, na.rm = TRUE) - quantile(fc_pool, 0.025, na.rm = TRUE)
    )
}
resample_df <- bind_rows(resample_stats)

# Calculate observed central 95% widths on the same matched target set.
obs_plasmix_width <- quantile(base_df$Plasmix_logFC, 0.975, na.rm=TRUE) - quantile(base_df$Plasmix_logFC, 0.025, na.rm=TRUE)
obs_ukb_width <- quantile(base_df$UKB_logFC, 0.975, na.rm=TRUE) - quantile(base_df$UKB_logFC, 0.025, na.rm=TRUE)

# Reshape the simulation output for plotting.
sim_plot_df <- resample_df %>%
    pivot_longer(c(Direct_Width, Pool_Width), names_to = "Method", values_to = "Central_95_Width") %>%
    mutate(Method = factor(Method, levels = c("Direct_Width", "Pool_Width")))

# Overlay the direct and pseudo-pool resampling distributions.
p_b <- ggplot(sim_plot_df, aes(x = Central_95_Width, fill = Method)) +
    geom_histogram(alpha = 0.9, position = "identity", bins = 100, color = "white", linewidth = 0.1) +
    geom_vline(xintercept = obs_ukb_width, color = "#065EAD", linetype = "dashed", linewidth = 0.5) +
    geom_vline(xintercept = obs_plasmix_width, color = "#56106E", linetype = "dashed", linewidth = 0.5) +
    annotate("text", x = obs_ukb_width, y = Inf, label = "UKB observed", vjust = 1.5, hjust = 1, angle = 90, color = "#065EAD", size = 2.5, fontface = "bold") +
    annotate("text", x = obs_plasmix_width, y = Inf, label = "Plasmix observed", vjust = -0.5, hjust = 1, angle = 90, color = "#56106E", size = 2.5, fontface = "bold") +
    scale_fill_manual(
        values = c("Direct_Width" = "#a8ddb5", "Pool_Width" = "#fdbb84"),
        labels = c("Direct_Width" = "Direct resampling (54M/51F)", "Pool_Width" = "Pseudo-pool resampling (54M/51F)")
    ) +
    labs(x = expression("Central 95% interval width ("*log[2]~"scale)"), y = "Frequency (5000 iterations)") +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.01))) +
    plasmix_theme +
    theme(legend.position = c(0.9, 1), legend.justification = c(1, 1), legend.background = element_blank(), panel.grid.major = element_blank())

# 7. Panel c: GO enrichment of excess-shift proteins ----
# Define the matched-target universe and excess-shift foreground.
bg_uniprot <- base_df %>% pull(UniProt) %>% unique()
excess_uniprot <- base_df %>% filter(is_excess_shift) %>% pull(UniProt) %>% unique()

# Map foreground and background UniProt identifiers to Entrez identifiers using the human annotation database.
id_map_ext <- clusterProfiler::bitr(
    unique(c(bg_uniprot, excess_uniprot)),
    fromType = "UNIPROT",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db::org.Hs.eg.db
)
bg_entrez <- id_map_ext %>% filter(UNIPROT %in% bg_uniprot) %>% pull(ENTREZID) %>% unique()
excess_entrez <- id_map_ext %>% filter(UNIPROT %in% excess_uniprot) %>% pull(ENTREZID) %>% unique()

cat("GO foreground proteins (Excess-shift):", length(excess_entrez), "\n")
cat("GO background proteins (Valid Matches):", length(bg_entrez), "\n")

# Run BP, CC and MF enrichment and collapse semantically redundant terms.
run_go_ext <- function(ont_use) {
    # Test over-representation against the matched-target universe.
    res <- clusterProfiler::enrichGO(
        gene = excess_entrez,
        universe = bg_entrez,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        ont = ont_use,
        pvalueCutoff = 0.05
    )
    if (is.null(res) || nrow(res) == 0) return(NULL)
    # Retain the most significant representative among semantically similar GO terms.
    res_sim <- clusterProfiler::simplify(res, cutoff = 0.65, by = "p.adjust", select_fun = min)
    as.data.frame(res_sim) %>% mutate(Ontology = ont_use)
}

enrich_ext_all <- map_dfr(c("BP", "CC", "MF"), run_go_ext)
enrich_ext_summary <- enrich_ext_all %>%
    mutate(Description = paste0(toupper(substr(Description, 1, 1)), substr(Description, 2, nchar(Description)))) %>%
    mutate(
        GR_val = as.numeric(sub("/.*", "", GeneRatio)) / as.numeric(sub(".*/", "", GeneRatio)),
        BR_val = as.numeric(sub("/.*", "", BgRatio)) / as.numeric(sub(".*/", "", BgRatio)),
        FoldEnrichment = GR_val / BR_val
    ) %>%
    group_by(Ontology) %>%
    arrange(p.adjust) %>%
    mutate(rank = row_number()) %>%
    ungroup()

ontology_labs_ext <- c(BP = "Biological process", CC = "Cellular component", MF = "Molecular function")

# Plot fold enrichment with protein count and adjusted P value.
p_c <- ggplot(enrich_ext_summary, aes(x = FoldEnrichment, y = reorder(Description, FoldEnrichment))) +
    geom_segment(aes(x = 0, xend = FoldEnrichment, y = Description, yend = Description), color = "gray90", linewidth = 0.8) +
    geom_point(aes(size = Count, fill = -log10(p.adjust)), alpha = 0.8, color = "black", shape = 21) +
    scale_fill_viridis_c(option = "rocket", name = expression(bold(-log[10](P[adj]))), direction = -1, end = 0.8, breaks = c(1.5, 2.0, 2.5)) +
    facet_wrap(. ~ Ontology, scales = "free", space = "free_y", nrow = 2, labeller = labeller(Ontology = as_labeller(ontology_labs_ext))) +
    scale_size_continuous(range = c(1, 5), name = "Protein count", breaks = c(5, 10, 20, 40)) +
    scale_y_discrete(labels = function(x) str_wrap(x, width = 70)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1)), limits = c(0, NA)) +
    labs(x = "Fold enrichment", y = NULL) +
    plasmix_theme +
    theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.35),
        axis.line.x = element_blank(), axis.line.y = element_blank(),
        strip.background = element_rect(fill = "gray20"),
        strip.text = element_text(color = "white", face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.key.size = unit(0.4, "lines"),
        legend.text = element_text(margin = margin(t = 2, r = 0, b = 0, l = 0), vjust = 0.5),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(0, 0, 5, -100),
        legend.title = element_text(face = "bold", vjust = 0.5)
    )

# 8. Assemble and export ----
final_fig <- ggpubr::ggarrange(p_a, p_b, p_c, nrow = 1, widths = c(0.9, 1, 1.2), labels = c("a", "b", "c"), font.label = label_style, label.x = 0, label.y = 1, hjust = -0.2, vjust = 1.2)
ggsave("figures/ed5_plasmix_sex_gradient.pdf", final_fig, width = 10, height = 3)
ggsave("figures/ed5_plasmix_sex_gradient.png", final_fig, width = 10, height = 3, dpi = 600, bg = "white")

source_data <- list(
    "Matched_UKB_Plasmix" = base_df,
    "UKB_resampling" = resample_df,
    "Excess_shift_GO" = enrich_ext_summary
)
write.xlsx(source_data, "tables/SourceData_EDFigure5.xlsx", overwrite = TRUE, keepNA = TRUE, na.string = "NA")
fwrite(enrich_ext_summary, "results/plasmix_ukb_enrichment.tsv.gz", sep = "\t", na = "NA")
message("Extended Data Figure 5 completed: figures/ed5_plasmix_sex_gradient.pdf；source data: tables/SourceData_EDFigure5.xlsx")
