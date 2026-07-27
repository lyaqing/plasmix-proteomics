# Prepare protein physicochemical annotations

# 0. Setup ----
source("scripts/_project_setup.R")
use_packages(c("data.table", "tidyverse", "bio3d", "pbapply", "Peptides", "reticulate"))
set.seed(2026)
options(stringsAsFactors = FALSE)
# 1. Project proteins and UniProt annotations ----
profile_columns <- c("UniProtID")
long_df <- fread("data/protein_profiles_long.tsv.gz", select = profile_columns)
target_uniprot_ids <- unique(long_df$UniProtID)
target_uniprot_ids <- target_uniprot_ids[
    !is.na(target_uniprot_ids) & target_uniprot_ids != ""
]

uniprot_raw <- read_tsv(upstream_path("references/uniprotkb_AND_model_organism_9606_2026_04_01.tsv.gz"), show_col_types = FALSE) %>%
    rename(GeneName = `Gene Names (primary)`) %>%
    mutate(across(everything(), ~replace_na(as.character(.), ""))) %>%
    filter(Entry %in% target_uniprot_ids)

nonstd_seq_qc <- uniprot_raw %>%
    mutate(Clean_Seq = str_replace_all(Sequence, "\\s+", "")) %>%
    summarize(
        N_total = n(), N_nonstd = sum(str_detect(Clean_Seq, "[^ACDEFGHIKLMNPQRSTVWY]"), na.rm = TRUE),
        Percent_nonstd = 100 * N_nonstd / N_total
    )
print(nonstd_seq_qc)

nonstd_seq_table <- uniprot_raw %>%
    mutate(Clean_Seq = str_replace_all(Sequence, "\\s+", "")) %>%
    filter(str_detect(Clean_Seq, "[^ACDEFGHIKLMNPQRSTVWY]")) %>%
    transmute(
        Entry,
        GeneName,
        ProteinName = `Protein names`,
        Length = nchar(Clean_Seq),
        Nonstandard_symbols = str_extract_all(
            Clean_Seq,
            "[^ACDEFGHIKLMNPQRSTVWY]"
        ) %>%
            map_chr(~paste(sort(unique(.x)), collapse = "; "))
    )
print(nonstd_seq_table)

# Recorded notebook output (reference only; regenerate by running this script):
# N_total = 12,008; N_nonstd = 16; Percent_nonstd = 0.133%.

# 2. Defensive annotation and sequence helpers ----

# Peptides-derived indices are set to NA for sequences containing non-standard amino acids.
has_nonstd_aa <- function(seq) {
  is.na(seq) || nchar(seq) == 0 || str_detect(seq, "[^ACDEFGHIKLMNPQRSTVWY]")
}
safe_charge <- function(seq) {
  if (has_nonstd_aa(seq)) return(NA_real_)
  tryCatch(charge(seq, pH = 7.4, pKscale = "Lehninger"), error = function(e) NA_real_)
}
safe_pI <- function(seq) {
  if (has_nonstd_aa(seq)) return(NA_real_)
  tryCatch(pI(seq, pKscale = "Lehninger"), error = function(e) NA_real_)
}
safe_instaIndex <- function(seq) {
  if (has_nonstd_aa(seq) || nchar(seq) < 2) return(NA_real_)
  tryCatch(instaIndex(seq), error = function(e) NA_real_)
}
safe_aIndex <- function(seq) {
  if (has_nonstd_aa(seq)) return(NA_real_)
  tryCatch(aIndex(seq), error = function(e) NA_real_)
}

is_blank <- function(x) { is.na(x) || str_trim(x) == "" }

# Parse UniProt feature intervals and merge overlaps before calculating sequence coverage.
parse_feature_ranges <- function(x, key) {
  m <- str_match_all(x, paste0("\\b", key, "\\b\\s+(\\d+)\\.\\.(\\d+)"))[[1]]
  if (nrow(m) == 0) {
    return(matrix(numeric(0), ncol = 2, dimnames = list(NULL, c("start", "end"))))
  }
  ranges <- cbind(start = as.integer(m[, 2]), end = as.integer(m[, 3]))
  ranges <- ranges[!is.na(ranges[, "start"]) & !is.na(ranges[, "end"]) & ranges[, "end"] >= ranges[, "start"], , drop = FALSE]
  ranges
}

calc_union_length <- function(ranges) {
  if (nrow(ranges) == 0) return(0L)
  ranges <- ranges[order(ranges[, "start"], ranges[, "end"]), , drop = FALSE]
  total <- 0L
  s <- ranges[1, "start"]
  e <- ranges[1, "end"]
  if (nrow(ranges) > 1) {
    for (i in 2:nrow(ranges)) {
      if (ranges[i, "start"] <= e + 1L) {
        e <- max(e, ranges[i, "end"])
      } else {
        total <- total + e - s + 1L
        s <- ranges[i, "start"]
        e <- ranges[i, "end"]
      }
    }
  }
  total + e - s + 1L
}

calc_feature_fraction <- function(x, key, length_aa, empty_value) {
  if (is_blank(x)) return(empty_value)
  if (is.na(length_aa) || length_aa <= 0) return(NA_real_)
  ranges <- parse_feature_ranges(x, key)
  if (nrow(ranges) == 0) return(NA_real_)
  calc_union_length(ranges) / length_aa
}

calc_feature_count <- function(x, key) {
  if (is_blank(x)) return(0L)
  ranges <- parse_feature_ranges(x, key)
  if (nrow(ranges) == 0) return(NA_integer_)
  nrow(ranges)
}

count_feature_records <- function(x, key) {
  if (is_blank(x)) return(0L)
  str_count(x, paste0("\\b", key, "\\b"))
}

extract_feature_notes <- function(x) {
  if (is_blank(x)) return(NA_character_)
  m <- str_match_all(x, '/note="([^"]+)"')[[1]]
  if (nrow(m) == 0) return(NA_character_)
  paste(sort(unique(m[, 2])), collapse = "; ")
}

count_secretory_glyco <- function(x) {
  if (is_blank(x)) return(0L)
  records <- str_split(x, "(?=\\bCARBOHYD\\b)")[[1]] |> str_trim()
  records <- records[str_detect(records, "^CARBOHYD")]
  if (length(records) == 0) return(0L)
  records_l <- str_to_lower(records)
  sum(
    !str_detect(records_l, "glycation") &
      !str_detect(records_l, "in vitro") &
      !str_detect(records_l, "in variant")
  )
}

detect_subunit_assembly <- function(x) {
  if (is_blank(x)) return(0L)
  txt <- str_to_lower(x)
  has_oligomer <- str_detect(txt, "\\b(?:homo|hetero)?(?:dimer|trimer|tetramer|pentamer|hexamer|heptamer|octamer|decamer)(?:s|ic|ize|izes|ization)?\\b|\\b(?:homo|hetero|hom|heter)?oligomer(?:s|ic|ize|izes|ization)?\\b|\\bmultimer(?:s|ic|ize|izes|ization)?\\b")
  has_complex <- str_detect(txt, "\\b(?:sub)?complex(?:es)?\\b")
  has_composed_assembly <- str_detect(txt, "\\b(?:composed of|consists of)\\b.*\\b(?:chains?|subunits?|units?|molecules?)\\b")
  as.integer(has_oligomer | has_complex | has_composed_assembly)
}

count_named_isoforms <- function(x) {
  if (is_blank(x)) return(1L)
  if (!str_detect(x, "Named isoforms=")) return(1L)
  as.integer(str_extract(x, "(?<=Named isoforms=)\\d+"))
}

count_interactors <- function(x) {
  if (is_blank(x)) return(0L)
  str_count(x, ";") + 1L
}

# 3. Sequence and UniProt-derived features ----

uniprot_clean <- uniprot_raw %>%
  rowwise() %>%
  mutate(
    Length_AA = as.numeric(Length),
    Mass_Da = as.numeric(str_remove_all(Mass, ",")),
    Glyco_Count = count_secretory_glyco(Glycosylation),
    Glyco_Density = Glyco_Count / Length_AA,
    Disulfide_Count = count_feature_records(`Disulfide bond`, "DISULFID"),
    Lipidation_Count = count_feature_records(Lipidation, "LIPID"),
    Has_Lipidation = as.integer(Lipidation_Count > 0),
    Lipidation_Notes = extract_feature_notes(Lipidation),
    Propeptide_Count = count_feature_records(Propeptide, "PROPEP"),
    Has_Propeptide = as.integer(Propeptide_Count > 0),
    Propeptide_Fraction = calc_feature_fraction(Propeptide, "PROPEP", Length_AA, empty_value = 0),
    Propeptide_Notes = extract_feature_notes(Propeptide),
    Transmembrane_Count = count_feature_records(Transmembrane, "TRANSMEM"),
    Is_Transmembrane = as.integer(Transmembrane_Count > 0),
    Is_Secreted = {
      loc <- str_to_lower(replace_na(`Subcellular location [CC]`, ""))
      as.integer(str_detect(loc, "\\bsecreted\\b|extracellular space"))
    },
    Is_Membrane = {
      loc <- str_to_lower(replace_na(`Subcellular location [CC]`, ""))
      loc_no_basement <- str_replace_all(loc, "basement membrane", "")
      as.integer(str_detect(loc_no_basement, "\\bmembrane\\b|cell surface|lipid-anchor|gpi-anchor") | Transmembrane_Count > 0)
    },
    Is_Complex = detect_subunit_assembly(`Subunit structure`),
    Isoform_Count = count_named_isoforms(`Alternative products (isoforms)`),
    Hub_Score = count_interactors(`Interacts with`),
    Helix_Fraction = calc_feature_fraction(Helix, "HELIX", Length_AA, empty_value = NA_real_),
    Beta_Fraction = calc_feature_fraction(`Beta strand`, "STRAND", Length_AA, empty_value = NA_real_),
    CompBias_Fraction = calc_feature_fraction(`Compositional bias`, "COMPBIAS", Length_AA, empty_value = 0),
    CompBias_Count = calc_feature_count(`Compositional bias`, "COMPBIAS"),
    CompBias_Notes = extract_feature_notes(`Compositional bias`)
  ) %>%
  ungroup() %>%
  mutate(Clean_Seq = str_replace_all(Sequence, "\\s+", "")) %>%
  rowwise() %>%
  mutate(
      Net_Charge_7.4 = safe_charge(Clean_Seq), pI_Value = safe_pI(Clean_Seq),
      Delta_pI_7.4 = ifelse(is.na(pI_Value), NA, abs(pI_Value - 7.4)),
      Sequence_Instability = safe_instaIndex(Clean_Seq), Aliphatic_Index = safe_aIndex(Clean_Seq)
  ) %>%
  ungroup() %>%
  select(
      Entry, GeneName, Length_AA, Mass_Da, Glyco_Count, Glyco_Density, Disulfide_Count, Has_Lipidation,
      Lipidation_Count, Lipidation_Notes, Has_Propeptide, Propeptide_Count, Propeptide_Fraction,
      Propeptide_Notes, Transmembrane_Count, Is_Transmembrane, Is_Secreted, Is_Membrane, Is_Complex,
      Isoform_Count, Hub_Score, Helix_Fraction, Beta_Fraction, CompBias_Fraction, CompBias_Count,
      CompBias_Notes, Net_Charge_7.4, pI_Value, Delta_pI_7.4, Sequence_Instability, Aliphatic_Index
  )

# 4. Integrate HPA abundance and annotation fields ----
# The HPA-derived fields are retained in the shared `physchem_matrix` so that Figure 1, its HPA coverage Extended Data figure, and later abundance-related analyses can reuse the same frozen annotation table without reading the original HPA release again.
# The following fields are descriptive annotations rather than correlation/model inputs:
# - `BloodConc_log10_pgml`
# - `HPA_Protein_Class`
# - `HPA_Subcellular`
# `Log10_Abundance` remains the model-facing abundance feature recorded in the physicochemical dictionary.

# Aggregate all HPA records per UniProt before applying the secreted/membrane keyword rules.
hpa_raw <- fread(upstream_path("references/proteinatlas_v25.0.tsv.gz"))
hpa_clean <- hpa_raw %>%
    group_by(Uniprot) %>%
    summarize(
        Conc_MS = median(
            `Blood concentration - Conc. blood MS [pg/L]`,
            na.rm = TRUE
        ),
        Conc_IM = median(
            `Blood concentration - Conc. blood IM [pg/L]`,
            na.rm = TRUE
        ),
        HPA_Protein_Class = paste(
            sort(unique(
                `Protein class`[
                    !is.na(`Protein class`) & `Protein class` != ""
                ]
            )),
            collapse = "; "
        ),
        .groups = "drop"
    ) %>%
    mutate(
        Baseline_Abundance_pgL = coalesce(Conc_MS, Conc_IM),
        Log10_Abundance = ifelse(
            is.finite(Baseline_Abundance_pgL) &
                Baseline_Abundance_pgL > 0,
            log10(Baseline_Abundance_pgL),
            NA_real_
        ),
        BloodConc_log10_pgml = ifelse(
            is.finite(Baseline_Abundance_pgL) &
                Baseline_Abundance_pgL > 0,
            log10(Baseline_Abundance_pgL / 1000),
            NA_real_
        ),
        Abundance_Source = case_when(
            is.finite(Conc_MS) ~ "MS",
            !is.finite(Conc_MS) & is.finite(Conc_IM) ~ "IM",
            TRUE ~ "Unknown"
        ),
        HPA_Subcellular = case_when(
            str_detect(
                HPA_Protein_Class,
                regex("membrane", ignore_case = TRUE)
            ) &
                str_detect(
                    HPA_Protein_Class,
                    regex("secreted", ignore_case = TRUE)
                ) ~ "Secreted & Membrane",
            str_detect(
                HPA_Protein_Class,
                regex("membrane", ignore_case = TRUE)
            ) ~ "Membrane",
            str_detect(
                HPA_Protein_Class,
                regex("secreted", ignore_case = TRUE)
            ) ~ "Secreted",
            HPA_Protein_Class == "" ~ "Unknown",
            TRUE ~ "Intracellular"
        )
    ) %>%
    select(Uniprot, Log10_Abundance, BloodConc_log10_pgml, Abundance_Source, HPA_Protein_Class, HPA_Subcellular)

anno_matrix_v1 <- uniprot_clean %>%
    left_join(hpa_clean, by = c("Entry" = "Uniprot"))
target_entries <- anno_matrix_v1$Entry

# 5. AlphaFold structural features ----
# This section requires the local AlphaFold PDB archive and the Python module `freesasa` in the environment used by `reticulate`. The notebook checks the dependency but does not install software automatically.

# AlphaFold coordinates provide compactness, confidence and solvent-accessibility features.
if (!py_module_available("freesasa")) {
    stop("Python module 'freesasa' is required. Install it in the reticulate environment before running this section.")
}
fs <- import("freesasa")
af_features_list <- pblapply(target_entries, function(entry) {
  pdb_file <- upstream_path(
      "references", "UP000005640_9606_HUMAN_v6", paste0("AF-", entry, "-F1-model_v6.pdb.gz")
  )
  if (file.exists(pdb_file)) {
    pdb <- suppressWarnings(read.pdb(pdb_file))
    b_factors <- pdb$atom$b[pdb$calpha]
    plddt_low_frac <- sum(b_factors < 50) / length(b_factors)
    rg_value <- rgyr(pdb)
    tmp_pdb <- tempfile(fileext = ".pdb")
    con <- gzfile(pdb_file, "rt")
    pdb_lines <- readLines(con, warn = FALSE)
    close(con)
    writeLines(pdb_lines, tmp_pdb)
    structure <- fs$Structure(tmp_pdb)
    result <- fs$calc(structure)
    SASA_Total <- result$totalArea()
    sasa_classes <- fs$classifyResults(result, structure)
    SASA_Hydrophobic <- sasa_classes$Apolar
    if (is.null(SASA_Hydrophobic)) {
      SASA_Hydrophobic <- sasa_classes$apolar
    }
    SASA_Hydro_Ratio <- ifelse(SASA_Total > 0, SASA_Hydrophobic / SASA_Total, NA)
    unlink(tmp_pdb)
    return(data.frame(
      Entry = entry,
      pLDDT_Fraction_Low = plddt_low_frac,
      Radius_of_Gyration = rg_value,
      SASA_Total = SASA_Total,
      SASA_Hydrophobic = SASA_Hydrophobic,
      SASA_Hydro_Ratio = SASA_Hydro_Ratio,
      stringsAsFactors = FALSE
    ))
  } else {
    return(NULL)
  }
})
af_features <- bind_rows(af_features_list)

# 6. Assemble and export the annotation matrix ----

physchem_matrix <- anno_matrix_v1 %>% left_join(af_features, by = "Entry")

fwrite(as.data.table(physchem_matrix), "data/physchem_matrix.tsv.gz", sep = "\t", na = "NA")

# 7. Feature correlation, decisions and dictionary ----

annotation_only_columns <- c("BloodConc_log10_pgml", "HPA_Protein_Class", "HPA_Subcellular")

features_for_cor <- physchem_matrix %>%
    select(where(is.numeric)) %>%
    select(-any_of(c("CompBias_Count", "Transmembrane_Count", "BloodConc_log10_pgml")))

stopifnot(
    all(vapply(features_for_cor, is.numeric, logical(1))),
    !any(annotation_only_columns %in% names(features_for_cor))
)

cor_matrix <- cor(features_for_cor, method = "spearman", use = "pairwise.complete.obs")

spearman_p_matrix <- function(df) {
    mat <- as.data.frame(df)
    out <- matrix(NA_real_, ncol(mat), ncol(mat), dimnames = list(colnames(mat), colnames(mat)))
    diag(out) <- 0
    if (ncol(mat) < 2) return(out)
    for (i in seq_len(ncol(mat) - 1)) {
        for (j in (i + 1):ncol(mat)) {
            ok <- is.finite(mat[[i]]) & is.finite(mat[[j]])
            if (sum(ok) >= 3) {
                p <- suppressWarnings(cor.test(mat[[i]][ok], mat[[j]][ok], method = "spearman", exact = FALSE)$p.value)
                out[i, j] <- p
                out[j, i] <- p
            }
        }
    }
    out
}
cor_pmat <- spearman_p_matrix(features_for_cor)
fwrite(
    as.data.table(as.data.frame(cor_matrix), keep.rownames = "Feature"), "results/physchem_cor_matrix.tsv",
    sep = " ", na = "NA"
)
fwrite(
    as.data.table(as.data.frame(cor_pmat), keep.rownames = "Feature"), "results/physchem_cor_pmat.tsv",
    sep = " ", na = "NA"
)
feature_decision <- tribble(
  ~Category, ~Property, ~Feature, ~Definition, ~Computation, ~Representation, ~Retained, ~Redundant_With, ~Decision_Rationale,

  "Structure", "Sequence length", "Length_AA",
  "Length of the UniProt canonical protein sequence.",
  "Parsed directly from the UniProt Length field.",
  "continuous", "No", "Mass_Da",
  "Excluded because sequence length was nearly collinear with molecular mass; molecular mass was retained as the more directly physical descriptor of protein size.",

  "Structure", "Molecular mass", "Mass_Da",
  "Molecular mass of the UniProt canonical protein sequence.",
  "Parsed from the UniProt Mass field after removing comma separators.",
  "continuous", "Yes", NA_character_,
  "Retained as the primary descriptor of absolute protein size.",

  "Structure", "Radius of gyration", "Radius_of_Gyration",
  "Spatial radius of gyration of the AlphaFold-predicted protein structure.",
  "Calculated from AlphaFold PDB coordinates using bio3d::rgyr.",
  "continuous", "Yes", NA_character_,
  "Retained to represent spatial extent and compactness, which is not equivalent to molecular mass.",

  "Structure", "Subunit assembly", "Is_Complex",
  "Whether UniProt Subunit structure annotations indicate subunit assembly, oligomerization or complex formation.",
  "Set to 1 when the Subunit structure text contains oligomeric, multimeric, complex, or composed-of-chain/subunit assembly evidence; otherwise set to 0.",
  "binary", "Yes", NA_character_,
  "Retained as a marker of complex or assembly state that may affect target definition and platform-specific detection.",

  "Structure", "Isoform count", "Isoform_Count",
  "Number of named UniProt isoforms, reflecting isoform complexity and potential target-identity ambiguity.",
  "Extracted from the UniProt Alternative products field using Named isoforms=N; proteins without a named-isoform record were assigned 1.",
  "count", "Yes", NA_character_,
  "Retained because isoform complexity may affect whether different platforms measure identical molecular entities, epitopes, aptamer-binding regions or peptides.",

  "Structure", "Interaction partners", "Hub_Score",
  "Number of UniProt-curated interaction partners listed for the protein.",
  "Calculated as the number of semicolon-delimited entries in the UniProt Interacts with field; empty fields were assigned 0.",
  "count", "Yes", NA_character_,
  "Retained as a proxy for documented interaction-partner burden, which may affect complex state, epitope accessibility, indirect capture or extraction behavior.",

  "Structure", "Alpha-helix fraction", "Helix_Fraction",
  "Fraction of the protein sequence covered by UniProt HELIX annotations.",
  "Calculated as the union length of HELIX coordinate intervals divided by full-length protein length; empty HELIX fields were treated as missing because absence of annotation does not imply absence of helices.",
  "fraction", "Yes", NA_character_,
  "Retained as a secondary-structure coverage feature.",

  "Structure", "Beta-strand fraction", "Beta_Fraction",
  "Fraction of the protein sequence covered by UniProt STRAND annotations.",
  "Calculated as the union length of STRAND coordinate intervals divided by full-length protein length; empty STRAND fields were treated as missing because absence of annotation does not imply absence of beta strands.",
  "fraction", "Yes", NA_character_,
  "Retained as a secondary-structure coverage feature.",

  "Structure", "Aliphatic index", "Aliphatic_Index",
  "Sequence-based aliphatic index, reflecting the relative volume of aliphatic side chains.",
  "Calculated from the full-length UniProt sequence using Peptides::aIndex; sequences containing non-standard amino-acid symbols were set to missing.",
  "continuous", "Yes", NA_character_,
  "Retained as a sequence-level proxy for aliphatic side-chain content and hydrophobic packing tendency.",

  "Surface", "Total SASA", "SASA_Total",
  "Total solvent-accessible surface area of the AlphaFold-predicted protein structure.",
  "Calculated from AlphaFold PDB coordinates using freesasa.",
  "continuous", "No", "Mass_Da",
  "Excluded because absolute SASA was strongly size-dependent; molecular mass was retained for protein size, and SASA hydrophobic ratio was retained for surface chemistry.",

  "Surface", "Apolar SASA", "SASA_Hydrophobic",
  "Apolar component of the solvent-accessible surface area.",
  "Calculated from AlphaFold PDB coordinates using freesasa residue/atom classification.",
  "continuous", "No", "SASA_Total",
  "Excluded because apolar SASA was nearly collinear with total SASA and protein size; SASA hydrophobic ratio was retained as the surface hydrophobicity descriptor.",

  "Surface", "Surface hydrophobicity", "SASA_Hydro_Ratio",
  "Fraction of solvent-accessible surface area classified as hydrophobic or apolar.",
  "Calculated as apolar SASA divided by total SASA from freesasa.",
  "fraction", "Yes", NA_character_,
  "Retained as the primary descriptor of exposed surface hydrophobicity.",

  "Surface", "Membrane-associated", "Is_Membrane",
  "Whether UniProt annotations indicate membrane-associated, cell-surface, lipid-anchor, GPI-anchor or transmembrane evidence.",
  "Set to 1 when the UniProt subcellular-location text contained membrane-associated or cell-surface terms after removing basement membrane matches, or when a TRANSMEM feature was present; otherwise set to 0.",
  "binary", "Yes", NA_character_,
  "Retained as a broad membrane/surface-associated annotation relevant to platform accessibility and hydrophobic-context effects.",

  "Surface", "Lipidation", "Has_Lipidation",
  "Whether UniProt annotates any lipid modification site for the protein.",
  "Set to 1 when the UniProt Lipidation field contained at least one LIPID feature; otherwise set to 0.",
  "binary", "Yes", NA_character_,
  "Retained as a binary lipid-modification marker; lipidation count was nearly identical and less interpretable because lipidation subtypes are heterogeneous.",

  "Surface", "Lipidation count", "Lipidation_Count",
  "Number of UniProt-annotated lipid modification features.",
  "Calculated as the number of LIPID feature records in the UniProt Lipidation field.",
  "count", "No", "Has_Lipidation",
  "Excluded because it was nearly collinear with the binary lipidation indicator, while the biological meaning of increasing lipidation count is subtype-dependent.",

  "Surface", "Transmembrane", "Is_Transmembrane",
  "Whether UniProt annotates at least one transmembrane region.",
  "Set to 1 when the UniProt Transmembrane field contained at least one TRANSMEM feature; otherwise set to 0.",
  "binary", "Yes", NA_character_,
  "Retained as a strict indicator of annotated transmembrane regions, complementary to the broader membrane-associated feature.",

  "Charge", "Net charge (pH 7.4)", "Net_Charge_7.4",
  "Estimated net charge of the full-length protein sequence at physiological pH.",
  "Calculated from the full-length UniProt sequence using Peptides::charge at pH 7.4 with the Lehninger pK scale; sequences containing non-standard amino-acid symbols were set to missing.",
  "continuous", "Yes", NA_character_,
  "Retained as the direct charge-state descriptor at physiological pH.",

  "Charge", "Isoelectric point", "pI_Value",
  "Estimated isoelectric point of the full-length protein sequence.",
  "Calculated from the full-length UniProt sequence using Peptides::pI with the Lehninger pK scale; sequences containing non-standard amino-acid symbols were set to missing.",
  "continuous", "No", "Net_Charge_7.4",
  "Excluded because pI was highly correlated with net charge at pH 7.4; net charge was retained as the more direct physiological charge descriptor.",

  "Charge", "Delta pI (pH 7.4)", "Delta_pI_7.4",
  "Absolute distance between the estimated isoelectric point and physiological pH 7.4.",
  "Calculated as abs(pI - 7.4), using the Peptides-derived pI value.",
  "continuous", "Yes", NA_character_,
  "Retained to represent proximity to charge-neutralization conditions near physiological pH, complementary to signed net charge at pH 7.4.",

  "Disorder", "Low-pLDDT fraction", "pLDDT_Fraction_Low",
  "Fraction of residues with low AlphaFold pLDDT, used as a proxy for low-confidence or disordered regions.",
  "Calculated as the fraction of C-alpha residues with AlphaFold pLDDT below 50.",
  "fraction", "Yes", NA_character_,
  "Retained as a structure-derived disorder proxy.",

  "Disorder", "Sequence instability", "Sequence_Instability",
  "Guruprasad dipeptide-based instability index of the full-length protein sequence.",
  "Calculated from the full-length UniProt sequence using Peptides::instaIndex; sequences shorter than two residues or containing non-standard amino-acid symbols were set to missing.",
  "continuous", "Yes", NA_character_,
  "Retained as a sequence-derived instability descriptor.",

  "Disorder", "Compositional bias", "CompBias_Fraction",
  "Fraction of the protein sequence covered by UniProt COMPBIAS annotations.",
  "Calculated as the union length of COMPBIAS coordinate intervals divided by full-length protein length; proteins without COMPBIAS annotations were assigned 0.",
  "fraction", "Yes", NA_character_,
  "Retained as the primary representation of compositional-bias burden.",

  "Secretory", "Glycosylation density", "Glyco_Density",
  "Length-normalized density of UniProt-annotated secretory glycosylation sites.",
  "Calculated as filtered CARBOHYD feature count divided by full-length protein length, after excluding glycation, in vitro glycation and variant-specific glycosylation annotations.",
  "density", "Yes", NA_character_,
  "Retained as the primary glycosylation descriptor because it reduces length dependence.",

  "Secretory", "Glycosylation sites", "Glyco_Count",
  "Number of UniProt-annotated secretory glycosylation sites.",
  "Calculated as the number of CARBOHYD feature records after excluding glycation, in vitro glycation and variant-specific glycosylation annotations.",
  "count", "No", "Glyco_Density",
  "Excluded because it was nearly collinear with glycosylation density; density was retained as the length-normalized representation.",

  "Secretory", "Disulfide bonds", "Disulfide_Count",
  "Number of UniProt-annotated disulfide features.",
  "Calculated as the number of DISULFID feature records in the UniProt Disulfide bond field.",
  "count", "Yes", NA_character_,
  "Retained as a count of secretory maturation and structural-constraint features.",

  "Secretory", "Secreted", "Is_Secreted",
  "Whether UniProt subcellular-location annotations indicate secreted or extracellular-space localization.",
  "Set to 1 when the UniProt subcellular-location text contained Secreted or extracellular space; otherwise set to 0.",
  "binary", "Yes", NA_character_,
  "Retained as a conservative secreted/extracellular localization marker.",

  "Secretory", "Propeptide", "Has_Propeptide",
  "Whether UniProt annotates a propeptide in the canonical precursor sequence.",
  "Set to 1 when the UniProt Propeptide field contained at least one PROPEP feature; otherwise set to 0.",
  "binary", "Yes", NA_character_,
  "Retained as a precursor-processing marker; this does not imply that the mature circulating protein retains the propeptide region.",

  "Secretory", "Propeptide count", "Propeptide_Count",
  "Number of UniProt-annotated propeptide features.",
  "Calculated as the number of PROPEP feature records in the UniProt Propeptide field.",
  "count", "No", "Has_Propeptide",
  "Excluded because it was nearly collinear with the binary propeptide indicator and is less interpretable for mature or platform-specific measured analytes.",

  "Secretory", "Propeptide fraction", "Propeptide_Fraction",
  "Fraction of the canonical precursor sequence covered by UniProt PROPEP annotations.",
  "Calculated as the union length of PROPEP coordinate intervals divided by full-length protein length.",
  "fraction", "No", "Has_Propeptide",
  "Excluded because it was nearly collinear with the binary propeptide indicator and propeptide length does not necessarily describe the mature measured analyte.",

  "Abundance", "Circulating abundance", "Log10_Abundance",
  "Baseline circulating protein abundance in blood.",
  "Calculated as log10-transformed HPA blood concentration, using MS-based concentration when available and immunoassay-based concentration otherwise.",
  "continuous", "Yes", NA_character_,
  "Retained as a baseline abundance descriptor."
)

safe_stat <- function(x, func, ...) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(NA_real_)
    func(x, ...)
}
statistics <- physchem_matrix %>%
    select(any_of(feature_decision$Feature)) %>%
    pivot_longer(everything(), names_to = "Feature", values_to = "Value") %>%
    group_by(Feature) %>%
    summarize(
        N_total          = n(),
        N_nonmissing     = sum(!is.na(Value)),
        N_missing        = sum(is.na(Value)),
        Missing_Fraction = mean(is.na(Value)),
        N_zero           = sum(Value == 0, na.rm = TRUE),
        N_nonzero        = sum(Value != 0, na.rm = TRUE),
        N_unique         = n_distinct(Value, na.rm = TRUE),
        Min              = safe_stat(Value, min),
        Q05              = safe_stat(Value, quantile, probs = 0.05, names = FALSE),
        Q25              = safe_stat(Value, quantile, probs = 0.25, names = FALSE),
        Median           = safe_stat(Value, median),
        Q75              = safe_stat(Value, quantile, probs = 0.75, names = FALSE),
        Q95              = safe_stat(Value, quantile, probs = 0.95, names = FALSE),
        Max              = safe_stat(Value, max),
        Unique_values_if_few = if_else(N_unique <= 10, paste(sort(unique(Value[!is.na(Value)])), collapse = "; "), NA_character_),
        .groups = "drop"
    )

physchem_dict <- feature_decision %>%
    left_join(statistics, by = "Feature") %>%
    mutate(Target_Feat = str_extract(Redundant_With, "^[^ /]+")) %>%
    rowwise() %>%
    mutate(
        Redundancy_Rho = if (!is.na(Target_Feat) && Target_Feat %in% colnames(cor_matrix) && Feature %in% rownames(cor_matrix)) {
            cor_matrix[Feature, Target_Feat]
        } else {
            NA_real_
        },
        Redundancy_P = if (!is.na(Target_Feat) && Target_Feat %in% colnames(cor_pmat) && Feature %in% rownames(cor_pmat)) {
            cor_pmat[Feature, Target_Feat]
        } else {
            NA_real_
        }
    ) %>%
    ungroup() %>%
    select(-Target_Feat)
fwrite(as.data.table(physchem_dict), "data/physchem_dictionary.tsv", sep = "\t", na = "NA")

# 8. Final validation summary ----

annotation_outputs <- c(
    "data/physchem_matrix.tsv.gz", "data/physchem_dictionary.tsv", "results/physchem_cor_matrix.tsv",
    "results/physchem_cor_pmat.tsv"
)
stopifnot(all(file.exists(annotation_outputs)))

required_hpa_columns <- c(
    "Log10_Abundance", "BloodConc_log10_pgml", "Abundance_Source", "HPA_Protein_Class", "HPA_Subcellular"
)
stopifnot(all(required_hpa_columns %in% names(physchem_matrix)))

cat("Annotated proteins:", nrow(physchem_matrix), "\n")
cat("Correlation features:", ncol(features_for_cor), "\n")
cat("Proteins with AlphaFold-derived features:", sum(!is.na(physchem_matrix$Radius_of_Gyration)), "\n")
cat("Retained model features:", sum(physchem_dict$Retained == "Yes", na.rm = TRUE), "\n")
cat("All expected outputs are present.\n")

# Recorded notebook output (reference only; regenerate by running this script):
# Annotated proteins: 12,008
# Correlation features: 30
# Proteins with AlphaFold-derived features: 11,959
# Retained model features: 22
# All expected outputs are present.
