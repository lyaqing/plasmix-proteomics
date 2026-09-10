# Shared feature-to-protein mapping rules

# Analysis profiles retain named AAG/NLS analytes and combine SOMAmers within protein.
# Value is already log2: its arithmetic mean is a geometric mean on the linear scale.
# A SOM protein measurement is detected when a contributing assay exceeds its own LoD;
# the existing majority-of-technical-replicates criterion is subsequently applied.
aggregate_som_profiles <- function(data) {
    data <- as_tibble(data)
    som <- data %>% filter(Platform == "SOM", !is.na(UniProtID), UniProtID != "", !grepl("[|;]", UniProtID))
    keys <- intersect(c("Platform", "Batch", "Sample", "ColName", "ProcessLevel", "DataTier", "UniProtID"), names(som))
    if (!nrow(som)) return(data)
    som <- som %>% mutate(UniqueID = UniProtID,
                          IsDetectedMeasurement = is.finite(Value) & is.finite(LOD) & Value > LOD)
    text_cols <- intersect(c("RawName", "AssayID", "TargetName"), names(som))
    som <- som %>% group_by(across(all_of(keys))) %>%
        summarize(UniqueID = first(UniProtID),
                  Value = if (any(is.finite(Value))) mean(Value[is.finite(Value)]) else NA_real_,
                  LOD = NA_real_, IsDetectedMeasurement = any(IsDetectedMeasurement),
                  across(all_of(text_cols), ~ paste(sort(unique(.x[!is.na(.x)])), collapse = ";")),
                  across(any_of("Include"), ~ any(.x %in% TRUE)), .groups = "drop")
    bind_rows(data %>% filter(Platform != "SOM"), som)
}

analysis_feature_metadata <- function(feature_metadata) {
    data <- as_tibble(feature_metadata)
    som <- data %>% filter(Platform == "SOM", !Is_Protein_Group, !Is_Unknown) %>%
        group_by(Platform, UniProtID) %>%
        summarize(across(where(is.character) & !any_of(c("UniqueID", "Distinction_Key", "Suffix")),
                         ~ paste(sort(unique(.x[!is.na(.x)])), collapse = ";")),
                  UniqueID = first(UniProtID), Distinction_Key = first(UniProtID), Suffix = first(UniProtID),
                  Is_Protein_Group = FALSE, Is_Unknown = FALSE, .groups = "drop")
    bind_rows(data %>% filter(Platform != "SOM"), som)
}

normalize_protein_ids <- function(x) {
    vapply(x, function(value) {
        if (is.na(value) || value == "") return(NA_character_)
        ids <- trimws(unlist(strsplit(as.character(value), "[:;_|,]")))
        ids <- sort(unique(ids[ids != ""]))
        if (length(ids)) paste(ids, collapse = "|") else NA_character_
    }, character(1))
}

add_feature_key <- function(data) {
    data %>% mutate(Distinction_Key = case_when(
        Platform == "SOM" ~ AssayID,
        Platform == "OLK" ~ coalesce(UniProtID, TargetName),
        Platform %in% c("NLS", "AAG") ~ TargetName,
        Platform == "DIA" ~ coalesce(UniProtID, TargetName),
        TRUE ~ coalesce(AssayID, TargetName, UniProtID)
    ))
}

build_feature_mapping <- function(feature_rows, uniprot_info) {
    add_feature_key(feature_rows) %>%
        distinct(Platform, AssayID, TargetName, UniProtID, Distinction_Key) %>%
        left_join(uniprot_info, by = "UniProtID") %>%
        mutate(Is_Protein_Group = grepl("[|;]", UniProtID), Is_Unknown = is.na(UniProtID) | UniProtID == "") %>%
        group_by(Platform, UniProtID) %>%
        mutate(Is_Tau_Primary = UniProtID == "P10636" & TargetName %in% c("MAPT", "tTau"),
               Is_Tau_Variant = UniProtID == "P10636" & grepl("pTau", TargetName)) %>%
        arrange(desc(Is_Tau_Primary), Is_Tau_Variant, Distinction_Key) %>%
        mutate(Rank = row_number(), Suffix = Distinction_Key, .Assays_Per_Target = n_distinct(Distinction_Key),
               UniqueID = case_when(
                   Is_Unknown ~ Suffix,
                   Platform == "DIA" ~ UniProtID,
                   .Assays_Per_Target == 1 ~ UniProtID,
                   TRUE ~ paste0(UniProtID, "_", Suffix)
               )) %>%
        ungroup() %>% select(Platform, AssayID, TargetName, UniProtID, Distinction_Key, Protein_Full_Name,
                             Is_Tau_Primary, Is_Tau_Variant, Rank, Suffix, UniqueID, Is_Protein_Group, Is_Unknown)
}

# Use character() for all single-accession features; the default additionally requires one SOMAmer per accession within each batch.
get_batch_analysis_features <- function(feature_metadata, profile_data, strict_platforms = "SOM") {
    required_metadata <- c("Platform", "AssayID", "TargetName", "UniProtID", "UniqueID", "Is_Protein_Group", "Is_Unknown")
    required_profiles <- c("Platform", "Batch", "UniqueID")
    missing_metadata <- setdiff(required_metadata, names(feature_metadata))
    missing_profiles <- setdiff(required_profiles, names(profile_data))
    if (length(missing_metadata)) stop("Feature metadata is missing: ", paste(missing_metadata, collapse = ", "), call. = FALSE)
    if (length(missing_profiles)) stop("Protein profiles are missing: ", paste(missing_profiles, collapse = ", "), call. = FALSE)

    valid_metadata <- as_tibble(feature_metadata) %>%
        mutate(.Protein_Group = coalesce(as.logical(Is_Protein_Group), FALSE), .Unknown = coalesce(as.logical(Is_Unknown), FALSE),
               .Assay_Key = case_when(
                   Platform %in% c("SOM", "OLK") ~ as.character(AssayID),
                   Platform %in% c("NLS", "AAG") ~ as.character(TargetName),
                   Platform == "DIA" ~ as.character(UniProtID),
                   TRUE ~ coalesce(as.character(AssayID), as.character(TargetName), as.character(UniProtID))
               )) %>%
        filter(!.Protein_Group, !.Unknown, !is.na(.Assay_Key), .Assay_Key != "", !is.na(UniProtID), UniProtID != "") %>%
        distinct(Platform, UniqueID, UniProtID, .Assay_Key)

    profile_features <- as_tibble(profile_data) %>% distinct(Platform, Batch, UniqueID)
    standard_pairs <- profile_features %>%
        filter(!(Platform %in% strict_platforms)) %>%
        semi_join(valid_metadata %>% distinct(Platform, UniqueID), by = c("Platform", "UniqueID"))
    strict_pairs <- profile_features %>%
        filter(Platform %in% strict_platforms) %>%
        inner_join(valid_metadata %>% filter(Platform %in% strict_platforms), by = c("Platform", "UniqueID")) %>%
        group_by(Platform, Batch, UniProtID) %>%
        mutate(.Assays_Per_Batch_Target = n_distinct(.Assay_Key)) %>%
        ungroup() %>%
        filter(.Assays_Per_Batch_Target == 1) %>%
        select(Platform, Batch, UniqueID)
    bind_rows(standard_pairs, strict_pairs) %>% distinct(Platform, Batch, UniqueID)
}

filter_batch_analysis_features <- function(data, feature_metadata, strict_platforms = "SOM") {
    feature_pairs <- get_batch_analysis_features(feature_metadata, data, strict_platforms = strict_platforms)
    as_tibble(data) %>% semi_join(feature_pairs, by = c("Platform", "Batch", "UniqueID"))
}
