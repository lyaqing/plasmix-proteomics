# Shared Plasmix figure theme and color palettes.

theme_plasmix <- function() {
    list(
        ggplot2::theme_bw(base_family = "sans"),
        ggplot2::theme(
            # Plot background and margins
            plot.margin = ggplot2::margin(5, 5, 5, 5),
            plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
            panel.background = ggplot2::element_rect(fill = "transparent", color = NA),
            # Axes and panel border
            panel.border = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_line(color = "black", linewidth = 0.35),
            axis.line.y = ggplot2::element_line(color = "black", linewidth = 0.35),
            axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.35),
            axis.title = ggplot2::element_text(size = 9, color = "black"),
            axis.text = ggplot2::element_text(size = 8, color = "black"),
            # Grid lines
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
            # Facet strips
            strip.text = ggplot2::element_text(size = 8.5, face = "bold"),
            strip.background = ggplot2::element_rect(fill = "grey95", color = NA),
            # Legend
            legend.title = ggplot2::element_text(size = 7.5, face = "bold", vjust = 0.5),
            legend.text = ggplot2::element_text(size = 7.5, vjust = 0.5),
            legend.key.size = grid::unit(0.5, "lines"),
            legend.background = ggplot2::element_blank(),
            legend.spacing.x = grid::unit(0.5, "cm"),
            legend.spacing.y = grid::unit(0.1, "cm"),
            legend.margin = ggplot2::margin(5, 5, 5, 5),
            # Titles
            plot.title = ggplot2::element_text(face = "bold", size = 8.5, hjust = 0.5),
            plot.subtitle = ggplot2::element_text(face = "plain", size = 8.5, hjust = 0.5)
        ),
        ggplot2::coord_cartesian(clip = "off")
    )
}

# The default theme uses L-shaped axes. To draw a full panel border, add:
# ggplot2::theme(panel.border = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.35), axis.line.x = ggplot2::element_blank(), axis.line.y = ggplot2::element_blank())
# With clipping disabled above, linewidth = 0.35 visually matches the axis lines.

sample_order <- c("M", "Y", "P", "X", "F", "N")

sample_color <- c(
    M = "#3171b8", Y = "#398364", P = "#6d3390", X = "#F2B342", F = "#ae231c", N = "#6e6e6e",
    A = "#4cc3d9", E = "#f16745", A1 = "#4cc3d9", E1 = "#f16745", A2 = "#4cc4d96e", E2 = "#f1674569",
    PM = "#7bc8a4", BLK = "#febabd", QC = "#f780a1", CAL = "#e63568", NC = "#febabd", SC = "#f780a1",
    IPC = "#e63568", GW = "#c4c1c1", OF1 = "#dfd1c1", OF2 = "#c6a185", OF3 = "#876e4e", OF4 = "#654b3a", OF5 = "#270f0b"
)

platform_color <- c(OLK = "#489FA7", SOM = "#B33E90", NLS = "#F2B341", AAG = "#CC373A", DIA = "#155289")

batch_color <- c(
    OLK_P1_B1 = "#B8E2DE", OLK_P1_B2 = "#72BEB7", OLK_P2_B1 = "#389191", OLK_P2_B2 = "#226b6b",
    SOM_P1_B1 = "#D8AEC8", SOM_P1_B2 = "#d162b5", SOM_P2_B1 = "#b34b99", SOM_P2_B2 = "#873573",
    NLS_P1_B1 = "#F2B341", AAG_P1_B1 = "#db8f8f", AAG_P2_B1 = "#ce4c4e",
    DIA_P1_B1 = "#B9DBF4", DIA_P2_B1 = "#56A2C0", DIA_P3_B1 = "#95AAD3",
    DIA_P4_B1 = "#155289", DIA_P5_B1 = "#5C5D9E", DIA_P5_B2 = "#2e1667"
)

protocol_color <- c(
    "Olink Explore 384" = "#489FA7", "Olink Explore HT" = "#28827A",
    "SomaScan 11K" = "#D8AEC8", "Illumina Protein Prep" = "#B04D97",
    "NULISA CNS 120" = "#F2B341", "AAgAtlas" = "#CC373A", "AAgAtlas IgA/IgG" = "#CC373A",
    "Top14-E480" = "#3A68AE", "Top14-Astral" = "#155289", "Neat-Astral" = "#56A2C0",
    "SiO-Astral" = "#5C5D9E", "Mag-timsTOF" = "#584482"
)
