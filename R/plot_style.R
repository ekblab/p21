## ************************************************************************************************
## p21 bulk RNA-seq - Shared figure style
##
## Visual standard:
##   GraphPad Prism-inspired: white background, black axes, no gridlines,
##   clear data points, understated colors, unboxed legends on the right.
##   Based on the OTC miRNA QC and analysis reports (October 2026).
##
## Usage in subsequent analysis scripts:
##   source(here::here("R", "plot_style.R"))
##   ggplot(...) + p21_theme()
## ************************************************************************************************

## **** Nature Publishing Group (NPG) palette, matching ggsci::pal_npg("nrc")
## Named colors are stored locally so plotting does not require ggsci at runtime.
npg_colors <- c(red = "#E64B35", cyan = "#4DBBD5", green = "#00A087",
                navy = "#3C5488", coral = "#F39B7F", lavender = "#8491B4",
                mint = "#91D1C2", crimson = "#DC0000", brown = "#7E6148",
                beige = "#B09C85")

## Assignments are consistent across QC, DE, interaction and functional plots.
cell_colors <- c(HT1080 = npg_colors[["navy"]], MCF10A = npg_colors[["red"]])
genotype_colors <- c(wt = npg_colors[["navy"]], p21_neg = npg_colors[["red"]])
genotype_shapes <- c(wt = 16, p21_neg = 17)
radiation_colors <- c(control = npg_colors[["brown"]], `7gy` = npg_colors[["navy"]],
                      UVA = npg_colors[["red"]], UVB = npg_colors[["green"]])
processing_colors <- c(EKB = npg_colors[["navy"]], UKD = npg_colors[["red"]])

## **** Theme for scatter, PCA, dot and box plots
p21_theme <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size, base_family = "sans") +
    ggplot2::theme(
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.grid = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.5),
      axis.ticks = ggplot2::element_line(color = "black", linewidth = 0.45),
      axis.ticks.length = grid::unit(2.5, "mm"),
      axis.title = ggplot2::element_text(color = "black", size = base_size),
      axis.text = ggplot2::element_text(color = "black", size = base_size - 1),
      plot.title = ggplot2::element_text(face = "bold", color = "black", size = base_size + 1,
                                         hjust = 0, margin = ggplot2::margin(b = 4)),
      plot.subtitle = ggplot2::element_text(color = "#555555", size = base_size - 1,
                                            margin = ggplot2::margin(b = 9)),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", color = "black", size = base_size),
      panel.spacing = grid::unit(1.0, "lines"),
      legend.position = "right",
      legend.title = ggplot2::element_text(face = "plain", color = "black", size = base_size - 1),
      legend.text = ggplot2::element_text(color = "black", size = base_size - 1),
      legend.key = ggplot2::element_blank(),
      legend.background = ggplot2::element_blank(),
      legend.box.background = ggplot2::element_blank()
    )
}
