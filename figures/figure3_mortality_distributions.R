# Figure 3: age-at-death distributions before and after the Soviet collapse.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "R", "project_utils.R"))
assert_packages(c("ggplot2", "patchwork"))
mortality_results <- load_single_object(
  file.path(root, "data", "derived", "mortality_results.RData"),
  "mortality_results"
)

make_panel <- function(data, sex) {
  treated <- data[data$g == "Eastern Europe", ]
  # Match the control curves displayed in the submitted paper figure.
  control <- data[data$g == "Western Europe" & !data$abb %in% c("France", "UK"), ]
  ggplot2::ggplot() +
    ggplot2::geom_line(
      data = treated,
      ggplot2::aes(x = age, y = d, group = abb, color = abb)
    ) +
    ggplot2::geom_line(
      data = control,
      ggplot2::aes(x = age, y = d, group = abb),
      color = "black", alpha = 0.8
    ) +
    ggplot2::facet_grid(t ~ g) +
    ggplot2::guides(color = ggplot2::guide_legend(title = NULL, nrow = 1)) +
    ggplot2::labs(x = "Age (years)", y = "Density", title = sex) +
    ggplot2::theme_bw() +
    ggplot2::theme(text = ggplot2::element_text(size = 14), legend.position = "top")
}

figure <- patchwork::wrap_plots(
  make_panel(mortality_results$female_curves, "Female"),
  make_panel(mortality_results$male_curves, "Male"),
  ncol = 1
)
output <- file.path(root, "output", "figures", "figure3.pdf")
ensure_directory(dirname(output))
ggplot2::ggsave(output, figure, width = 10, height = 12, device = grDevices::cairo_pdf)
cat(sprintf("Saved %s\n", output))
