# Figure 7: ternary plots of electricity-generation compositions.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
suppressPackageStartupMessages(library(ggtern))
load(file.path(root, "data", "derived", "energy_results_1995_2020.RData"))

points <- energy_results$points
points$t <- factor(points$t, levels = c("1995", "2020"))
endpoints <- data.frame(
  V1 = c(energy_results$estimate$gdd$start[1], energy_results$estimate$gdd$end[1]),
  V2 = c(energy_results$estimate$gdd$start[2], energy_results$estimate$gdd$end[2]),
  V3 = c(energy_results$estimate$gdd$start[3], energy_results$estimate$gdd$end[3]),
  t = factor(c("1995", "2020"), levels = c("1995", "2020")),
  g = "Liberalized states"
)

figure <- suppressWarnings(ggtern::ggtern(points) +
  ggplot2::geom_point(
    ggtern::aes(x = V1, y = V2, z = V3, color = t, shape = t),
    size = 4, alpha = 0.6
  ) +
  ggplot2::geom_point(
    data = endpoints,
    ggtern::aes(x = V1, y = V2, z = V3, color = t, shape = t),
    size = 8
  ) +
  ggplot2::facet_wrap(ggplot2::vars(g)) +
  ggplot2::xlab("Fossil") + ggplot2::ylab("Nuclear") + ggtern::zlab("Renewable") +
  ggplot2::scale_color_discrete(name = NULL) +
  ggplot2::scale_shape_manual(values = c("1995" = 16, "2020" = 18), name = NULL) +
  ggtern::theme_rgbw() + ggtern::theme_nomask() +
  ggplot2::theme(
    legend.position = "inside",
    legend.position.inside = c(0, 1),
    legend.justification = c(0, 1.55),
    legend.box.just = "left",
    strip.text = ggplot2::element_text(size = 20)
  ))

output <- file.path(root, "output", "figures", "figure7.pdf")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(output, figure, width = 15, height = 7, device = grDevices::cairo_pdf)
cat(sprintf("Saved %s\n", output))
