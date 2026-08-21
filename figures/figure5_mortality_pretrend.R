# Figure 5: pre-treatment comparison for the parallel-trends assessment.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
load(file.path(root, "data", "derived", "mortality_results.RData"))

grid <- mortality_results$grid
plot_data <- data.frame(
  density = c(
    unlist(mortality_results$placebo$Female$density),
    unlist(mortality_results$placebo$Male$density)
  ),
  age = rep(grid, 4),
  sex = rep(c("Female", "Male"), each = 2 * length(grid)),
  endpoint = factor(rep(rep(c("Start", "End"), each = length(grid)), 2),
                    levels = c("Start", "End"))
)
figure <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(x = age, y = density, color = endpoint, linetype = endpoint)
) +
  ggplot2::geom_line() +
  ggplot2::facet_wrap(ggplot2::vars(sex)) +
  ggplot2::guides(color = ggplot2::guide_legend(title = NULL),
                  linetype = ggplot2::guide_legend(title = NULL)) +
  ggplot2::labs(x = "Age (years)", y = "Density") +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 16))

output <- file.path(root, "output", "figures", "figure5.pdf")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(output, figure, width = 10, height = 4, device = grDevices::cairo_pdf)
cat(sprintf("Saved %s\n", output))
