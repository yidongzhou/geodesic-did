# Figure 4: estimated GATT for female and male age-at-death distributions.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "R", "project_utils.R"))
assert_packages("ggplot2")
mortality_results <- load_single_object(
  file.path(root, "data", "derived", "mortality_results.RData"),
  "mortality_results"
)

grid <- mortality_results$grid
plot_data <- data.frame(
  density = c(
    unlist(mortality_results$main$Female$density),
    unlist(mortality_results$main$Male$density)
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

output <- file.path(root, "output", "figures", "figure4.pdf")
ensure_directory(dirname(output))
ggplot2::ggsave(output, figure, width = 10, height = 4, device = grDevices::cairo_pdf)
cat(sprintf("Saved %s; distances = %.6f (female), %.6f (male)\n", output,
            mortality_results$main$Female$distance,
            mortality_results$main$Male$distance))
