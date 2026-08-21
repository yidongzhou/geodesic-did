# Figure 2: Monte Carlo errors for three outcome spaces.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
load(file.path(root, "data", "derived", "em.RData"))
load(file.path(root, "data", "derived", "emmv.RData"))
load(file.path(root, "data", "derived", "en.RData"))
n_values <- c(50, 200, 1000)
n_repetitions <- nrow(em)

labels <- c(
  "Univariate probability distribution",
  "Multivariate probability distribution",
  "Network"
)
plot_data <- data.frame(
  error = c(em, emmv, en),
  n = factor(rep(rep(n_values, each = n_repetitions), 3), levels = n_values),
  outcome = factor(
    rep(labels, each = n_repetitions * length(n_values)),
    levels = labels
  )
)

figure <- ggplot2::ggplot(plot_data, ggplot2::aes(x = n, y = error)) +
  ggplot2::geom_boxplot(outlier.alpha = 0.5) +
  ggplot2::facet_wrap(ggplot2::vars(outcome), scales = "free_y") +
  ggplot2::labs(x = "n", y = "Error") +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 16))

output <- file.path(root, "output", "figures", "figure2.pdf")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
ggplot2::ggsave(output, figure, width = 15, height = 6, device = grDevices::cairo_pdf)
cat(sprintf("Saved %s\n", output))
