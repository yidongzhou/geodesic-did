# Figure 2: Monte Carlo errors for three outcome spaces.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "R", "project_utils.R"))
assert_packages("ggplot2")

em <- load_single_object(file.path(root, "data", "derived", "em.RData"), "em")
emmv <- load_single_object(file.path(root, "data", "derived", "emmv.RData"), "emmv")
en <- load_single_object(file.path(root, "data", "derived", "en.RData"), "en")
n_values <- c(50, 200, 1000)
stopifnot(identical(dim(em), c(500L, 3L)),
          identical(dim(emmv), c(500L, 3L)),
          identical(dim(en), c(500L, 3L)))

slopes <- vapply(
  list(em, emmv, en),
  function(errors) unname(coef(lm(log(colMeans(errors)) ~ log(n_values)))[2]),
  numeric(1)
)
stopifnot(isTRUE(all.equal(slopes, c(-0.4122140, -0.4730584, -0.5091448),
                                  tolerance = 1e-6)))

labels <- c(
  "Univariate probability distribution",
  "Multivariate probability distribution",
  "Network"
)
plot_data <- data.frame(
  error = c(em, emmv, en),
  n = factor(rep(rep(n_values, each = 500), 3), levels = n_values),
  outcome = factor(rep(labels, each = 500 * length(n_values)), levels = labels)
)

figure <- ggplot2::ggplot(plot_data, ggplot2::aes(x = n, y = error)) +
  ggplot2::geom_boxplot(outlier.alpha = 0.5) +
  ggplot2::facet_wrap(ggplot2::vars(outcome), scales = "free_y") +
  ggplot2::labs(x = "n", y = "Error") +
  ggplot2::theme_bw() +
  ggplot2::theme(text = ggplot2::element_text(size = 16))

output <- file.path(root, "output", "figures", "figure2.pdf")
ensure_directory(dirname(output))
ggplot2::ggsave(output, figure, width = 15, height = 6, device = grDevices::cairo_pdf)
cat(sprintf("Saved %s; slopes = %s\n", output,
            paste(sprintf("%.3f", slopes), collapse = ", ")))
