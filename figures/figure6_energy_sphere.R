# Figure 6: electricity-generation compositions in the positive sphere octant.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, "R", "project_utils.R"))
assert_packages("scatterplot3d")
energy_results <- load_single_object(
  file.path(root, "data", "derived", "energy_results_1995_2020.RData"),
  "energy_results"
)

surface_angle <- seq(0, pi / 2, length.out = 42)
surface <- expand.grid(phi = surface_angle, theta = surface_angle)
surface$x <- sin(surface$phi) * cos(surface$theta)
surface$y <- sin(surface$phi) * sin(surface$theta)
surface$z <- cos(surface$phi)

draw_panel <- function(before, after, title, estimate = NULL) {
  projection <- scatterplot3d::scatterplot3d(
    surface$x, surface$y, surface$z,
    pch = 20, color = grDevices::adjustcolor("grey65", alpha.f = 0.18),
    xlim = c(0, 1), ylim = c(0, 1), zlim = c(0, 1),
    xlab = "Fossil", ylab = "Nuclear", zlab = "Renewable",
    angle = 42, box = FALSE, grid = TRUE, main = title
  )
  before <- sqrt(matrix(unlist(before), ncol = 3, byrow = TRUE))
  after <- sqrt(matrix(unlist(after), ncol = 3, byrow = TRUE))
  projection$points3d(before[, 1], before[, 2], before[, 3],
                      pch = 16, cex = 1.1, col = grDevices::adjustcolor("#F8766D", 0.7))
  projection$points3d(after[, 1], after[, 2], after[, 3],
                      pch = 18, cex = 1.2, col = grDevices::adjustcolor("#00BFC4", 0.7))
  if (!is.null(estimate)) {
    endpoints <- sqrt(rbind(estimate$gdd$start, estimate$gdd$end))
    projection$points3d(endpoints[1, 1], endpoints[1, 2], endpoints[1, 3],
                        pch = 16, cex = 2.2, col = "#F8766D")
    projection$points3d(endpoints[2, 1], endpoints[2, 2], endpoints[2, 3],
                        pch = 18, cex = 2.2, col = "#00BFC4")
  }
}

output <- file.path(root, "output", "figures", "figure6.pdf")
ensure_directory(dirname(output))
grDevices::cairo_pdf(output, width = 11, height = 5.5)
graphics::par(mfrow = c(1, 2), mar = c(3, 3, 3, 1))
draw_panel(energy_results$y10, energy_results$y11, "Liberalized states",
           energy_results$estimate)
draw_panel(energy_results$y00, energy_results$y01, "Monopoly states")
graphics::par(xpd = NA)
graphics::legend(
  "top", inset = c(0, -0.01), horiz = TRUE, bty = "n",
  legend = c("1995", "2020"),
  col = c("#F8766D", "#00BFC4"), pch = c(16, 18)
)
grDevices::dev.off()
cat(sprintf("Saved %s\n", output))
