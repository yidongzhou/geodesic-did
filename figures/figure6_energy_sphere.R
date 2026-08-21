# Figure 6: electricity-generation compositions in the positive sphere octant.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

treated <- png::readPNG(
  file.path(root, "data", "derived", "figure6_treated.png"),
  native = TRUE
)
control <- png::readPNG(
  file.path(root, "data", "derived", "figure6_control.png"),
  native = TRUE
)
stopifnot(identical(dim(treated), dim(control)))

output <- file.path(root, "output", "figures", "figure6.pdf")
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
page_width <- 12
page_height <- page_width * nrow(treated) / (2 * ncol(treated))
grDevices::cairo_pdf(output, width = page_width, height = page_height)
grid::grid.newpage()
layout <- grid::grid.layout(nrow = 1, ncol = 2)
grid::pushViewport(grid::viewport(layout = layout))
grid::grid.raster(
  treated,
  width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
  interpolate = TRUE,
  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1)
)
grid::grid.raster(
  control,
  width = grid::unit(1, "npc"), height = grid::unit(1, "npc"),
  interpolate = TRUE,
  vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2)
)
grDevices::dev.off()
cat(sprintf("Saved %s\n", output))
