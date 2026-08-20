# Generate Figures 2--7 from the archived RData files.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
figure_scripts <- file.path(
  root, "figures",
  c(
    "figure2_simulations.R",
    "figure3_mortality_distributions.R",
    "figure4_mortality_gatt.R",
    "figure5_mortality_pretrend.R",
    "figure6_energy_sphere.R",
    "figure7_energy_ternary.R"
  )
)

old_working_directory <- setwd(root)
on.exit(setwd(old_working_directory), add = TRUE)
for (script in figure_scripts) {
  status <- system2(
    file.path(R.home("bin"), "Rscript"),
    file.path("figures", basename(script))
  )
  if (status != 0L) stop(sprintf("Figure script failed: %s", basename(script)))
}
cat("Generated Figures 2--7 in output/figures/.\n")
