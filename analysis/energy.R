# 23 states that passed the liberalization law: 1996-2000
# AR, AZ, CA, CT, DE, IL, MA, MD, ME, MT, NH, NJ, NM, NV, NY, OH, OK, OR, PA, RI, TX, VA, WV
# 27 states that did not pass the liberalization law:
# AL, AK, CO, FL, GA, HI, IA, ID, IN, KS, KY, LA, MI, MN, MO, MS, ND, NE, NC, SC, SD, TN, UT, VT, WA, WI, WY
# Eight treated states later suspended deregulation (AZ, AR, CA, MT, NM, NV,
# OK, and VA); they remain in the treated group, matching the paper.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, 'R', 'gdd.R'))
source(file.path(root, 'R', 'lcm.R'))
state1 <- c("AZ", "AR", "CA", "CT", "DE", "IL", "ME", "MD", "MA", "MT", "NV",
            "NH", "NJ", "NM", "NY", "OH", "OK", "OR", "PA", "RI", "TX", "VA",
            "WV")
state0 <- c("AL", "AK", "CO", "FL", "GA", "HI", "ID", "IN", "IA", "KS", "KY",
            "LA", "MI", "MN", "MS", "MO", "NE", "NC", "ND", "SC", "SD", "TN",
            "UT", "VT", "WA", "WI", "WY")
tx <- ifelse(state.abb %in% state1, 1, 0)
energy <- readxl::read_excel(file.path(root, 'data', 'raw', 'annual_generation_state2023.xls'), skip = 1)
energy <- as.data.frame(energy)
# Follow the analysis convention of using absolute reported generation values.
energy$`GENERATION (Megawatthours)` <- abs(energy$`GENERATION (Megawatthours)`)
energy <- energy[energy$`TYPE OF PRODUCER` == 'Total Electric Power Industry' & energy$STATE %in% state.abb, c(1, 2, 4, 5)]

start <- as.integer(Sys.getenv("GDID_START_YEAR", unset = "1995"))
end <- as.integer(Sys.getenv("GDID_END_YEAR", unset = "2020"))
energys <- energy[energy$YEAR == start, ]
energye <- energy[energy$YEAR == end, ]
y00 <- y01 <- y10 <- y11 <- list()
# Fossil includes coal and petroleum; renewable includes conventional hydro,
# wind, wood and wood-derived fuels, geothermal, and solar.
for (i in 1:length(state.abb)) {
  energysi <- energys[energys$STATE == state.abb[i], ]
  compsi <- c(sum(energysi[energysi$`ENERGY SOURCE` %in% c("Coal", "Petroleum"), 4]),
              sum(energysi[energysi$`ENERGY SOURCE` %in% c("Nuclear"), 4]),
              sum(energysi[energysi$`ENERGY SOURCE` %in% c("Hydroelectric Conventional", "Wind",
                                                           "Wood and Wood Derived Fuels",
                                                           "Geothermal",
                                                           "Solar Thermal and Photovoltaic"), 4]))
  energyei <- energye[energye$STATE == state.abb[i], ]
  compei <- c(sum(energyei[energyei$`ENERGY SOURCE` %in% c("Coal", "Petroleum"), 4]),
              sum(energyei[energyei$`ENERGY SOURCE` %in% c("Nuclear"), 4]),
              sum(energyei[energyei$`ENERGY SOURCE` %in% c("Hydroelectric Conventional", "Wind",
                                                           "Wood and Wood Derived Fuels",
                                                           "Geothermal",
                                                           "Solar Thermal and Photovoltaic"), 4]))
  if(tx[i]) {
    y10[[state.abb[i]]] <- compsi / sum(compsi)
    y11[[state.abb[i]]] <- compei / sum(compei)
  } else {
    y00[[state.abb[i]]] <- compsi / sum(compsi)
    y01[[state.abb[i]]] <- compei / sum(compei)
  }
}
res <- gdd(y00 = y00, y01 = y01, y10 = y10, y11 = y11, optns = list(type = 'composition'))

df <- as.data.frame(rbind(matrix(unlist(y00), ncol = 3, byrow = TRUE),
                          matrix(unlist(y01), ncol = 3, byrow = TRUE),
                          matrix(unlist(y10), ncol = 3, byrow = TRUE),
                          matrix(unlist(y11), ncol = 3, byrow = TRUE)))
df$t <- c(rep(as.character(c(start, end)), c(length(y00), length(y01))),
          rep(as.character(c(start, end)), c(length(y10), length(y11))))
df$g <- rep(c('Monopoly states', 'Liberalized states'), c(2 * length(y00), 2 * length(y10)))
df$abb <- c(names(y00), names(y01), names(y10), names(y11))

energy_results <- list(
  start_year = start,
  end_year = end,
  y00 = y00,
  y01 = y01,
  y10 = y10,
  y11 = y11,
  estimate = res,
  points = df,
  distance = acos(sum(sqrt(res$gdd$start) * sqrt(res$gdd$end))),
  treated_states = state1,
  control_states = state0,
  categories = c("Fossil", "Nuclear", "Renewable")
)
save(
  energy_results,
  file = file.path(root, "data", "derived",
                   sprintf("energy_results_%d_%d.RData", start, end)),
  version = 2
)
cat(sprintf("Energy GATT distance for %d--%d: %.7f\n", start, end,
            energy_results$distance))
