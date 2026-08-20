# 23 states that passed the liberalization law: 1996-2000
# AR, AZ, CA, CT, DE, IL, MA, MD, ME, MT, NH, NJ, NM, NV, NY, OH, OK, OR, PA, RI, TX, VA, WV
# 27 states that did not pass the liberalization law:
# AL, AK, CO, FL, GA, HI, IA, ID, IN, KS, KY, LA, MI, MN, MO, MS, ND, NE, NC, SC, SD, TN, UT, VT, WA, WI, WY
# 8 states that passed but later suspended deregulation:
# AZ, AR, CA, MT, NM, NV, OK, VA
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(root, 'R', 'gdd.R'))
source(file.path(root, 'R', 'lcm.R'))
library(ggplot2)
library(ggtern)
library(plotly)
# Discard implicit top-level plot printing; named outputs are saved explicitly.
grDevices::pdf(NULL)
state1 <- c("AZ", "AR", "CA", "CT", "DE", "IL", "ME", "MD", "MA", "MT", "NV",
            "NH", "NJ", "NM", "NY", "OH", "OK", "OR", "PA", "RI", "TX", "VA",
            "WV")
statee <- c("AZ", "AR", "CA", "MT", "NM", "NV", "OK", "VA")
state0 <- c("AL", "AK", "CO", "FL", "GA", "HI", "ID", "IN", "IA", "KS", "KY",
            "LA", "MI", "MN", "MS", "MO", "NE", "NC", "ND", "SC", "SD", "TN",
            "UT", "VT", "WA", "WI", "WY")
tx <- ifelse(state.abb %in% state1, 1, 0)
energy <- readxl::read_excel(file.path(root, 'data', 'raw', 'annual_generation_state2023.xls'), skip = 1)# 1990-2023
# https://www.eia.gov/electricity/
energy <- as.data.frame(energy)
unique(energy$`ENERGY SOURCE`)
table(energy[energy$`GENERATION (Megawatthours)` < 0 & energy$STATE %in% state.abb, c(1, 4)])
energy$`GENERATION (Megawatthours)` <- abs(energy$`GENERATION (Megawatthours)`)
energy <- energy[energy$`TYPE OF PRODUCER` == 'Total Electric Power Industry' & energy$STATE %in% state.abb, c(1, 2, 4, 5)]

td1 <- matrix(0, nrow = length(unique(energy$YEAR)), ncol = length(unique(energy$`ENERGY SOURCE`)[-1]))
td0 <- matrix(0, nrow = length(unique(energy$YEAR)), ncol = length(unique(energy$`ENERGY SOURCE`)[-1]))
for(i in 1:length(unique(energy$YEAR))) {
  for(j in 1:length(unique(energy$`ENERGY SOURCE`)[-1])) {
    td1ij <- rep(0, length(state1))
    idx <- state1 %in% energy[energy$`ENERGY SOURCE` == unique(energy$`ENERGY SOURCE`)[-1][j] &
                                energy$YEAR == unique(energy$YEAR)[i] & energy$STATE %in% state1, 2]
    td1ij[idx] <- energy[energy$`ENERGY SOURCE` == unique(energy$`ENERGY SOURCE`)[-1][j] &
                      energy$YEAR == unique(energy$YEAR)[i] & energy$STATE %in% state1, 4]
    td1[i, j] <- mean(td1ij / energy[energy$`ENERGY SOURCE` == 'Total' &
                                       energy$YEAR == unique(energy$YEAR)[i] & energy$STATE %in% state1, 4])
    td0ij <- rep(0, length(state0))
    idx <- state0 %in% energy[energy$`ENERGY SOURCE` == unique(energy$`ENERGY SOURCE`)[-1][j] &
                                energy$YEAR == unique(energy$YEAR)[i] & energy$STATE %in% state0, 2]
    td0ij[idx] <- energy[energy$`ENERGY SOURCE` == unique(energy$`ENERGY SOURCE`)[-1][j] &
                           energy$YEAR == unique(energy$YEAR)[i] & energy$STATE %in% state0, 4]
    td0[i, j] <- mean(td0ij / energy[energy$`ENERGY SOURCE` == 'Total' &
                                       energy$YEAR == unique(energy$YEAR)[i] & energy$STATE %in% state0, 4])
  }
}
unique(energy$`ENERGY SOURCE`)[-1][order(td1[34, ], decreasing = TRUE)]
unique(energy$`ENERGY SOURCE`)[-1][order(td0[1, ], decreasing = TRUE)]

ggplot(data = data.frame(y = c(td1, td0),
                         year = rep(rep(unique(energy$YEAR), length(unique(energy$`ENERGY SOURCE`)[-1])), 2),
                         source = rep(rep(unique(energy$`ENERGY SOURCE`)[-1], each = length(unique(energy$YEAR))), 2),
                         g = rep(c('Monopoly states', 'Liberalized states'), each = length(td1))),
       aes(x = year, y = y, color = source)) +
  geom_line() +
  facet_wrap(vars(g)) +
  geom_vline(xintercept = c(1990, 1995, 2020), linetype = 'dashed') +
  labs(x = 'Year', y = 'Energy source')

start <- as.integer(Sys.getenv("GDID_START_YEAR", unset = "1995"))
end <- as.integer(Sys.getenv("GDID_END_YEAR", unset = "2020"))
energys <- energy[energy$YEAR == start, ]
energye <- energy[energy$YEAR == end, ]
# table(energye[energye$`GENERATION (Megawatthours)` < 0 & energye$STATE %in% state.abb, c(1, 4)])
y00 <- y01 <- y10 <- y11 <- list()
for (i in 1:length(state.abb)) {
  energysi <- energys[energys$STATE == state.abb[i], ]
  compsi <- c(sum(energysi[energysi$`ENERGY SOURCE` %in% c("Coal", "Petroleum"), 4]),# "Natural Gas", "Other Gases", "Other Biomass"
              sum(energysi[energysi$`ENERGY SOURCE` %in% c("Nuclear"), 4]),
              sum(energysi[energysi$`ENERGY SOURCE` %in% c("Hydroelectric Conventional", "Wind",
                                                           "Wood and Wood Derived Fuels",
                                                           "Geothermal", # "Pumped Storage", "Other"
                                                           "Solar Thermal and Photovoltaic"), 4]))
  energyei <- energye[energye$STATE == state.abb[i], ]
  compei <- c(sum(energyei[energyei$`ENERGY SOURCE` %in% c("Coal", "Petroleum"), 4]),# "Natural Gas", "Other Gases", "Other Biomass"
              sum(energyei[energyei$`ENERGY SOURCE` %in% c("Nuclear"), 4]),#
              sum(energyei[energyei$`ENERGY SOURCE` %in% c("Hydroelectric Conventional", "Wind",
                                                           "Wood and Wood Derived Fuels",
                                                           "Geothermal", # "Pumped Storage", "Other"
                                                           "Solar Thermal and Photovoltaic"), 4]))
  if(tx[i]) {
  # if(tx[i] & !(state.abb[i] %in% statee)) {
    y10[[state.abb[i]]] <- compsi / sum(compsi)
    y11[[state.abb[i]]] <- compei / sum(compei)
  } else {
    y00[[state.abb[i]]] <- compsi / sum(compsi)
    y01[[state.abb[i]]] <- compei / sum(compei)
  }
}
# brct expects vectors on the sphere (sqrt of compositions)
rbind(brct(y10, list(type = 'composition')), brct(y11, list(type = 'composition')))
rbind(brct(y00, list(type = 'composition')), brct(y01, list(type = 'composition')))

res <- gdd(y00 = y00, y01 = y01, y10 = y10, y11 = y11, optns = list(type = 'composition'))
acos(sum(sqrt(res$gdd$start) * sqrt(res$gdd$end)))# geodesic angle
rbind(res$gdd$start, res$gdd$end)
acos(sum(sqrt(res$omega) * sqrt(res$gdd$end)))# 10 -> 11

df <- as.data.frame(rbind(matrix(unlist(y00), ncol = 3, byrow = TRUE),
                          matrix(unlist(y01), ncol = 3, byrow = TRUE),
                          matrix(unlist(y10), ncol = 3, byrow = TRUE),
                          matrix(unlist(y11), ncol = 3, byrow = TRUE)))
df$t <- c(rep(as.character(c(start, end)), c(length(y00), length(y01))),
          rep(as.character(c(start, end)), c(length(y10), length(y11))))
df$g <- rep(c('Monopoly states', 'Liberalized states'), c(2 * length(y00), 2 * length(y10)))
df$abb <- c(names(y00), names(y01), names(y10), names(y11))
ggtern(data = df) +
  geom_point(aes(x = V1, y = V2, z = V3, color = t, shape = t), size = 4, alpha = 0.6) + # vjust = 'top',
  # geom_text(aes(x = V1, y = V2, z = V3, color = t, label = abb), size = 2.5, position = position_jitter_tern(x = 0.05, y = 0.05, z = 0.1)) +
  geom_point(data = data.frame(V1 = c(res$gdd$start[1], res$gdd$end[1]),
                               V2 = c(res$gdd$start[2], res$gdd$end[2]),
                               V3 = c(res$gdd$start[3], res$gdd$end[3]),
                               t = as.character(c(start, end)),
                               g = rep('Liberalized states', 2)),
             aes(x = V1, y = V2, z = V3, color = t, shape = t), size = 8) +
  xlab('Fossil') + ylab('Nuclear') + zlab('Renewable') +
  scale_color_discrete(name = '') +
  scale_shape_manual(values = setNames(c(16, 18), c(start, end)), name = '') +
  theme_rgbw() +
  facet_wrap(vars(g)) +
  theme(legend.position = "inside",
        legend.position.inside = c(0, 1),
        legend.justification = c(0, 1.55), legend.box.just = "left",
        # legend.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 20)) +
  theme_nomask()
dir.create(file.path(root, "output", "legacy"), recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(root, 'output', 'legacy', sprintf('energyt_%d_%d.pdf', start, end)),
       width = 15, height = 7)

# Generate points on the surface of the 3D unit sphere
n_points <- 100
phi <- seq(0, pi/2, length.out = n_points)
theta <- seq(0, pi/2, length.out = n_points)
phi_grid <- rep(phi, each = n_points)
theta_grid <- rep(theta, n_points)
x <- sin(phi_grid) * cos(theta_grid)
y <- sin(phi_grid) * sin(theta_grid)
z <- cos(phi_grid)

# Create data frame for plotting
sphere_data <- data.frame(x = x, y = y, z = z)

# Create a plotly 3D scatter plot
plot_ly(sphere_data, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
        marker = list(size = 3, color = 'gray', opacity = 0.2), showlegend = FALSE) %>%
  add_trace(data = as.data.frame(matrix(sqrt(unlist(y10)), ncol = 3, byrow = TRUE)), x = ~V1, y = ~V2, z = ~V3,
            type = "scatter3d", mode = "markers", marker = list(color = "#F8766D", symbol = "circle", size = 10, opacity = 0.6)) %>%
  add_trace(data = as.data.frame(matrix(sqrt(unlist(y11)), ncol = 3, byrow = TRUE)), x = ~V1, y = ~V2, z = ~V3,
            type = "scatter3d", mode = "markers", marker = list(color = "#00BFC4", symbol = "diamond", size = 10, opacity = 0.6)) %>%
  add_trace(data = as.data.frame(t(sqrt(res$gdd$start))), x = ~V1, y = ~V2, z = ~V3,
            type = "scatter3d", mode = "markers", marker = list(color = "#F8766D", symbol = "circle", size = 20, opacity = 1)) %>%
  add_trace(data = as.data.frame(t(sqrt(res$gdd$end))), x = ~V1, y = ~V2, z = ~V3,
            type = "scatter3d", mode = "markers", marker = list(color = "#00BFC4", symbol = "diamond", size = 20, opacity = 1)) %>%
  layout(scene = list(xaxis = list(title = "Fossil", titlefont = list(size = 20), tickfont = list(size = 15)),
                      yaxis = list(title = "Nuclear", titlefont = list(size = 20), tickfont = list(size = 15)),
                      zaxis = list(title = "Renewable", titlefont = list(size = 20), tickfont = list(size = 15))))

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
plot_ly(sphere_data, x = ~x, y = ~y, z = ~z, type = 'scatter3d', mode = 'markers',
             marker = list(size = 3, color = 'gray', opacity = 0.2), showlegend = FALSE) %>%
  add_trace(data = as.data.frame(matrix(sqrt(unlist(y00)), ncol = 3, byrow = TRUE)), x = ~V1, y = ~V2, z = ~V3,
            type = "scatter3d", mode = "markers", marker = list(color = "#F8766D", symbol = "circle", size = 10, opacity = 0.6)) %>%
  add_trace(data = as.data.frame(matrix(sqrt(unlist(y01)), ncol = 3, byrow = TRUE)), x = ~V1, y = ~V2, z = ~V3,
            type = "scatter3d", mode = "markers", marker = list(color = "#00BFC4", symbol = "diamond", size = 10, opacity = 0.6)) %>%
  layout(scene = list(xaxis = list(title = "Fossil", titlefont = list(size = 20), tickfont = list(size = 15)),
                      yaxis = list(title = "Nuclear", titlefont = list(size = 20), tickfont = list(size = 15)),
                      zaxis = list(title = "Renewable", titlefont = list(size = 20), tickfont = list(size = 15))))
grDevices::dev.off()
