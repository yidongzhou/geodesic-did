# https://www.mortality.org/Home/Index
# Belarus, Estonia, Latvia, Lithuania, Russia, Ukraine
# Bulgaria, Czechia, Hungary, Poland, Slovakia # Romania.
# Austria, Belgium, Croatia, Denmark, Finland, France, Germany,
# Greece, Iceland, Ireland, Italy, Luxembourg, Netherlands, Norway,
# Portugal, Slovenia, Spain, Sweden, Switzerland, U.K.
# 11 countries vs 20 countries

# Germany is constructed in memory by summing the East and West German `dx`
# columns within period and age, as documented in data/README.md.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(gsub("~+~", " ", sub("^--file=", "", script_arg[1]), fixed = TRUE), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
female_dir <- Sys.getenv("HMD_FEMALE_DIR")
male_dir <- Sys.getenv("HMD_MALE_DIR")
if (!dir.exists(female_dir) || !dir.exists(male_dir)) {
  stop("Set HMD_FEMALE_DIR and HMD_MALE_DIR to the HMD 1x5 life-table folders.")
}
source(file.path(root, 'R', 'gdd.R'))
source(file.path(root, 'R', 'lcm.R'))
library(ggplot2)
dir.create(file.path(root, "output", "legacy"), recursive = TRUE, showWarnings = FALSE)

yf00 <- list()
yf01 <- list()
yf10 <- list()
yf11 <- list()
ym00 <- list()
ym01 <- list()
ym10 <- list()
ym11 <- list()
df00 <- list()
df01 <- list()
df10 <- list()
df11 <- list()
dm00 <- list()
dm01 <- list()
dm10 <- list()
dm11 <- list()
tabb <- c('BLR', 'EST', 'LVA', 'LTU', 'RUS', 'UKR')#, 'HUN', 'POL', 'SVK'
cabb <- c('AUT', 'BEL', 'DNK', 'FIN', 'FRATNP', 'DEU', 'GRC', 'ISL', 'IRL', # total population is used
          'ITA', 'LUX', 'NLD', 'NOR', 'PRT', 'SVN', 'ESP', 'SWE', 'CHE', 'GBR_NP')# 'HRV' 2001-2004, 'BGR', 'CZE'
fn <- c(
  BLR = "Belarus", EST = "Estonia", LVA = "Latvia", LTU = "Lithuania",
  RUS = "Russia", UKR = "Ukraine", AUT = "Austria", BEL = "Belgium",
  DNK = "Denmark", FIN = "Finland", FRATNP = "France", DEU = "Germany",
  GRC = "Greece", ISL = "Iceland", IRL = "Ireland", ITA = "Italy",
  LUX = "Luxembourg", NLD = "Netherlands", NOR = "Norway",
  PRT = "Portugal", SVN = "Slovenia", ESP = "Spain", SWE = "Sweden",
  CHE = "Switzerland", GBR_NP = "UK"
)
codes <- c(cabb, tabb)
file_listf <- paste0(codes, '.fltper_1x5.txt')
file_listm <- paste0(codes, '.mltper_1x5.txt')
required_codes <- unlist(lapply(codes, function(code) {
  if (code == "DEU") c("DEUTE", "DEUTW") else code
}), use.names = FALSE)
required_file_listf <- paste0(required_codes, '.fltper_1x5.txt')
required_file_listm <- paste0(required_codes, '.mltper_1x5.txt')
missing_files <- c(
  file.path(female_dir, required_file_listf)[!file.exists(file.path(female_dir, required_file_listf))],
  file.path(male_dir, required_file_listm)[!file.exists(file.path(male_dir, required_file_listm))]
)
if (length(missing_files) > 0) {
  stop("Missing required HMD life-table files:\n", paste(missing_files, collapse = "\n"))
}
read_hmd_life_table <- function(directory, code, sex) {
  suffix <- if (sex == "female") ".fltper_1x5.txt" else ".mltper_1x5.txt"
  read_one <- function(country_code) {
    read.table(file.path(directory, paste0(country_code, suffix)),
               skip = 1, header = TRUE)
  }
  if (code != "DEU") return(read_one(code))

  east <- read_one("DEUTE")
  west <- read_one("DEUTW")
  if (!identical(as.character(east$Year), as.character(west$Year)) ||
      !identical(as.character(east$Age), as.character(west$Age))) {
    stop("East and West German HMD tables do not have matching Year/Age rows.")
  }
  east$dx <- east$dx + west$dx
  east
}
m <- 1001
dSup <- seq(0, 100, length.out = m)
qSup <- seq(0, 1, length.out = m)
adjust <- 2
bw <- 3
# 5 years
for(i in seq_along(file_listf)){
  abbi <- codes[i]
  fni <- unname(fn[[abbi]])
  if(file_listf[i] %in% paste0(cabb, '.fltper_1x5.txt')) {
    lt <- read_hmd_life_table(female_dir, abbi, "female")
    lt0 <- lt[lt$Year == '1985-1989', ][1:100, 'dx']
    df00[[fni]] <- frechet::CreateDensity(freq = lt0, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
                                                       # userBwMu = bw))$y
    # lt0 <- lt[lt$Year == '1985-1989', ][1:101, 'dx']
    # df00[[fni]] <- density(x = rep(0:100, lt0), adjust = adjust, from = 0, to = 100, n = m)$y
    yf00[[fni]] <- fdadensity::dens2quantile(dens = df00[[fni]], dSup = dSup, qSup = qSup)
    lt1 <- lt[lt$Year == '1990-1994', ][1:100, 'dx']
    df01[[fni]] <- frechet::CreateDensity(freq = lt1, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
                                                       # userBwMu = bw))$y
    # lt1 <- lt[lt$Year == '1990-1994', ][1:101, 'dx']
    # df01[[fni]] <- density(x = rep(0:100, lt1), adjust = adjust, from = 0, to = 100, n = m)$y
    yf01[[fni]] <- fdadensity::dens2quantile(dens = df01[[fni]], dSup = dSup, qSup = qSup)
  } else if(file_listf[i] %in% paste0(tabb, '.fltper_1x5.txt')) {
    lt <- read_hmd_life_table(female_dir, abbi, "female")
    lt0 <- lt[lt$Year == '1985-1989', ][1:100, 'dx']
    df10[[fni]] <- frechet::CreateDensity(freq = lt0, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
                                                       # userBwMu = bw))$y
    # lt0 <- lt[lt$Year == '1985-1989', ][1:101, 'dx']
    # df10[[fni]] <- density(x = rep(0:100, lt0), adjust = adjust, from = 0, to = 100, n = m)$y
    yf10[[fni]] <- fdadensity::dens2quantile(dens = df10[[fni]], dSup = dSup, qSup = qSup)
    lt1 <- lt[lt$Year == '1990-1994', ][1:100, 'dx']
    df11[[fni]] <- frechet::CreateDensity(freq = lt1, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
                                                       # userBwMu = bw))$y
    # lt1 <- lt[lt$Year == '1990-1994', ][1:101, 'dx']
    # df11[[fni]] <- density(x = rep(0:100, lt1), adjust = adjust, from = 0, to = 100, n = m)$y
    yf11[[fni]] <- fdadensity::dens2quantile(dens = df11[[fni]], dSup = dSup, qSup = qSup)
  }
  if(file_listm[i] %in% paste0(cabb, '.mltper_1x5.txt')) {
    lt <- read_hmd_life_table(male_dir, abbi, "male")
    lt0 <- lt[lt$Year == '1985-1989', ][1:100, 'dx']
    dm00[[fni]] <- frechet::CreateDensity(freq = lt0, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
                                                       # userBwMu = bw))$y
    # lt0 <- lt[lt$Year == '1985-1989', ][1:101, 'dx']
    # dm00[[fni]] <- density(x = rep(0:100, lt0), adjust = adjust, from = 0, to = 100, n = m)$y
    ym00[[fni]] <- fdadensity::dens2quantile(dens = dm00[[fni]], dSup = dSup, qSup = qSup)
    lt1 <- lt[lt$Year == '1990-1994', ][1:100, 'dx']
    dm01[[fni]] <- frechet::CreateDensity(freq = lt1, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
                                                       # userBwMu = bw))$y
    # lt1 <- lt[lt$Year == '1990-1994', ][1:101, 'dx']
    # dm01[[fni]] <- density(x = rep(0:100, lt1), adjust = adjust, from = 0, to = 100, n = m)$y
    ym01[[fni]] <- fdadensity::dens2quantile(dens = dm01[[fni]], dSup = dSup, qSup = qSup)
  } else if(file_listm[i] %in% paste0(tabb, '.mltper_1x5.txt')) {
    lt <- read_hmd_life_table(male_dir, abbi, "male")
    lt0 <- lt[lt$Year == '1985-1989', ][1:100, 'dx']
    dm10[[fni]] <- frechet::CreateDensity(freq = lt0, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
                                                       # userBwMu = bw))$y
    # lt0 <- lt[lt$Year == '1985-1989', ][1:101, 'dx']
    # dm10[[fni]] <- density(x = rep(0:100, lt0), adjust = adjust, from = 0, to = 100, n = m)$y
    ym10[[fni]] <- fdadensity::dens2quantile(dens = dm10[[fni]], dSup = dSup, qSup = qSup)
    lt1 <- lt[lt$Year == '1990-1994', ][1:100, 'dx']
    dm11[[fni]] <- frechet::CreateDensity(freq = lt1, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
                                                       # userBwMu = bw))$y
    # lt1 <- lt[lt$Year == '1990-1994', ][1:101, 'dx']
    # dm11[[fni]] <- density(x = rep(0:100, lt1), adjust = adjust, from = 0, to = 100, n = m)$y
    ym11[[fni]] <- fdadensity::dens2quantile(dens = dm11[[fni]], dSup = dSup, qSup = qSup)
  }
}

dff <- data.frame(d = c(unlist(df00), unlist(df01), unlist(df10), unlist(df11)),
                 age = rep(dSup, 2 * (length(cabb) + length(tabb))),
                 t = c(rep(c('1985-1989', '1990-1994'), each = m * length(cabb)),
                       rep(c('1985-1989', '1990-1994'), each = m * length(tabb))),
                 g = rep(c('Western Europe', 'Eastern Europe'), c(2 * m * length(cabb), 2 * m * length(tabb))),
                 abb = rep(c(names(df00), names(df01), names(df10), names(df11)), each = m))
pdf(file = file.path(root, "output", "legacy", "fd.pdf"), width = 10, height = 6)
ggplot() +# [dff$abb %in% c("Bulgaria", "Czechia", "Hungary", "Poland", "Slovakia"), ]
  # geom_smooth(aes(x = age, y = d, group = abb), span = 0.3, se = FALSE, color = 'black', linewidth = 0.3, alpha = 0.5) +
  geom_line(data = dff[dff$abb %in%  fn[names(fn) %in% tabb], ], aes(x = age, y = d, group = abb, color = abb)) +
  guides(color = guide_legend(title = NULL, nrow = 1)) +
  # scale_color_manual(values = c("Western Europe" = "black",
  #                               "Belarus" = "red",
  #                               "Estonia" = "blue",
  #                               "Lithuania" = "green",
  #                               "Latvia" = "purple",
  #                               "Russia" = "orange",
  #                               "Ukraine" = "brown"), name = '') +
  geom_line(data = dff[dff$abb %in%  fn[names(fn) %in% cabb], ], aes(x = age, y = d, group = abb)) +
  facet_grid(t ~ g) +
  labs(x = "Age (years)", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = 'top')
dev.off()
dfm <- data.frame(d = c(unlist(dm00), unlist(dm01), unlist(dm10), unlist(dm11)),
                 age = rep(dSup, 2 * (length(cabb) + length(tabb))),
                 t = c(rep(c('1985-1989', '1990-1994'), each = m * length(cabb)),
                       rep(c('1985-1989', '1990-1994'), each = m * length(tabb))),
                 g = rep(c('Western Europe', 'Eastern Europe'), c(2 * m * length(cabb), 2 * m * length(tabb))),
                 abb = rep(c(names(dm00), names(dm01), names(dm10), names(dm11)), each = m))
pdf(file = file.path(root, "output", "legacy", "md.pdf"), width = 10, height = 6)
ggplot() +
  geom_line(data = dfm[dfm$abb %in%  fn[names(fn) %in% tabb], ], aes(x = age, y = d, group = abb,
                color = abb)) +
  guides(color = guide_legend(title = NULL, nrow = 1)) +
  geom_line(data = dfm[dfm$abb %in%  fn[names(fn) %in% cabb], ], aes(x = age, y = d, group = abb)) +
  facet_grid(t ~ g) +
  labs(x = "Age (years)", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 20), legend.position = 'top')
dev.off()

resf <- gdd(y00 = yf00, y01 = yf01, y10 = yf10, y11 = yf11, optns = list(type = 'measure', lower = 0, upper = 100))
resm <- gdd(y00 = ym00, y01 = ym01, y10 = ym10, y11 = ym11, optns = list(type = 'measure', lower = 0, upper = 100))
sqrt(mean((resf$gdd$end - resf$gdd$start)^2))# 6: 2.11, 9: 1.35, 11: 1.32,
sqrt(mean((resm$gdd$end - resm$gdd$start)^2))# 6: 4.26, 9: 2.86, 11: 2.65,
mean(abs(resf$gdd$end - resf$gdd$start))
mean(abs(resm$gdd$end - resm$gdd$start))
100 * mean(abs(resm$gdd$end - resm$gdd$start))
df <- list(start = density(x = resf$gdd$start, adjust = adjust, from = 0, to = 100, n = m)$y,
           end = density(x = resf$gdd$end, adjust = adjust, from = 0, to = 100, n = m)$y)
dm <- list(start = density(x = resm$gdd$start, adjust = adjust, from = 0, to = 100, n = m)$y,
           end = density(x = resm$gdd$end, adjust = adjust, from = 0, to = 100, n = m)$y)
main_effect <- list(
  Female = list(result = resf, density = df,
                distance = sqrt(mean((resf$gdd$end - resf$gdd$start)^2))),
  Male = list(result = resm, density = dm,
              distance = sqrt(mean((resm$gdd$end - resm$gdd$start)^2)))
)

# density
pdf(file = file.path(root, "output", "legacy", "dse.pdf"), width = 10, height = 4)
ggplot(data = data.frame(d = c(unlist(df), unlist(dm)),
                         age = rep(dSup, 2 * 2),
                         s = rep(c('Female', 'Male'), each = 2 * m),
                         g = factor(rep(rep(c('Start', 'End'), each = m), 2), levels = c('Start', 'End')))) +
  geom_line(aes(x = age, y = d, color = g, linetype = g)) +
  guides(color = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) +
  facet_wrap(vars(s)) +
  labs(x = "Age (years)", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()

# 1980-1984 vs 1985-1989
for(i in seq_along(file_listf)){
  abbi <- codes[i]
  fni <- unname(fn[[abbi]])
  if(file_listf[i] %in% paste0(cabb, '.fltper_1x5.txt')) {
    lt <- read_hmd_life_table(female_dir, abbi, "female")
    if(abbi == "GRC") {
      lt0 <- lt[lt$Year == '1981-1984', ][1:100, 'dx']
    } else if(abbi == "SVN") {
      lt0 <- lt[lt$Year == '1983-1984', ][1:100, 'dx']
    } else {
      lt0 <- lt[lt$Year == '1980-1984', ][1:100, 'dx']
    }
    df00[[fni]] <- frechet::CreateDensity(freq = lt0, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
    # userBwMu = bw))$y
    # lt0 <- lt[lt$Year == '1980-1984', ][1:101, 'dx']
    # df00[[fni]] <- density(x = rep(0:100, lt0), adjust = adjust, from = 0, to = 100, n = m)$y
    yf00[[fni]] <- fdadensity::dens2quantile(dens = df00[[fni]], dSup = dSup, qSup = qSup)
    lt1 <- lt[lt$Year == '1985-1989', ][1:100, 'dx']
    df01[[fni]] <- frechet::CreateDensity(freq = lt1, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
    # userBwMu = bw))$y
    # lt1 <- lt[lt$Year == '1985-1989', ][1:101, 'dx']
    # df01[[fni]] <- density(x = rep(0:100, lt1), adjust = adjust, from = 0, to = 100, n = m)$y
    yf01[[fni]] <- fdadensity::dens2quantile(dens = df01[[fni]], dSup = dSup, qSup = qSup)
  } else if(file_listf[i] %in% paste0(tabb, '.fltper_1x5.txt')) {
    lt <- read_hmd_life_table(female_dir, abbi, "female")
    lt0 <- lt[lt$Year == '1980-1984', ][1:100, 'dx']
    df10[[fni]] <- frechet::CreateDensity(freq = lt0, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
    # userBwMu = bw))$y
    # lt0 <- lt[lt$Year == '1980-1984', ][1:101, 'dx']
    # df10[[fni]] <- density(x = rep(0:100, lt0), adjust = adjust, from = 0, to = 100, n = m)$y
    yf10[[fni]] <- fdadensity::dens2quantile(dens = df10[[fni]], dSup = dSup, qSup = qSup)
    lt1 <- lt[lt$Year == '1985-1989', ][1:100, 'dx']
    df11[[fni]] <- frechet::CreateDensity(freq = lt1, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
    # userBwMu = bw))$y
    # lt1 <- lt[lt$Year == '1985-1989', ][1:101, 'dx']
    # df11[[fni]] <- density(x = rep(0:100, lt1), adjust = adjust, from = 0, to = 100, n = m)$y
    yf11[[fni]] <- fdadensity::dens2quantile(dens = df11[[fni]], dSup = dSup, qSup = qSup)
  }
  if(file_listm[i] %in% paste0(cabb, '.mltper_1x5.txt')) {
    lt <- read_hmd_life_table(male_dir, abbi, "male")
    if(abbi == "GRC") {
      lt0 <- lt[lt$Year == '1981-1984', ][1:100, 'dx']
    } else if(abbi == "SVN") {
      lt0 <- lt[lt$Year == '1983-1984', ][1:100, 'dx']
    } else {
      lt0 <- lt[lt$Year == '1980-1984', ][1:100, 'dx']
    }
    dm00[[fni]] <- frechet::CreateDensity(freq = lt0, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
    # userBwMu = bw))$y
    # lt0 <- lt[lt$Year == '1980-1984', ][1:101, 'dx']
    # dm00[[fni]] <- density(x = rep(0:100, lt0), adjust = adjust, from = 0, to = 100, n = m)$y
    ym00[[fni]] <- fdadensity::dens2quantile(dens = dm00[[fni]], dSup = dSup, qSup = qSup)
    lt1 <- lt[lt$Year == '1985-1989', ][1:100, 'dx']
    dm01[[fni]] <- frechet::CreateDensity(freq = lt1, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
    # userBwMu = bw))$y
    # lt1 <- lt[lt$Year == '1985-1989', ][1:101, 'dx']
    # dm01[[fni]] <- density(x = rep(0:100, lt1), adjust = adjust, from = 0, to = 100, n = m)$y
    ym01[[fni]] <- fdadensity::dens2quantile(dens = dm01[[fni]], dSup = dSup, qSup = qSup)
  } else if(file_listm[i] %in% paste0(tabb, '.mltper_1x5.txt')) {
    lt <- read_hmd_life_table(male_dir, abbi, "male")
    lt0 <- lt[lt$Year == '1980-1984', ][1:100, 'dx']
    dm10[[fni]] <- frechet::CreateDensity(freq = lt0, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
    # userBwMu = bw))$y
    # lt0 <- lt[lt$Year == '1980-1984', ][1:101, 'dx']
    # dm10[[fni]] <- density(x = rep(0:100, lt0), adjust = adjust, from = 0, to = 100, n = m)$y
    ym10[[fni]] <- fdadensity::dens2quantile(dens = dm10[[fni]], dSup = dSup, qSup = qSup)
    lt1 <- lt[lt$Year == '1985-1989', ][1:100, 'dx']
    dm11[[fni]] <- frechet::CreateDensity(freq = lt1, bin = 0:100,
                                          optns = list(outputGrid = seq(0, 100, length.out = m)))$y#,
    # userBwMu = bw))$y
    # lt1 <- lt[lt$Year == '1985-1989', ][1:101, 'dx']
    # dm11[[fni]] <- density(x = rep(0:100, lt1), adjust = adjust, from = 0, to = 100, n = m)$y
    ym11[[fni]] <- fdadensity::dens2quantile(dens = dm11[[fni]], dSup = dSup, qSup = qSup)
  }
}

resf <- gdd(y00 = yf00, y01 = yf01, y10 = yf10, y11 = yf11, optns = list(type = 'measure', lower = 0, upper = 100))
resm <- gdd(y00 = ym00, y01 = ym01, y10 = ym10, y11 = ym11, optns = list(type = 'measure', lower = 0, upper = 100))
sqrt(mean((resf$gdd$end - resf$gdd$start)^2))
sqrt(mean((resm$gdd$end - resm$gdd$start)^2))
mean(abs(resf$gdd$end - resf$gdd$start))
mean(abs(resm$gdd$end - resm$gdd$start))
100 * mean(abs(resm$gdd$end - resm$gdd$start))
df <- list(start = density(x = resf$gdd$start, adjust = adjust, from = 0, to = 100, n = m)$y,
           end = density(x = resf$gdd$end, adjust = adjust, from = 0, to = 100, n = m)$y)
dm <- list(start = density(x = resm$gdd$start, adjust = adjust, from = 0, to = 100, n = m)$y,
           end = density(x = resm$gdd$end, adjust = adjust, from = 0, to = 100, n = m)$y)
placebo_effect <- list(
  Female = list(result = resf, density = df,
                distance = sqrt(mean((resf$gdd$end - resf$gdd$start)^2))),
  Male = list(result = resm, density = dm,
              distance = sqrt(mean((resm$gdd$end - resm$gdd$start)^2)))
)

# density
pdf(file = file.path(root, "output", "legacy", "dse8089.pdf"), width = 10, height = 4)
ggplot(data = data.frame(d = c(unlist(df), unlist(dm)),
                         age = rep(dSup, 2 * 2),
                         s = rep(c('Female', 'Male'), each = 2 * m),
                         g = factor(rep(rep(c('Start', 'End'), each = m), 2), levels = c('Start', 'End')))) +
  geom_line(aes(x = age, y = d, color = g, linetype = g)) +
  guides(color = guide_legend(title = NULL), linetype = guide_legend(title = NULL)) +
  facet_wrap(vars(s)) +
  labs(x = "Age (years)", y = "Density") +
  theme_bw() +
  theme(text = element_text(size = 20))
dev.off()

mortality_results <- list(
  grid = dSup,
  female_curves = dff,
  male_curves = dfm,
  main = main_effect,
  placebo = placebo_effect,
  treated_codes = tabb,
  control_codes = cabb
)
stopifnot(
  isTRUE(all.equal(main_effect$Female$distance, 2.114461517, tolerance = 1e-7)),
  isTRUE(all.equal(main_effect$Male$distance, 4.258858158, tolerance = 1e-7)),
  isTRUE(all.equal(placebo_effect$Female$distance, 0.586718557, tolerance = 1e-7)),
  isTRUE(all.equal(placebo_effect$Male$distance, 0.932983344, tolerance = 1e-7))
)
save(mortality_results,
     file = file.path(root, "data", "derived", "mortality_results.RData"),
     version = 2)
cat(sprintf(
  "GATT distances: female %.6f, male %.6f; placebo: female %.6f, male %.6f\n",
  main_effect$Female$distance, main_effect$Male$distance,
  placebo_effect$Female$distance, placebo_effect$Male$distance
))
