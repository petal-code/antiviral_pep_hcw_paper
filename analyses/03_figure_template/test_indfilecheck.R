x <- readRDS(here("outputs", "simulation","full_obv80", "WestAfrica_p007_r01.rds"))
names(x)
x$n_prevented
x$n_hcw_deaths
table(x$tdf$class)
str(x$prevented_completed)