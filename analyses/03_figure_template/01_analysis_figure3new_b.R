library(here)
library(ggplot2)
library(dplyr)

sdb <- readRDS(here("data-processed", "SDB_communityDeath_blended.rds"))

# With conflict: direct transform
sdb$coverage_conflict <- sdb$value * 80 / max(sdb$value)
sdb$dpc_conflict       <- 1 + 4 * (1 - sdb$value / max(sdb$value))

# Find the peak coverage day, restricted to day < 200
sub <- sdb[sdb$day < 200, ]
peak_row <- sub[which.max(sub$coverage_conflict), ]
peak_day <- peak_row$day

cat("peak_day =", peak_day, "\n")
cat("peak_coverage =", peak_row$coverage_conflict, "\n")

# Hold DPC at 1 up until peak_day (with conflict)
sdb$dpc_conflict[sdb$day <= peak_day] <- 1

# Without conflict: coverage same up to peak_day, then flat
sdb$coverage_noconflict <- sdb$coverage_conflict
sdb$coverage_noconflict[sdb$day > peak_day] <- peak_row$coverage_conflict

# Without conflict: DPC held at 1 throughout
sdb$dpc_noconflict <- 1

# ---- Plot coverage + DPC with secondary axis, axis colors matching line colors ----
scale_factor <- 5 / 80  # maps DPC range (1-5) onto coverage's visual scale (0-80)

cov_color <- "#E08214"
dpc_color <- "black"

p <- ggplot(sdb, aes(x = day)) +
  geom_line(aes(y = coverage_conflict), color = cov_color, linetype = "dashed", linewidth = 1.2) +
  geom_line(aes(y = coverage_noconflict), color = cov_color, linetype = "solid", linewidth = 1.2) +
  geom_line(aes(y = dpc_conflict / scale_factor), color = dpc_color, linetype = "dashed", linewidth = 1.0) +
  geom_line(aes(y = dpc_noconflict / scale_factor), color = dpc_color, linetype = "solid", linewidth = 1.0) +
  scale_y_continuous(
    name   = "Coverage (%)",
    sec.axis = sec_axis(~ . * scale_factor, name = "DPC (days)")
  ) +
  labs(x = "Days post ") +
  theme_minimal() +
  theme(
    axis.title.y       = element_text(color = cov_color),
    axis.text.y        = element_text(color = cov_color),
    axis.title.y.right  = element_text(color = dpc_color),
    axis.text.y.right   = element_text(color = dpc_color)
  )

print(p)