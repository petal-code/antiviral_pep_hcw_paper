# =============================================================================
# 00_copy_sample_runs.R
#
# Copies 10 sample RDS files per arm (5 per scenario) from
# outputs/simulation/{arm}/ to output_samples/{arm}/ for sharing or
# quick inspection without copying the full simulation output.
# =============================================================================
library(here)

SRC_BASE <- here("outputs", "simulation")
DST_BASE <- here("output_samples")
N_PER_SCENARIO <- 5L   # 5 WestAfrica + 5 DRC = 10 per arm

arm_dirs <- list.dirs(SRC_BASE, full.names = FALSE, recursive = FALSE)

for (arm in arm_dirs) {
  src_dir <- file.path(SRC_BASE, arm)
  dst_dir <- file.path(DST_BASE, arm)
  dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
  
  files <- list.files(src_dir, pattern = "\\.rds$", full.names = FALSE)
  
  for (sc in c("WestAfrica", "DRC")) {
    sc_files <- sort(files[grepl(sprintf("^%s_", sc), files)])
    sample_files <- head(sc_files, N_PER_SCENARIO)
    
    for (f in sample_files) {
      file.copy(file.path(src_dir, f), file.path(dst_dir, f), overwrite = TRUE)
    }
  }
  
  message(sprintf("%s: copied %d files", arm,
                  length(list.files(dst_dir, pattern = "\\.rds$"))))
}

message("Done.")