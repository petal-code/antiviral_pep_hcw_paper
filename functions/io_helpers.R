# io_helpers.R
# -----------------------------------------------------------------------------
# Generic filesystem helpers for locating analysis artefacts on disk. These have
# NO dependency on any model code (or on each other's globals), so this file is
# safe to source first from any analysis as a base utility layer.
#
# Contents
#   list_files_matching() : files in a directory whose name matches a string,
#                           returned newest-first. Never errors on "no match"
#                           (returns character(0)) -- use when you want the N
#                           most recent, or just to count.
#   find_latest_file()    : the single most recent matching file. Thin wrapper
#                           around list_files_matching(); errors by default if
#                           nothing matches, so `readRDS(find_latest_file(...))`
#                           fails loudly with a clear message rather than on a
#                           mysterious NA path.
#
# Typical use -- pull the most recent West Africa fit out of outputs/:
#   f <- find_latest_file(here::here("outputs"), "WestAfrica")
#   fit <- readRDS(f)
#
# "Most recent" by name vs by mtime
# ---------------------------------
# Default sorting is by = "name" (lexicographic on the file NAME). This assumes
# the timestamp is embedded in the name in a sortable form, e.g.
#   fiber_ABC_SMC_Worst_WestAfrica_2026-05-26.rds
# (YYYY-MM-DD and YYYYMMDD_HHMMSS both sort chronologically as text). Name-based
# sorting is deliberately the default because it is robust to a fresh `git
# clone`/checkout, which resets every tracked file's modification time to the
# clone time -- so for committed artefacts mtime carries no chronological
# signal. Pass by = "mtime" when names are NOT timestamped and you genuinely
# want "most recently written" (correct for files generated in the current
# session, unreliable for freshly-checked-out ones).
# -----------------------------------------------------------------------------


# List files in `dir` whose NAME matches `pattern`, sorted so the newest is
# first.
#   pattern     : substring matched against the file name (fixed = TRUE, the
#                 default) or a regular expression (fixed = FALSE). NULL or ""
#                 matches every file.
#   by          : "name" (default) sorts lexicographically by file name; "mtime"
#                 sorts by last-modified time. See the header note on which to
#                 use.
#   fixed       : TRUE treats `pattern` as a literal substring; FALSE as a regex.
#   ignore_case : case-insensitive matching.
#   recursive   : recurse into sub-directories of `dir`.
#   full_names  : return full paths (default) or bare file names.
# Returns a character vector, newest-first, length 0 if nothing matches. Only
# regular files are returned (sub-directories are never included).
list_files_matching <- function(dir, pattern = NULL,
                                by = c("name", "mtime"),
                                fixed = TRUE, ignore_case = FALSE,
                                recursive = FALSE, full_names = TRUE) {
  by <- match.arg(by)
  if (!dir.exists(dir)) {
    stop("Directory does not exist: ", dir, call. = FALSE)
  }

  files <- list.files(dir, full.names = TRUE, recursive = recursive)
  files <- files[!dir.exists(files)]            # keep regular files only
  if (length(files) == 0L) return(character(0))

  if (!is.null(pattern) && nzchar(pattern)) {
    nm <- basename(files)
    keep <- if (fixed) {
      if (ignore_case) grepl(tolower(pattern), tolower(nm), fixed = TRUE)
      else             grepl(pattern, nm, fixed = TRUE)
    } else {
      grepl(pattern, nm, ignore.case = ignore_case)
    }
    files <- files[keep]
  }
  if (length(files) == 0L) return(character(0))

  ord <- if (by == "mtime") {
    order(file.info(files)$mtime, decreasing = TRUE)
  } else {
    order(basename(files), decreasing = TRUE)
  }
  files <- files[ord]

  if (!full_names) files <- basename(files)
  files
}


# Most recent single file in `dir` matching `pattern`. Arguments are passed
# straight through to list_files_matching(); see there for their meaning.
#   error : TRUE (default) stops with an informative message if nothing matches;
#           FALSE returns NA_character_ instead.
# Returns a single file path (or NA_character_ when error = FALSE and there is
# no match).
find_latest_file <- function(dir, pattern = NULL,
                             by = c("name", "mtime"),
                             fixed = TRUE, ignore_case = FALSE,
                             recursive = FALSE, error = TRUE) {
  hits <- list_files_matching(dir, pattern = pattern, by = by, fixed = fixed,
                              ignore_case = ignore_case, recursive = recursive,
                              full_names = TRUE)
  if (length(hits) == 0L) {
    if (isTRUE(error)) {
      stop("No file matching ", if (is.null(pattern)) "*" else paste0("'", pattern, "'"),
           " in ", dir, call. = FALSE)
    }
    return(NA_character_)
  }
  hits[[1L]]
}
