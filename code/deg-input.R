# Capsule-owned adapter shared by the four downstream DEG capsules.
# Keep copies synchronized; do not place this file in the managed MOSuite tree.

deg_unique_file <- function(paths, label) {
  if (length(paths) != 1L) {
    stop(
      "Expected exactly one ",
      label,
      " input; found ",
      length(paths),
      if (length(paths)) paste0(": ", paste(paths, collapse = ", ")) else ".",
      " Remove extra matching inputs and attach only the intended result.",
      call. = FALSE
    )
  }
  paths[[1L]]
}

deg_validate_table <- function(dat, feature_id_colname = NULL) {
  if (!is.data.frame(dat) || ncol(dat) < 2L || nrow(dat) == 0L) {
    stop(
      "DEG table must contain features and contrast statistics.",
      call. = FALSE
    )
  }
  if (anyDuplicated(names(dat))) {
    stop("DEG table contains duplicate column names.", call. = FALSE)
  }
  if (is.null(feature_id_colname)) {
    feature_id_colname <- names(dat)[[1L]]
  }
  if (!feature_id_colname %in% names(dat)) {
    stop("Feature ID column not found: ", feature_id_colname, call. = FALSE)
  }
  ids <- as.character(dat[[feature_id_colname]])
  if (anyNA(ids) || any(!nzchar(trimws(ids))) || anyDuplicated(ids)) {
    stop("DEG feature IDs must be unique and non-empty.", call. = FALSE)
  }
  contrasts <- sub("_logFC$", "", grep("_logFC$", names(dat), value = TRUE))
  if (!length(contrasts) || any(!nzchar(contrasts))) {
    stop("DEG table requires contrast-prefixed _logFC columns.", call. = FALSE)
  }
  for (contrast in contrasts) {
    p_columns <- intersect(paste0(contrast, c("_pval", "_adjpval")), names(dat))
    if (!length(p_columns)) {
      stop("Missing _pval or _adjpval for contrast: ", contrast, call. = FALSE)
    }
    columns <- intersect(
      paste0(
        contrast,
        c("_FC", "_logFC", "_SE", "_tstat", "_pval", "_adjpval")
      ),
      names(dat)
    )
    for (column in columns) {
      values <- dat[[column]]
      if (!is.numeric(values) || any(!is.finite(values))) {
        stop(
          "DEG statistic must be finite and numeric: ",
          column,
          call. = FALSE
        )
      }
      if (column %in% p_columns && any(values < 0 | values > 1)) {
        stop("DEG p-values must be between 0 and 1: ", column, call. = FALSE)
      }
    }
  }
  # Put the chosen ID first: MOSuite's wide-table joins use the first column.
  dat <- dat[,
    c(
      feature_id_colname,
      setdiff(names(dat), feature_id_colname)
    ),
    drop = FALSE
  ]
  dat[[feature_id_colname]] <- ids
  dat
}

load_deg_input <- function(
  feature_id_colname = NULL,
  data_dir = file.path("..", "data")
) {
  files <- list.files(data_dir, recursive = TRUE, full.names = TRUE)
  rds <- files[grepl("[.]rds$", files, ignore.case = TRUE)]

  if (length(rds)) {
    path <- deg_unique_file(rds, "MOO RDS")
    input_type <- "moo"
  } else {
    tables <- files[tolower(basename(files)) == "deg_analysis.csv"]
    if (!length(tables)) {
      stop(
        "No supported DEG input found. Attach one MOSuite MOO RDS or a ",
        "CSV named DEG_Analysis.csv.",
        call. = FALSE
      )
    }
    path <- deg_unique_file(tables, "DEG_Analysis.csv")
    input_type <- "table"
  }

  message("Reading ", input_type, " DEG input from ", path)
  if (input_type == "moo") {
    moo <- readr::read_rds(path)
    if (
      !S7::S7_inherits(moo, S7::S7_object) ||
        !any(grepl("(^|::)multiOmicDataSet$", class(moo)))
    ) {
      stop("RDS input is not a multiOmicDataSet.", call. = FALSE)
    }
    # Capsule pins include both the older MOSuite and newer MOObject classes.
    if (!S7::S7_inherits(moo, MOSuite::multiOmicDataSet)) {
      moo <- MOSuite::multiOmicDataSet(
        sample_metadata = moo@sample_meta,
        anno_dat = moo@annotation,
        counts_lst = moo@counts,
        analyses_lst = moo@analyses
      )
    }
    S7::validate(moo)
    if (
      !is.list(moo@analyses$diff) ||
        !length(moo@analyses$diff) ||
        is.null(names(moo@analyses$diff)) ||
        anyNA(names(moo@analyses$diff)) ||
        any(!nzchar(trimws(names(moo@analyses$diff)))) ||
        anyDuplicated(names(moo@analyses$diff)) ||
        !all(vapply(moo@analyses$diff, is.data.frame, logical(1)))
    ) {
      stop(
        "MOO must contain named per-contrast analyses$diff tables.",
        call. = FALSE
      )
    }
    return(moo)
  }
  dat <- utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  deg_validate_table(dat, feature_id_colname)
}

deg_table_to_moo <- function(dat, feature_id_colname = NULL) {
  dat <- deg_validate_table(dat, feature_id_colname)
  id <- names(dat)[[1L]]
  contrasts <- sub("_logFC$", "", grep("_logFC$", names(dat), value = TRUE))
  statistics <- c("FC", "logFC", "SE", "tstat", "pval", "adjpval")
  diff <- lapply(contrasts, function(contrast) {
    columns <- intersect(paste0(contrast, "_", statistics), names(dat))
    result <- dat[, c(id, columns), drop = FALSE]
    names(result)[-1L] <- substring(columns, nchar(contrast) + 2L)
    result
  })
  names(diff) <- contrasts
  # This is a results container, not reconstructed expression/count data.
  MOSuite::multiOmicDataSet(
    sample_metadata = data.frame(Sample = character()),
    anno_dat = dat[, id, drop = FALSE],
    counts_lst = list(raw = dat[FALSE, id, drop = FALSE]),
    analyses_lst = list(
      diff = diff,
      deg_input = list(
        type = "table",
        raw_counts_available = FALSE,
        source_table = dat
      )
    )
  )
}
