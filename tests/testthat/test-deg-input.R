deg_repo_root <- normalizePath(file.path(testthat::test_path(), "..", ".."))
suppressMessages(devtools::load_all(
  file.path(deg_repo_root, "code", "MOSuite"), quiet = TRUE
))
source(file.path(deg_repo_root, "code", "deg-input.R"))

deg_fixture <- function() {
  dat <- data.frame(GeneName = paste0("gene", seq_len(12)), check.names = FALSE)
  for (contrast in c("B_X-A_X", "B_X-C_X")) {
    logfc <- if (contrast == "B_X-A_X") {
      c(2, -2, 1, 3, 0, -3, rep(0.2, 6))
    } else {
      c(0.2, -2, 2, 0.2, 0, -3, rep(0.2, 6))
    }
    p <- c(0.01, 0.02, 0.01, 0.05, 0.2, 0.01, rep(0.5, 6))
    dat[[paste0(contrast, "_FC")]] <- ifelse(logfc < 0, -2^abs(logfc), 2^logfc)
    dat[[paste0(contrast, "_logFC")]] <- logfc
    dat[[paste0(contrast, "_SE")]] <- 0.25
    dat[[paste0(contrast, "_tstat")]] <- logfc / 0.25
    dat[[paste0(contrast, "_pval")]] <- p
    dat[[paste0(contrast, "_adjpval")]] <- p
  }
  dat[["A_X_Mean"]] <- 7
  dat[["A_X_SE"]] <- 0.3
  dat[["sample1"]] <- 7.5
  dat
}

test_that("DEG validation preserves OMIX statistics, names and expression", {
  dat <- deg_fixture()
  expect_identical(deg_validate_table(dat), dat)
  reordered <- dat[, c("sample1", setdiff(names(dat), "sample1"))]
  expect_identical(
    names(deg_validate_table(reordered, "GeneName"))[1], "GeneName"
  )
  bad <- dat
  bad$GeneName[1] <- bad$GeneName[2]
  expect_error(deg_validate_table(bad), "unique")
  bad <- dat
  bad[["B_X-A_X_adjpval"]][1] <- NA
  expect_error(deg_validate_table(bad), "finite")
  bad <- dat
  bad[["B_X-A_X_adjpval"]][1] <- 1.1
  expect_error(deg_validate_table(bad), "between 0 and 1")
  bad <- dat
  bad[["B_X-A_X_tstat"]] <- "not numeric"
  expect_error(deg_validate_table(bad), "numeric")
  expect_error(deg_validate_table(dat[, c("GeneName", "sample1")]), "_logFC")
  expect_error(
    deg_validate_table(dat[, c("GeneName", "B_X-A_X_logFC")]), "_pval"
  )
  expect_error(deg_validate_table(dat, "missing"), "not found")
  bad <- dat
  names(bad)[2] <- "GeneName"
  expect_error(deg_validate_table(bad), "duplicate column")
})

test_that("discovery uses MOO first and ignores OMIX companion CSV files", {
  root <- tempfile("deg-discovery-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE), add = TRUE)
  dat <- deg_fixture()
  write.csv(dat, file.path(root, "DEG_Analysis.csv"), row.names = FALSE)
  write.csv(
    data.frame(Sample = "sample1"), file.path(root, "Sample_Metadata.csv"),
    row.names = FALSE
  )
  write.csv(
    data.frame(Package = "limma"), file.path(root, "r-packages.csv"),
    row.names = FALSE
  )
  expect_equal(load_deg_input(data_dir = root), dat)
  moo <- deg_table_to_moo(dat)
  saveRDS(moo, file.path(root, "moo-diff.rds"))
  expect_s7_class(load_deg_input(data_dir = root), multiOmicDataSet)
  expect_equal(load_deg_input("table", data_dir = root), dat)
  saveRDS(moo, file.path(root, "other.rds"))
  expect_error(load_deg_input(data_dir = root), "exactly one")
  expect_s7_class(
    load_deg_input(input_file = file.path(root, "moo-diff.rds")),
    multiOmicDataSet
  )
  dir.create(file.path(root, "other"))
  write.csv(
    dat, file.path(root, "other", "DEG_Analysis.csv"), row.names = FALSE
  )
  expect_error(load_deg_input("table", data_dir = root), "exactly one")
  expect_error(
    load_deg_input(input_file = file.path(root, "missing.csv")),
    "does not exist"
  )
  saveRDS(list(not_a_moo = TRUE), file.path(root, "invalid.rds"))
  expect_error(
    load_deg_input(input_file = file.path(root, "invalid.rds")),
    "not a multiOmicDataSet"
  )
  tsv <- file.path(root, "selected.tsv")
  write.table(dat, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  expect_equal(load_deg_input("table", tsv), dat)
  expect_error(load_deg_input("moo", tsv), "RDS")
})

test_that("analysis-only MOO invents no counts or standard-deviation aliases", {
  dat <- deg_fixture()
  moo <- deg_table_to_moo(dat)
  expect_s7_class(moo, multiOmicDataSet)
  expect_silent(S7::validate(moo))
  expect_equal(nrow(moo@counts$raw), 0)
  expect_equal(nrow(moo@sample_meta), 0)
  expect_false(moo@analyses$deg_input$raw_counts_available)
  expect_identical(names(moo@analyses$diff), c("B_X-A_X", "B_X-C_X"))
  for (contrast in names(moo@analyses$diff)) {
    expect_identical(
      moo@analyses$diff[[contrast]]$SE, dat[[paste0(contrast, "_SE")]]
    )
    expect_false("sd" %in% names(moo@analyses$diff[[contrast]]))
  }
})

test_that("MOO loading bridges classes and rejects malformed contrast lists", {
  dat <- deg_fixture()
  moo <- deg_table_to_moo(dat)
  path <- tempfile(fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  bad <- moo
  bad@analyses$diff <- list(contrast = "not a table")
  saveRDS(bad, path)
  expect_error(load_deg_input(input_file = path), "per-contrast")
  bad <- moo
  names(bad@analyses$diff) <- c("", "valid")
  saveRDS(bad, path)
  expect_error(load_deg_input(input_file = path), "per-contrast")
  legacy_class <- S7::new_class(
    "multiOmicDataSet",
    package = "MOSuite",
    properties = list(
      sample_meta = S7::class_data.frame,
      annotation = S7::class_data.frame,
      counts = S7::class_list,
      analyses = S7::class_list
    )
  )
  legacy <- legacy_class(
    sample_meta = moo@sample_meta, annotation = moo@annotation,
    counts = moo@counts, analyses = moo@analyses
  )
  saveRDS(legacy, path)
  restored <- load_deg_input(input_file = path)
  expect_s7_class(restored, multiOmicDataSet)
  expect_equal(restored@analyses, moo@analyses)
  expect_equal(restored@counts, moo@counts)
  expect_equal(restored@sample_meta, moo@sample_meta)
  expect_equal(restored@annotation, moo@annotation)
  if (requireNamespace("MOObject", quietly = TRUE)) {
    current <- MOObject::multiOmicDataSet(
      sample_metadata = moo@sample_meta, anno_dat = moo@annotation,
      counts_lst = moo@counts, analyses_lst = moo@analyses
    )
    saveRDS(current, path)
    restored <- load_deg_input(input_file = path)
    expect_s7_class(restored, multiOmicDataSet)
    expect_equal(restored@analyses, moo@analyses)
    expect_equal(restored@counts, moo@counts)
    expect_equal(restored@sample_meta, moo@sample_meta)
    expect_equal(restored@annotation, moo@annotation)
  }
})

deg_run_cli <- function(input, mode, extra = character()) {
  workspace <- tempfile("deg-cli-")
  dir.create(workspace)
  code <- file.path(workspace, "code")
  data <- file.path(workspace, "data")
  results <- file.path(workspace, "results")
  dir.create(code)
  dir.create(data)
  dir.create(results)
  dir.create(file.path(results, "figures"))
  dir.create(file.path(results, "moo"))
  file.copy(file.path(deg_repo_root, "code", "deg-input.R"), code)
  main <- readLines(file.path(deg_repo_root, "code", "main.R"))
  main <- gsub(
    "devtools::load_all('/code/MOSuite')",
    sprintf(
      "devtools::load_all('%s')", file.path(deg_repo_root, "code", "MOSuite")
    ),
    main,
    fixed = TRUE
  )
  if (basename(deg_repo_root) == "MOSuite-plot-volcano-enhanced") {
    # Enhanced returns plot data but does not normally emit a companion CSV.
    main <- c(main, paste0(
      "write.csv(volc_enh_dat, ",
      "\"../results/moo/volcano_enhanced_data.csv\", row.names = FALSE)"
    ))
  }
  writeLines(main, file.path(code, "main.R"))
  if (mode == "moo") {
    saveRDS(input, file.path(data, "moo-diff.rds"))
  } else {
    write.csv(input, file.path(data, "DEG_Analysis.csv"), row.names = FALSE)
    write.csv(
      data.frame(Sample = "sample1"), file.path(data, "Sample_Metadata.csv"),
      row.names = FALSE
    )
  }
  old <- setwd(code)
  on.exit(setwd(old), add = TRUE)
  log <- suppressWarnings(system2(
    file.path(R.home("bin"), "Rscript"),
    c("main.R", extra),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(log, "status")
  if (!is.null(status) && status != 0L) {
    unlink(workspace, recursive = TRUE)
    writeLines(log)
    stop("Capsule CLI failed with status ", status)
  }
  list(workspace = workspace, results = results)
}

test_that("table and MOO CLI routes produce the same downstream artifacts", {
  repo_name <- basename(deg_repo_root)
  dat <- deg_fixture()
  moo <- deg_table_to_moo(dat)
  if (repo_name == "MOSuite-filter-diff") {
    extra <- "--round_estimates=FALSE"
    output <- "DEG_Analysis_Filtered.csv"
    image <- file.path("diff", "filt", "barchart.png")
  } else if (repo_name == "MOSuite-plot-venn-diagram") {
    extra <- c("--image_width=800", "--image_height=600")
    output <- "venn_diagram_data.csv"
    image <- file.path("diff", "venn_diagram.png")
  } else {
    extra <- c(
      "--num_features_to_label=0", "--image_width=3", "--image_height=3",
      "--dpi=72"
    )
    stem <- if (repo_name == "MOSuite-plot-volcano-summary") {
      "volcano_summary"
    } else {
      "volcano_enhanced"
    }
    output <- paste0(stem, "_data.csv")
    image <- if (stem == "volcano_enhanced") {
      file.path("diff", "volcano_enhanced_B_X-A_X.png")
    } else {
      file.path("diff", paste0(stem, ".png"))
    }
  }
  table_run <- deg_run_cli(dat, "table", extra)
  on.exit(unlink(table_run$workspace, recursive = TRUE), add = TRUE)
  moo_run <- deg_run_cli(moo, "moo", extra)
  on.exit(unlink(moo_run$workspace, recursive = TRUE), add = TRUE)
  read_output <- function(run) {
    read.csv(file.path(run$results, "moo", output), check.names = FALSE)
  }
  table_out <- read_output(table_run)
  moo_out <- read_output(moo_run)
  if (repo_name == "MOSuite-plot-volcano-enhanced") {
    # The existing Enhanced return retains unused source columns, including
    # OMIX expression. Compare all MOO statistics and computed rank columns.
    expect_true(all(names(moo_out) %in% names(table_out)))
    expect_equal(table_out[, names(moo_out), drop = FALSE], moo_out)
  } else {
    expect_equal(table_out, moo_out)
  }
  if (repo_name == "MOSuite-plot-volcano-summary") {
    expect_setequal(table_out$Contrast, c("B_X-A_X", "B_X-C_X"))
  }
  for (run in list(table_run, moo_run)) {
    expect_gt(file.info(file.path(run$results, "figures", image))$size, 0)
    expect_true(file.exists(file.path(run$results, "r-packages.csv")))
  }
  if (repo_name == "MOSuite-filter-diff") {
    expect_identical(table_out$GeneName, c("gene1", "gene2", "gene3", "gene6"))
    filtered <- readRDS(file.path(
      table_run$results, "moo", "moo-diff-filt.rds"
    ))
    expect_equal(filtered@analyses$diff_filt, table_out, ignore_attr = TRUE)
    all_run <- deg_run_cli(dat, "table", c(extra, "--filtering_mode=all"))
    on.exit(unlink(all_run$workspace, recursive = TRUE), add = TRUE)
    expect_identical(read_output(all_run)$GeneName, c("gene2", "gene6"))
  }
})

test_that("single-contrast tables and MOO inputs produce the same results", {
  repo <- basename(deg_repo_root)
  dat <- deg_fixture()
  dat <- dat[, !startsWith(names(dat), "B_X-C_X_"), drop = FALSE]
  extra <- if (repo == "MOSuite-filter-diff") {
    "--round_estimates=FALSE"
  } else if (repo == "MOSuite-plot-venn-diagram") {
    c("--image_width=800", "--image_height=600")
  } else {
    c(
      "--num_features_to_label=0", "--image_width=3", "--image_height=3",
      "--dpi=72"
    )
  }
  output <- switch(repo,
    "MOSuite-filter-diff" = "DEG_Analysis_Filtered.csv",
    "MOSuite-plot-venn-diagram" = "venn_diagram_data.csv",
    "MOSuite-plot-volcano-summary" = "volcano_summary_data.csv",
    "MOSuite-plot-volcano-enhanced" = "volcano_enhanced_data.csv"
  )
  table_run <- deg_run_cli(dat, "table", extra)
  on.exit(unlink(table_run$workspace, recursive = TRUE), add = TRUE)
  moo <- deg_table_to_moo(dat)
  analyses <- moo@analyses
  analyses$unrelated <- list(preserved = TRUE)
  moo <- multiOmicDataSet(
    sample_metadata = data.frame(Sample = "sample1", Group = "A_X"),
    anno_dat = moo@annotation,
    counts_lst = list(raw = data.frame(
      GeneName = dat$GeneName, sample1 = seq_len(nrow(dat))
    )),
    analyses_lst = analyses
  )
  moo_run <- deg_run_cli(moo, "moo", extra)
  on.exit(unlink(moo_run$workspace, recursive = TRUE), add = TRUE)
  read_output <- function(run) {
    read.csv(file.path(run$results, "moo", output), check.names = FALSE)
  }
  result <- read_output(table_run)
  expected <- read_output(moo_run)
  expect_equal(result[, names(expected), drop = FALSE], expected)
  if (repo == "MOSuite-filter-diff") {
    expect_identical(result$GeneName, c("gene1", "gene2", "gene6"))
    restored <- readRDS(file.path(moo_run$results, "moo", "moo-diff-filt.rds"))
    expect_equal(restored@sample_meta, moo@sample_meta)
    expect_equal(restored@counts, moo@counts)
    expect_equal(restored@annotation, moo@annotation)
    expect_equal(restored@analyses$unrelated, moo@analyses$unrelated)
  }
})

if (basename(deg_repo_root) == "MOSuite-plot-venn-diagram") {
  test_that("Venn CLI supports identical memberships for both inputs", {
    dat <- deg_fixture()
    for (contrast in c("B_X-A_X", "B_X-C_X")) {
      dat[[paste0(contrast, "_logFC")]] <- 2
      dat[[paste0(contrast, "_adjpval")]] <- 0.01
    }
    extra <- c("--image_width=800", "--image_height=600")
    table_run <- deg_run_cli(dat, "table", extra)
    on.exit(unlink(table_run$workspace, recursive = TRUE), add = TRUE)
    moo_run <- deg_run_cli(deg_table_to_moo(dat), "moo", extra)
    on.exit(unlink(moo_run$workspace, recursive = TRUE), add = TRUE)
    read_output <- function(run) {
      read.csv(file.path(run$results, "moo", "venn_diagram_data.csv"))
    }
    result <- read_output(table_run)
    expect_equal(result, read_output(moo_run))
    expect_setequal(result$Gene, dat$GeneName)
    expect_equal(unique(result$Id), 1)
    expect_equal(unique(result$Size), nrow(dat))
    for (run in list(table_run, moo_run)) {
      image <- file.path(
        run$results, "figures", "diff", "venn_diagram.png"
      )
      expect_gt(file.info(image)$size, 0)
    }
  })
}
