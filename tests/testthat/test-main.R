test_that("Code Ocean panel uses named parameters accepted by main.R", {
  main_args <- extract_main_arguments(read_repo_file("code", "main.R"))
  panel_lines <- read_repo_file(".codeocean", "app-panel.json")
  panel_args <- extract_panel_param_names(panel_lines)

  expect_true(
    any(grepl('"named_parameters"[[:space:]]*:[[:space:]]*true', panel_lines)),
    info = "Code Ocean should pass parameters by name to main.R"
  )
  expect_same_values(panel_args, main_args)
})

test_that("summary volcano capsule keeps expected CLI parameter contract", {
  main_lines <- read_repo_file("code", "main.R")
  main_text <- paste(main_lines, collapse = "\n")

  # plot_volcano_summary() args the capsule intentionally does not expose as
  # their own CLI flags:
  #   - add_features is derived from whether custom_gene_list is non-empty
  #   - flip_vplot is hardcoded to FALSE
  #   - graphics_device/print_plots/save_plots/plots_subdir are capsule internals
  not_exposed_by_design <- c(
    "add_features",
    "flip_vplot",
    "graphics_device",
    "print_plots",
    "save_plots",
    "plots_subdir"
  )

  expected_args <- setdiff(
    mosuite_plot_volcano_summary_args(),
    not_exposed_by_design
  )

  expect_same_values(extract_main_arguments(main_lines), expected_args)
  expect_match(main_text, "plot_volcano_summary\\(")
  expect_match(
    main_text,
    "signif_colname = parse_optional_vector\\(args\\$signif_colname\\)"
  )
  expect_match(
    main_text,
    "change_colname = parse_optional_vector\\(args\\$change_colname\\)"
  )
  expect_match(main_text, "label_font_size = args\\$label_font_size")
  expect_match(main_text, "draw_connectors = args\\$draw_connectors")
  expect_match(
    main_text,
    "has_custom_gene_list <- nchar\\(trimws\\(args\\$custom_gene_list\\)\\) > 0"
  )
  expect_match(main_text, "add_features = has_custom_gene_list")
  expect_false("add_features" %in% extract_main_arguments(main_lines))
  expect_false("flip_vplot" %in% extract_main_arguments(main_lines))
  expect_match(main_text, "flip_vplot = FALSE")
  expect_match(main_text, "axis_lab_size = args\\$axis_lab_size")
  expect_match(main_text, "axis_tick_lab_size = args\\$axis_tick_lab_size")
})

test_that("Code Ocean panel preserves summary volcano defaults", {
  panel_lines <- read_repo_file(".codeocean", "app-panel.json")

  expect_true(is.na(extract_panel_default(panel_lines, "signif_colname")))
  expect_true(is.na(extract_panel_default(panel_lines, "change_colname")))
  expect_equal(
    extract_panel_default(panel_lines, "value_to_sort_the_output_dataset"),
    "t-statistic"
  )
  expect_equal(
    extract_panel_default(panel_lines, "num_features_to_label"),
    "20"
  )
  expect_equal(extract_panel_default(panel_lines, "label_font_size"), "5")
  expect_equal(
    extract_panel_default(panel_lines, "custom_label_color"),
    "black"
  )
  expect_equal(
    extract_panel_default(panel_lines, "color_of_non_significant_features"),
    "grey30"
  )
  expect_equal(
    extract_panel_default(
      panel_lines,
      "color_of_logfold_change_threshold_line"
    ),
    "forestgreen"
  )
  expect_equal(
    extract_panel_default(
      panel_lines,
      "color_of_features_meeting_only_signif_threshold"
    ),
    "royalblue"
  )
  expect_equal(
    extract_panel_default(
      panel_lines,
      "color_for_features_meeting_pvalue_and_foldchange_thresholds"
    ),
    "red2"
  )
  expect_equal(extract_panel_default(panel_lines, "axis_lab_size"), "24")
  expect_equal(extract_panel_default(panel_lines, "axis_tick_lab_size"), "16")
  expect_equal(extract_panel_default(panel_lines, "title"), "Volcano Plots")
  expect_equal(
    extract_panel_default(panel_lines, "plot_filename"),
    "volcano_summary.png"
  )
})

test_that("Code Ocean boolean controls are TRUE/FALSE lists", {
  panel_lines <- read_repo_file(".codeocean", "app-panel.json")

  expect_boolean_list_parameter(panel_lines, "label_features", "FALSE")
  expect_boolean_list_parameter(
    panel_lines,
    "label_significant_features_only",
    "TRUE"
  )
  expect_boolean_list_parameter(panel_lines, "draw_connectors", "FALSE")
  expect_boolean_list_parameter(panel_lines, "use_default_x_axis_limit", "TRUE")
  expect_boolean_list_parameter(panel_lines, "use_default_y_axis_limit", "TRUE")
  expect_boolean_list_parameter(panel_lines, "use_custom_lab", "FALSE")
  expect_boolean_list_parameter(panel_lines, "use_default_grid_layout", "TRUE")
})

test_that("main.R CLI creates volcano summary outputs", {
  setup <- setup_cli_workspace("mosuite_plot_volcano_summary_test_")
  on.exit(unlink(setup$workspace, recursive = TRUE), add = TRUE)

  old_wd <- getwd()
  setwd(setup$code_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system2("Rscript", args = c("main.R", common_cli_args))
  expect_equal(exit_code, 0, info = "main.R should execute without error")

  expect_outputs_created(setup$results_dir)
})

test_that("run wrapper executes and creates volcano summary outputs", {
  setup <- setup_cli_workspace("mosuite_plot_volcano_summary_run_test_")
  on.exit(unlink(setup$workspace, recursive = TRUE), add = TRUE)

  file.copy(
    file.path(setup$repo_root, "code", "run"),
    file.path(setup$code_dir, "run"),
    overwrite = TRUE
  )

  old_wd <- getwd()
  setwd(setup$code_dir)
  on.exit(setwd(old_wd), add = TRUE)

  exit_code <- system2("bash", args = c("run", common_cli_args))
  expect_equal(exit_code, 0, info = "run script should execute without error")

  expect_outputs_created(setup$results_dir)
})
