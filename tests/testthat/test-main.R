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

  expected_args <- c(
    "feature_id_colname",
    "signif_colname",
    "change_colname",
    "signif_threshold",
    "change_threshold",
    "value_to_sort_the_output_dataset",
    "num_features_to_label",
    "label_features",
    "custom_gene_list",
    "label_significant_features_only",
    "default_label_color",
    "custom_label_color",
    "label_font_size",
    "draw_connectors",
    "change_sig_name",
    "change_lfc_name",
    "title",
    "title_font_size",
    "use_custom_lab",
    "color_of_signif_threshold_line",
    "color_of_non_significant_features",
    "color_of_logfold_change_threshold_line",
    "color_of_features_meeting_only_signif_threshold",
    "color_for_features_meeting_pvalue_and_foldchange_thresholds",
    "use_default_x_axis_limit",
    "x_axis_limit",
    "use_default_y_axis_limit",
    "y_axis_limit",
    "point_size",
    "axis_lab_size",
    "axis_tick_lab_size",
    "add_deg_columns",
    "image_width",
    "image_height",
    "dpi",
    "use_default_grid_layout",
    "number_of_rows_in_grid_layout",
    "plot_filename"
  )

  removed_args <- c(
    "label_x_adj",
    "label_y_adj",
    "line_thickness",
    "label_font_type",
    "displace_feature_labels",
    "custom_gene_list_special_label_displacement",
    "special_label_displacement_x_axis",
    "special_label_displacement_y_axis",
    "aspect_ratio"
  )

  expect_same_values(extract_main_arguments(main_lines), expected_args)
  expect_false(any(removed_args %in% extract_main_arguments(main_lines)))
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
  expect_match(main_text, "has_custom_gene_list <- nchar\\(trimws\\(args\\$custom_gene_list\\)\\) > 0")
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
  expect_equal(extract_panel_default(panel_lines, "custom_label_color"), "black")
  expect_equal(
    extract_panel_default(panel_lines, "color_of_non_significant_features"),
    "grey30"
  )
  expect_equal(
    extract_panel_default(panel_lines, "color_of_logfold_change_threshold_line"),
    "forestgreen"
  )
  expect_equal(
    extract_panel_default(panel_lines, "color_of_features_meeting_only_signif_threshold"),
    "royalblue"
  )
  expect_equal(
    extract_panel_default(panel_lines, "color_for_features_meeting_pvalue_and_foldchange_thresholds"),
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
