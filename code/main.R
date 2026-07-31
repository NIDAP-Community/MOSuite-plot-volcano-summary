#!/usr/bin/env Rscript
rlang::global_entrace()
library(argparse)
library(glue)
library(readr)
library(stringr)
library(dplyr)
devtools::load_all('/code/MOSuite')

# set up capsule environment
setup_capsule_environment()

# parse CLI arguments
parser <- ArgumentParser()

parser$add_argument(
  "--feature_id_colname",
  type = "character",
  default = NULL,
  help = "Column name for feature IDs"
)
parser$add_argument(
  "--signif_colname",
  type = "character",
  default = "",
  help = "Column name of significance values; leave blank to auto-detect"
)
parser$add_argument(
  "--change_colname",
  type = "character",
  default = "",
  help = "Column name of fold change values; leave blank to auto-detect"
)
parser$add_argument(
  "--signif_threshold",
  type = "double",
  default = 0.05,
  help = "Significance cutoff for p-values"
)
parser$add_argument(
  "--change_threshold",
  type = "double",
  default = 1.0,
  help = "Fold change cutoff for significance"
)
parser$add_argument(
  "--value_to_sort_the_output_dataset",
  type = "character",
  default = "t-statistic",
  help = "How to sort output: 'fold-change', 'p-value', or 't-statistic'"
)
parser$add_argument(
  "--num_features_to_label",
  type = "integer",
  default = 20,
  help = "Number of top features to label"
)
parser$add_argument(
  "--label_features",
  type = "logical",
  default = FALSE,
  help = "Label only features from custom_gene_list; otherwise custom_gene_list labels are added to automatic labels"
)
parser$add_argument(
  "--custom_gene_list",
  type = "character",
  default = "",
  help = "Comma-separated feature names to label"
)
parser$add_argument(
  "--label_significant_features_only",
  type = "logical",
  default = TRUE,
  help = "Label only features that meet the significance threshold"
)
parser$add_argument(
  "--default_label_color",
  type = "character",
  default = "black",
  help = "Color for default feature labels"
)
parser$add_argument(
  "--custom_label_color",
  type = "character",
  default = "black",
  help = "Color for custom feature labels"
)
parser$add_argument(
  "--label_font_size",
  type = "double",
  default = 5,
  help = "Font size for labels"
)
parser$add_argument(
  "--draw_connectors",
  type = "logical",
  default = FALSE,
  help = "Draw connector lines from labels to their points"
)
parser$add_argument(
  "--change_sig_name",
  type = "character",
  default = "p-value",
  help = "Display name for significance values"
)
parser$add_argument(
  "--change_lfc_name",
  type = "character",
  default = "log2FC",
  help = "Display name for fold change values"
)
parser$add_argument(
  "--title",
  type = "character",
  default = "Volcano Plots",
  help = "Plot title"
)
parser$add_argument(
  "--title_font_size",
  type = "double",
  default = 24,
  help = "Font size for the plot title"
)
parser$add_argument(
  "--use_custom_lab",
  type = "logical",
  default = FALSE,
  help = "Use custom axis labels"
)
parser$add_argument(
  "--color_of_signif_threshold_line",
  type = "character",
  default = "black",
  help = "Color of significance threshold line"
)
parser$add_argument(
  "--color_of_non_significant_features",
  type = "character",
  default = "grey30",
  help = "Color of non-significant features"
)
parser$add_argument(
  "--color_of_logfold_change_threshold_line",
  type = "character",
  default = "forestgreen",
  help = "Color of fold change threshold line"
)
parser$add_argument(
  "--color_of_features_meeting_only_signif_threshold",
  type = "character",
  default = "royalblue",
  help = "Color for features meeting only significance threshold"
)
parser$add_argument(
  "--color_for_features_meeting_pvalue_and_foldchange_thresholds",
  type = "character",
  default = "red2",
  help = "Color for features meeting both thresholds"
)
parser$add_argument(
  "--flip_vplot",
  type = "logical",
  default = FALSE,
  help = "Flip fold change values"
)
parser$add_argument(
  "--use_default_x_axis_limit",
  type = "logical",
  default = TRUE,
  help = "Use default X-axis limit"
)
parser$add_argument(
  "--x_axis_limit",
  type = "double",
  default = 5,
  help = "Custom X-axis limit"
)
parser$add_argument(
  "--use_default_y_axis_limit",
  type = "logical",
  default = TRUE,
  help = "Use default Y-axis limit"
)
parser$add_argument(
  "--y_axis_limit",
  type = "double",
  default = 10,
  help = "Custom Y-axis limit"
)
parser$add_argument(
  "--point_size",
  type = "double",
  default = 2,
  help = "Size of points in plot"
)
parser$add_argument(
  "--axis_lab_size",
  type = "double",
  default = 24,
  help = "Size of axis labels"
)
parser$add_argument(
  "--axis_tick_lab_size",
  type = "double",
  default = 16,
  help = "Size of axis tick labels"
)
parser$add_argument(
  "--add_deg_columns",
  type = "character",
  default = "FC,logFC,tstat,pval,adjpval",
  help = "Additional DEG columns to include (comma-separated)"
)
parser$add_argument(
  "--image_width",
  type = "integer",
  default = 15,
  help = "Output image width"
)
parser$add_argument(
  "--image_height",
  type = "integer",
  default = 15,
  help = "Output image height"
)
parser$add_argument(
  "--dpi",
  type = "integer",
  default = 300,
  help = "Dots per inch of output image"
)
parser$add_argument(
  "--use_default_grid_layout",
  type = "logical",
  default = TRUE,
  help = "Use default grid layout"
)
parser$add_argument(
  "--number_of_rows_in_grid_layout",
  type = "integer",
  default = 1,
  help = "Number of rows in grid layout"
)
parser$add_argument(
  "--plot_filename",
  type = "character",
  default = "volcano_summary.png",
  help = "Plot output filename"
)

args <- parser$parse_args()

# load multiOmicDataSet from data directory
moo <- load_moo_from_data_dir()
has_custom_gene_list <- nchar(trimws(args$custom_gene_list)) > 0

# run MOSuite
summary_dat <- plot_volcano_summary(
  moo,
  feature_id_colname = args$feature_id_colname,
  change_colname = parse_optional_vector(args$change_colname),
  signif_colname = parse_optional_vector(args$signif_colname),
  signif_threshold = args$signif_threshold,
  change_threshold = args$change_threshold,
  value_to_sort_the_output_dataset = args$value_to_sort_the_output_dataset,
  num_features_to_label = args$num_features_to_label,
  add_features = has_custom_gene_list,
  label_features = args$label_features,
  custom_gene_list = args$custom_gene_list,
  label_significant_features_only = args$label_significant_features_only,
  default_label_color = args$default_label_color,
  custom_label_color = args$custom_label_color,
  label_font_size = args$label_font_size,
  draw_connectors = args$draw_connectors,
  change_sig_name = args$change_sig_name,
  change_lfc_name = args$change_lfc_name,
  title = args$title,
  title_font_size = args$title_font_size,
  use_custom_lab = args$use_custom_lab,
  color_of_signif_threshold_line = args$color_of_signif_threshold_line,
  color_of_non_significant_features = args$color_of_non_significant_features,
  color_of_logfold_change_threshold_line = args$color_of_logfold_change_threshold_line,
  color_of_features_meeting_only_signif_threshold = args$color_of_features_meeting_only_signif_threshold,
  color_for_features_meeting_pvalue_and_foldchange_thresholds = args$color_for_features_meeting_pvalue_and_foldchange_thresholds,
  flip_vplot = args$flip_vplot,
  use_default_x_axis_limit = args$use_default_x_axis_limit,
  x_axis_limit = args$x_axis_limit,
  use_default_y_axis_limit = args$use_default_y_axis_limit,
  y_axis_limit = args$y_axis_limit,
  point_size = args$point_size,
  axis_lab_size = args$axis_lab_size,
  axis_tick_lab_size = args$axis_tick_lab_size,
  add_deg_columns = parse_optional_vector(args$add_deg_columns),
  image_width = args$image_width,
  image_height = args$image_height,
  dpi = args$dpi,
  use_default_grid_layout = args$use_default_grid_layout,
  number_of_rows_in_grid_layout = args$number_of_rows_in_grid_layout,
  plot_filename = args$plot_filename
)

readr::write_csv(
  summary_dat,
  file.path(getOption("moo_plots_dir"), "..", "moo", "volcano_summary_data.csv")
)
