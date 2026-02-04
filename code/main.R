#!/usr/bin/env Rscript
rlang::global_entrace()
library(argparse)
library(glue)
library(MOSuite)
library(readr)
library(stringr)
library(dplyr)

# set up results directory
results_dir <- file.path('..','results')
plots_dir <- file.path(results_dir, 'figures')
options(moo_plots_dir = plots_dir, moo_save_plots = TRUE)

# log installed packages & versions
pkg_versions <- tibble::as_tibble(installed.packages())
write_csv(pkg_versions, file.path(results_dir, 'r-packages.csv'))

# parse CLI arguments
parser <- ArgumentParser()

parser$add_argument("--regex_moo", type="character", default=".*\\.rds$")
parser$add_argument("--feature_id_colname", type="character", default=NULL, help="Column name for feature IDs")
parser$add_argument("--signif_colname", type="character", default="pval", help="Column name of significance values")
parser$add_argument("--signif_threshold", type="double", default=0.05, help="Significance cutoff for p-values")
parser$add_argument("--change_threshold", type="double", default=1.0, help="Fold change cutoff for significance")
parser$add_argument("--value_to_sort_the_output_dataset", type="character", default="t-statistic", help="How to sort output: 'fold-change', 'p-value', or 't-statistic'")
parser$add_argument("--num_features_to_label", type="integer", default=30, help="Number of top features to label")
parser$add_argument("--add_features", type="logical", default=FALSE, help="Add custom_gene_list to labels")
parser$add_argument("--label_features", type="logical", default=FALSE, help="Label only features from custom_gene_list")
parser$add_argument("--custom_gene_list", type="character", default="", help="Comma-separated feature names to label")
parser$add_argument("--default_label_color", type="character", default="black", help="Color for default feature labels")
parser$add_argument("--custom_label_color", type="character", default="green3", help="Color for custom feature labels")
parser$add_argument("--label_x_adj", type="double", default=0.2, help="X-axis adjustment for labels")
parser$add_argument("--label_y_adj", type="double", default=0.2, help="Y-axis adjustment for labels")
parser$add_argument("--line_thickness", type="double", default=0.5, help="Thickness of lines in plot")
parser$add_argument("--label_font_size", type="double", default=4, help="Font size for labels")
parser$add_argument("--label_font_type", type="integer", default=1, help="Font type for labels")
parser$add_argument("--displace_feature_labels", type="logical", default=FALSE, help="Displace specific feature labels")
parser$add_argument("--custom_gene_list_special_label_displacement", type="character", default="", help="Features for special label displacement")
parser$add_argument("--special_label_displacement_x_axis", type="double", default=2, help="X displacement for special labels")
parser$add_argument("--special_label_displacement_y_axis", type="double", default=2, help="Y displacement for special labels")
parser$add_argument("--color_of_signif_threshold_line", type="character", default="blue", help="Color of significance threshold line")
parser$add_argument("--color_of_non_significant_features", type="character", default="black", help="Color of non-significant features")
parser$add_argument("--color_of_logfold_change_threshold_line", type="character", default="red", help="Color of fold change threshold line")
parser$add_argument("--color_of_features_meeting_only_signif_threshold", type="character", default="lightgoldenrod2", help="Color for features meeting only significance threshold")
parser$add_argument("--color_for_features_meeting_pvalue_and_foldchange_thresholds", type="character", default="red", help="Color for features meeting both thresholds")
parser$add_argument("--flip_vplot", type="logical", default=FALSE, help="Flip fold change values")
parser$add_argument("--use_default_x_axis_limit", type="logical", default=TRUE, help="Use default X-axis limit")
parser$add_argument("--x_axis_limit", type="double", default=5, help="Custom X-axis limit")
parser$add_argument("--use_default_y_axis_limit", type="logical", default=TRUE, help="Use default Y-axis limit")
parser$add_argument("--y_axis_limit", type="double", default=10, help="Custom Y-axis limit")
parser$add_argument("--point_size", type="double", default=2, help="Size of points in plot")
parser$add_argument("--add_deg_columns", type="character", default="FC,logFC,tstat,pval,adjpval", help="Additional DEG columns to include (comma-separated)")
parser$add_argument("--image_width", type="integer", default=15, help="Output image width")
parser$add_argument("--image_height", type="integer", default=15, help="Output image height")
parser$add_argument("--dpi", type="integer", default=300, help="Dots per inch of output image")
parser$add_argument("--use_default_grid_layout", type="logical", default=TRUE, help="Use default grid layout")
parser$add_argument("--number_of_rows_in_grid_layout", type="integer", default=1, help="Number of rows in grid layout")
parser$add_argument("--aspect_ratio", type="double", default=0, help="Aspect ratio of output")
parser$add_argument("--plot_filename", type="character", default="volcano_summary.png", help="Plot output filename")

args <- parser$parse_args()

parse_optional_vector <- function(x) {
    if (is.null(x) || identical(x, "") || length(x) == 0) {
        return(NULL)
    }
    return(trimws(unlist(strsplit(x, ","))))
}

parse_vector_with_default <- function(x, default) {
    parsed <- parse_optional_vector(x)
    if (is.null(parsed)) {
        return(default)
    }
    return(parsed)
}

# validate inputs
data_files <- list.files(file.path('../data'), recursive = TRUE, full.names = TRUE)
moo_files <- Filter(\(x) str_detect(x, regex(args$regex_moo, ignore_case = TRUE)), data_files)

if (length(moo_files) == 0) {
    stop(glue("No files matching regex: {args$regex_moo}"))
}
moo_filename <- moo_files[1]
moo <- read_rds(moo_filename)
message(glue('Reading multiOmicDataSet from {moo_filename}'))
if (!inherits(moo, 'MOSuite::multiOmicDataSet')) {
    stop(glue('The input is not a multiOmicDataSet. class: {class(moo)}'))
}

# run MOSuite
plot_volcano_summary(
    moo,
    feature_id_colname = args$feature_id_colname,
    signif_colname = args$signif_colname,
    signif_threshold = args$signif_threshold,
    change_threshold = args$change_threshold,
    value_to_sort_the_output_dataset = args$value_to_sort_the_output_dataset,
    num_features_to_label = args$num_features_to_label,
    add_features = args$add_features,
    label_features = args$label_features,
    custom_gene_list = args$custom_gene_list,
    default_label_color = args$default_label_color,
    custom_label_color = args$custom_label_color,
    label_x_adj = args$label_x_adj,
    label_y_adj = args$label_y_adj,
    line_thickness = args$line_thickness,
    label_font_size = args$label_font_size,
    label_font_type = args$label_font_type,
    displace_feature_labels = args$displace_feature_labels,
    custom_gene_list_special_label_displacement = args$custom_gene_list_special_label_displacement,
    special_label_displacement_x_axis = args$special_label_displacement_x_axis,
    special_label_displacement_y_axis = args$special_label_displacement_y_axis,
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
    add_deg_columns = parse_optional_vector(args$add_deg_columns),
    image_width = args$image_width,
    image_height = args$image_height,
    dpi = args$dpi,
    use_default_grid_layout = args$use_default_grid_layout,
    number_of_rows_in_grid_layout = args$number_of_rows_in_grid_layout,
    aspect_ratio = args$aspect_ratio,
    plot_filename = args$plot_filename
)
