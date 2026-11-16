# MOSuite-plot-volcano-summary

Code Ocean capsule: MOSuite - Plot Volcano Summary

- [Code Ocean Capsule](https://poc-nci.codeocean.io/capsule/9325734/tree) | [Latest Release](https://poc-nci.codeocean.io/capsule/6849617/tree/latest)
- [MOSuite R package docs](https://ccbr.github.io/MOSuite/)

## DEG inputs

Attach either a MOSuite DEG MOO (`.rds`) or the wide `DEG_Analysis.csv`
emitted by [OMIX DEG Analysis](https://github.com/NIDAP-Community/OMIX/tree/main/modules/OMIX-DEG-Analysis).
No counts or sample metadata are required for these result plots.

- **Input type / `--input_type`**: `auto` (default), `moo`, or `table`.
  Auto prefers exactly one RDS; if none is attached, it discovers one
  `DEG_Analysis.csv`, `.tsv`, or `.txt` recursively. Companion metadata,
  package inventories, and count tables are not selected.
- **Input file / `--input_file`**: optional explicit file selector. Use this
  for a differently named table or ambiguous attached inputs. Select `table`
  to use the table when both MOO and CSV are attached.
- Default settings require unique feature IDs and contrast-prefixed `FC`,
  `logFC`, `tstat`, `pval`, and `adjpval` statistics (e.g. `B-A_logFC`). Feature IDs
  default to the first column (`GeneName` in OMIX); override with
  `--feature_id_colname`. Leave significance/change column controls blank for
  automatic detection, or enter exact columns for selected contrasts.

Before plotting, the capsule checks the columns needed by the selected
contrasts, sorting, and summary export settings. Reduced tables are supported:
for example, `GeneName`, `B-A_logFC`, and `B-A_adjpval` work with
`--value_to_sort_the_output_dataset=fold-change` and
`--add_deg_columns=logFC,adjpval`. Use `--add_deg_columns=none` to export only
feature IDs and contrast memberships. Missing required statistics are reported
explicitly; they are not reconstructed. `SE` is optional unless requested for
export.

The same MOSuite data-frame method processes both routes; MOO statistics are
joined into a wide table first. Existing plot and
`results/moo/volcano_summary_data.csv` paths and thresholds are unchanged.
OMIX's standard errors are not renamed to standard deviations.

## Validation and maintenance

Run `Rscript tests/testthat.R` from the repository root. Tests cover MOO/table
CLI parity, input ambiguity, OMIX companion files, column validation and
contrast names containing underscores. Code Ocean runs are still required
before release; local tests are not platform validation.

`code/deg-input.R` is capsule-owned and shared verbatim with the Enhanced
Volcano, Venn and Filter Differential Expression capsules. Keep the four adapter
copies and common input tests synchronized; Summary-specific sorting/export
preflight tests are maintained here. SyncWeaver-managed `code/MOSuite` and its
source pin are unchanged.