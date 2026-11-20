# Code Ocean capsule - MOSuite - Plot Volcano Summary

## development version

- Accept MOSuite DEG MOO or the specifically named `DEG_Analysis.csv`, with
  automatic discovery that prefers MOO when both supported formats are
  present. (#15, @TJoshMeyer)
- Validate tabular DEG inputs and preserve existing plotting/output behavior.
- Check required statistics against selected contrasts, sorting, and export
  settings before plotting; allow reduced tables with compatible options and
  keep standard errors optional unless selected for export. (#15, @TJoshMeyer)
- Validate selected table columns and support mixed significance suffixes
  across contrasts when no common suffix is available. (#15, @TJoshMeyer)
- Preserve numeric-looking feature identifiers when reading DEG CSV files.
  (#15, @TJoshMeyer)

## v6.0

- Use MOSuite v0.5.0. (#10, @kelly-sovacool)
- Use physical inch dimensions with 10 x 10 inch defaults and 300 DPI for saved volcano output. (#9, @TJoshMeyer)

## v5.0

- Disable the nonfunctional flip volcano plot app panel option for MOSuite v0.4.1 (#6, @phoman14).
- Improved the Code Ocean parameter UI for the plot volcano summary capsule (#2, @phoman14).
- Use MOSuite v0.3.2. (#3, @kelly-sovacool)

## v4.0

- The MOSuite package is now available in `code/MOSuite`. (#1)
- Use MOSuite v0.3.1.

## v3.0

- Use MOSuite v0.3.0.
- Fix: remove `regex_moo` parameter from configuration. multiOmicDataSet input files are now required to follow the standardized naming pattern (`.*\.rds$`).

## v2.0

Save summarized dataframe to `results/moo/volcano_summary_data.csv`.

<https://poc-nci.codeocean.io/capsule/6849617/tree/v2> `f64a38f`

## v1.0

Initial release

<https://poc-nci.codeocean.io/capsule/6849617/tree/v1> `0a5bebc`
