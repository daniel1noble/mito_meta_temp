# mito_meta_temp

How developmental **temperature** affects mitochondrial respiratory function: a systematic review and meta-analysis.

## 1. Project organisation

```
data/                              raw extraction workbooks (read-only)
  mito_meta_data_merged_03072025.xlsx     extraction round 1
  Kris_mito_meta_data_for_temp_OC.xlsx    extraction round 2 (KHW)
R/
  func.R                           helper functions (tree checks, themes, tables)
  1_data_process.R                 merge -> subset temperature -> clean -> phylogeny -> effect sizes
docs/
  results.qmd                      Quarto source for the results and supplement
bib/
  refs.bib                         references
  proceedings-of-the-royal-society-b.csl   citation style (target journal: Proc. R. Soc. B)
  template.docx                    Word reference document
output/
  data/      temp_stress.csv       analysis dataset
             temp_stress_full.csv  every temperature row, with exclusion flags
  checks/    check_*.csv           data-audit reports (see section 3)
  phylo/                           trees and species lists
  figures/, fig_explore/, tables/, models/
```

Run `R/1_data_process.R` first; everything in `docs/results.qmd` reads from `output/`.

## 2. Data pipeline

`R/1_data_process.R` does the following.

1. **Merge.** Both workbooks use the same 78-column extraction template, so they are row-bound after
   asserting that the column names match and that no `study` ID is shared between them. Every column is read
   as text and coerced explicitly, because `readxl` otherwise infers different types for the same column in
   the two files, and because non-numeric entries in numeric columns (which are informative) would otherwise
   be silently turned into `NA`.
2. **Subset.** Keep `envirn_type == "temp"`.
3. **Harmonise.** Squish whitespace in all character keys (`"s58 "` and `"s58"` were separate levels),
   collapse synonymous tissue labels, recode `stage`, and fix species binomials so they match the Open Tree
   of Life (several were misspelled; *Dicentrarchus labrax* appeared under two different misspellings).
4. **Derive temperature moderators.** `t1` is the treatment temperature and `t2` the control, so
   `temp_diff = t1 - t2`, `temp_direction` (warmer / cooler) and `temp_magnitude = |temp_diff|`.
   `temp_regime` records the one heat-wave design (it had been recorded in the `nutrition_sum` column).
5. **Audit.** Write the check files described below. Nothing is dropped silently.
6. **Exclude.** Add three flag columns and filter on them:
   - `flag_missing_info` — the extractors' `exclude_missing_info == 0`, meaning the error statistic or
     sample size was missing or ambiguous.
   - `flag_non_mito` — measures that are not mitochondrial traits: lactate dehydrogenase (cytosolic,
     anaerobic), Na<sup>+</sup>/K<sup>+</sup>-ATPase, alkaline phosphatase, and plasma corticosterone.
   - `flag_bad_stats` — an effect size cannot be computed (zero SD, missing mean/SD/N, N < 2).
7. **Phylogeny.** The models use a time-calibrated tree from [TimeTree 5](https://timetree.org):
   `output/phylo/phylo_pruned_species.nwk` (29 species, ultrametric, root 462.4 MYA; downloaded
   16 September 2026). The script does not build this tree. It matches species names against the Open Tree
   of Life with `rotl`, writes the matched list to `output/phylo/phylo_pruned_species.txt`, and stops if
   the species in the data no longer match the tips of the TimeTree file. Whenever the species list
   changes, resubmit that list to TimeTree ("Load a List of Species", then export to Newick) and replace the
   `.nwk`. TimeTree has no data for *Chiloscyllium plagiosum* and substitutes its congener
   *C. punctatum*; the tip keeps the *plagiosum* label. The `rotl` tree (`phylo.tre`, `phylo_pruned.tre`,
   plotted in `phylo.png`) has no branch lengths and is only a name check — it is not used in the models.
8. **Effect sizes.** `SMD_H` (Bonett 2008, 2009) via `metafor::escalc`, then correct the sign so that a
   positive effect always means *higher* mitochondrial respiratory function.

### Sign convention

`mito_efficiency_dir` is `1` when an increase in the measured variable means higher mitochondrial
respiratory function and `0` when an increase means lower function; effect sizes in the latter group are
multiplied by −1. Antioxidants, ATP-linked respiration, RCR, coupling efficiency, citrate synthase and
cytochrome *c* oxidase are coded `1`; ROS, malondialdehyde/TBARS, protein carbonyls, other oxidative damage,
proton leak / state-4 respiration and uncoupling protein expression are coded `0`. The script asserts that
only `0` and `1` occur before flipping signs, rather than letting an unrecognised code pass through unflipped.

## 3. Data-audit reports

`R/1_data_process.R` writes one CSV per check to `output/checks/`. These are diagnostics, not
exclusions — review them before trusting the analysis dataset.

| File | What it reports |
------ | -----------------
`check_missing_stats.csv` | rows where an effect size cannot be computed (missing or non-positive SD, N < 2)
`check_sd_se.csv` | rows where `sd != se * sqrt(n)`, with the implied *n*. `sd == se` means a standard error was entered in the SD column without conversion
`check_dir_codes.csv` | `mito_efficiency_dir` values other than 0 or 1
`check_dir_inconsistent.csv` | the same measure within one study coded with both direction values
`check_cat_inconsistent.csv` | the same measure within one study assigned to two different `measurement_category` values
`check_shared_controls.csv` | control groups compared against more than one treatment group. Kept in the analysis (the correlation is absorbed by the study and sample-dependency random effects) but recorded, because a shared control that is *not* actually shared in the paper indicates a mis-extraction
`check_suspicious.csv` | identical means or identical SDs between the two groups — implausible for values digitised from figures
`check_flagged_exclude.csv` | rows the extractors flagged for exclusion, plus any non-0/1 entry in that column
`check_extreme_effects.csv` | \|SMD<sub>H</sub>\| > 5, for the sensitivity analysis
`check_paper_crosscheck.csv` | row-level problems found by re-reading every source paper against the extracted values (severity, affected `row_id`s, the paper evidence, the recommended fix, and whether the problem was independently confirmed). Not generated by the script — written by a one-off audit (September 2026)
`check_paper_study_notes.csv` | per-study notes from that audit: error statistic and sample sizes the paper reports, which group is the control, crossed second factors, and developmental timing

The `row_id`s in the two `check_paper_*` files are `merged_<n>` / `kris_n_<n>`, where `<n>` is the row's position in the combined temperature data (merged-workbook temperature rows first, in workbook order, then the round-2 rows). The number is the same as in the script's `row_key` (`merged_03072025_<n>` / `kris_temp_OC_<n>`).

## 4. Reproducing the results

```r
source("R/1_data_process.R")      # writes output/data, output/phylo, output/checks
quarto::quarto_render("docs/results.qmd")
```

Model fitting in `docs/results.qmd` is gated behind `rerun = FALSE`; set it to `TRUE` to refit rather than
read the cached objects in `output/models/`.
