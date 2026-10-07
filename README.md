# mito_meta_temp

How developmental **temperature** affects mitochondrial respiratory function: a systematic review and meta-analysis. Target journal: *Proceedings of the Royal Society B*.

## 1. Status (7 October 2026)

**The analysis dataset is not final — do not model it yet.** The data have been merged and checked, but the
corrections have not been applied and several inclusion decisions are open.

| Stage | State |
------- | -------
Merge the two extraction workbooks, subset temperature, clean | Done (`R/1_data_process.R`) |
Time-calibrated phylogeny | Done (TimeTree, 29 species) |
Automated checks on every row | Done (`output/checks/check_*.csv`) |
Read all 41 source papers against the extracted values | Done (September 2026) |
Apply the corrections to the extraction workbooks | **Not started** — the raw workbooks are unchanged |
Team decisions on what is eligible | **Open** — 12 decisions, see `data_checklist.xlsx` |
Rewrite the modelling code for one dataset | **Not started** — 8 `TODO` markers in `docs/results.qmd` |
References for the included studies | **Missing** from `bib/refs.bib` |

Of the 550 temperature rows, **523 need something**: 226 rows carry an error confirmed against the paper,
164 wait only on a decision, 87 need a check, and 46 have minor coding issues. Every critical and high-severity
finding has been independently confirmed; none was a false alarm. Start from
`output/checks/data_checklist.xlsx` (section 5).

### The dataset as it currently stands

| | |
-- | --
Temperature rows in the two workbooks | 550 (41 studies) |
Rows excluded | 5 flagged by the extractors, 37 non-mitochondrial measures, 2 with unusable statistics |
Analysis dataset | 507 effect sizes, 39 studies, 29 species, 6 classes |
By class | fish 215, birds 114, mammals 81, reptiles 66, amphibians 29, cartilaginous fish 2 |
By measure | metabolic capacity 215, antioxidant 141, respiration 86, oxidative damage 60, oxidative/nitrosative stress 5 |

A caveat that matters for interpretation: the studies span **egg incubation through to acclimation of
near-adult animals**. Roughly 113 effects come from embryonic manipulations, 162 from neonates, hatchlings
or larvae, and 232 from juveniles — some of those being short acclimations of well-grown animals. Life
stage at exposure is worth carrying as a moderator.

## 2. Project organisation

```
data/                              raw extraction workbooks (read-only; fixes are made here)
  mito_meta_data_merged_03072025.xlsx     extraction round 1
  Kris_mito_meta_data_for_temp_OC.xlsx    extraction round 2 (KHW)
R/
  func.R                           helper functions (tree checks, themes, heterogeneity tables)
  1_data_process.R                 merge -> subset temperature -> clean -> audit -> phylogeny -> effect sizes
docs/
  results.qmd                      Quarto source for the results and supplement
  results.docx                     rendered output (git-ignored)
bib/
  refs.bib                         references (background only so far — see Status)
  proceedings-of-the-royal-society-b.csl   citation style
  template.docx                    Word reference document
output/
  data/      temp_stress.csv       analysis dataset
             temp_stress_full.csv  every temperature row, with exclusion flags
  checks/    check_*.csv           audit reports written by the script (section 5)
             data_checklist.xlsx   the team's working checklist (section 5)
  phylo/                           trees and species lists
  tables/, fig_explore/            summary tables and exploratory figures
  figures/, models/                empty until the analysis is run
```

Run `R/1_data_process.R` first; everything in `docs/results.qmd` reads from `output/`.

## 3. Data pipeline

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
5. **Audit.** Write the check files described in section 5. Nothing is dropped silently.
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

## 4. Checking the data

Three layers, all in `output/checks/`.

### Written by the script, every run

Diagnostics, not exclusions — review them before trusting the analysis dataset.

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

### Paper cross-check (one-off, September 2026)

Every one of the 41 source papers was read against the extracted values. Not regenerated by the script.

| File | What it reports |
------ | -----------------
`check_paper_crosscheck.csv` | row-level problems: severity, affected rows, the evidence in the paper, the recommended fix, and whether the problem was independently confirmed
`check_paper_study_notes.csv` | per-study notes: the error statistic and sample sizes the paper reports, which group is the control, crossed second factors, and developmental timing

Row IDs in these two files are `merged_<n>` / `kris_n_<n>`, where `<n>` matches the script's `row_key`
(`merged_03072025_<n>` / `kris_temp_OC_<n>`).

### `data_checklist.xlsx` — start here

The team's working document, in four sheets:

- **Rows to check** — one line per row needing attention: problem codes, a plain explanation, the suggested
  change, the source workbook and row to edit, and yellow columns for who checked it and the outcome.
- **All data** — all 550 rows with every original column, plus each row's codes.
- **Key findings** — the confirmed errors by study in plain language, and the 12 team decisions.
- **Metadata** — how to use it, what each of the 30 codes means, and column definitions.

Counts update as the Outcome column is filled in. It is a snapshot (built 17 September 2026) and is **not**
regenerated by the script, so edits made in it are kept.

## 5. Reproducing the results

```r
source("R/1_data_process.R")      # writes output/data, output/phylo, output/checks, output/tables
quarto::quarto_render("docs/results.qmd")
```

Two things to know before rendering:

- `output/models/` is empty, so the first run has to fit the models: set `rerun <- TRUE` in
  `docs/results.qmd` (it appears at the top of each modelling chunk). With `rerun <- FALSE` the chunks
  expect cached model objects that do not exist yet.
- The modelling chunks were inherited from a project with several stressor datasets and still loop over a
  list. They run against a one-element list, but the summary tables and figure panels need rewriting for a
  single dataset — the 8 `TODO (temperature restructure)` markers show where.

## 6. Software

R 4.5.3, with `pacman::p_load()` installing what is missing. Main packages: `metafor` (5.0.1) for the
models, `orchaRd` (2.2.0) for orchard plots and heterogeneity, `rotl` and `ape`/`phytools` for the
phylogeny, `tidyverse`, `here`.

`orchaRd` is not on CRAN — install it first:

```r
remotes::install_github("daniel1noble/orchaRd", force = TRUE)
```

## 7. Source papers

PDFs are not kept in this repository. The round-1 studies are in the companion extraction folder
(`Papers/Included papers/`); the eight round-2 studies (`WildNew1`–`WildNew8`) are scattered across other
folders, and Tullis & Baillie (2005) is in the Zotero library. `check_paper_study_notes.csv` identifies each
paper by full title, authors, year and journal.
