#### --------------------------------------------------  ####
# 1. Merge data sources, subset developmental TEMPERATURE data,
#    clean, build phylogeny and calculate effect sizes.
#
#    Inputs   data/mito_meta_data_merged_03072025.xlsx   (extraction round 1)
#             data/Kris_mito_meta_data_for_temp_OC.xlsx  (extraction round 2, KHW)
#    Outputs  output/data/temp_stress.csv       analysis dataset
#             output/data/temp_stress_full.csv  all temperature rows + flags
#             output/checks/*.csv               data-audit reports
#             output/phylo/*                    trees
#### --------------------------------------------------  ####

	# Load the required libraries
		source("./R/func.R")
		check_and_install("pacman")
		pacman::p_load(tidyverse, flextable, latex2exp, metafor, orchaRd, readxl,
					   here, ggrepel, patchwork, rotl, ape, phytools, kutils, ggtree, janitor)

#### --------------------------------------------------  ####
# 1.1 Merge the two extraction files
#### --------------------------------------------------  ####

	# Both workbooks share the same 78-column extraction template, so they can be
	# row-bound directly. We assert that rather than assume it.
	#
	# Everything is read as text: readxl otherwise guesses column types per file
	# (e.g. exclude_missing_info comes back numeric from one workbook and character
	# from the other, which breaks the bind), and free-text entries such as "." in a
	# numeric column would be coerced to NA before we get a chance to see them.
	# Numeric columns are converted explicitly, once, below.
		read_sheet <- function(file) {
			x <- read_excel(here("data", file), sheet = "data", col_types = "text",
							.name_repair = "minimal")
			# the template has two columns literally named "units": the measurement
			# units (col 57) and the units of the dose column (col 63)
			nm <- names(x)
			nm[duplicated(nm)] <- paste0(nm[duplicated(nm)], "_dose")
			names(x) <- nm
			x
		}

		data_old <- read_sheet("mito_meta_data_merged_03072025.xlsx")
		data_new <- read_sheet("Kris_mito_meta_data_for_temp_OC.xlsx")

		stopifnot(identical(names(data_old), names(data_new)))

	# `source_file` keeps track of provenance; `study` must be unique across files
	# (no shared IDs: the new file uses WildNew1-8). Two papers in the new file were
	# already in the old extraction (Price et al. 2017 = s53; Schnell & Seebacher
	# 2008 = N17) and were deliberately NOT re-extracted -- see the `summary` sheet.
		stopifnot(length(intersect(unique(data_old$study), unique(data_new$study))) == 0)

		data <- bind_rows(
					data_old %>% mutate(source_file = "merged_03072025"),
					data_new %>% mutate(source_file = "kris_temp_OC")
				)

#### --------------------------------------------------  ####
# 1.2 Subset developmental temperature data
#### --------------------------------------------------  ####

	# `envirn_type` is the manipulation category. The extraction workbook also holds
	# rows for other manipulations (nutrition, cort, care deprivation, disturbance)
	# which are not analysed here.
		data <- data %>%
					mutate(envirn_type = str_squish(envirn_type)) %>%
					filter(envirn_type == "temp")

	# Trailing / leading whitespace in identifiers silently creates duplicate
	# factor levels (e.g. "s58 " vs "s58"), so squish every character key.
		data <- data %>%
					mutate(across(c(study, class, order, family, genus, species, common_name,
									stage, sex, envirn_type, admin, type, tissue, tissue_sum,
									mito_preparation, mito_ambiguity, measurement_category,
									respiration_category, antioxidant_category,
									oxidative_damage_category, units, measure_listed),
								  ~ str_squish(as.character(.x))))

	# Observation-level ID (used as the residual random effect) and a stable row key
	# that survives filtering, so audit reports can be traced back to a row.
		data <- data %>%
					mutate(row_key     = paste0(source_file, "_", row_number()),
						   observation = row_number())

#### --------------------------------------------------  ####
# 1.3 Harmonise moderator levels
#### --------------------------------------------------  ####

	# Free-text entry produced synonymous levels across extractors. Collapse them.
		data <- data %>%
					mutate(
						# tissue: "skeletal muscle"/"muscle", "whole animal"/"whole body" etc.
						tissue_sum = case_when(
							tissue_sum %in% c("brown adipose tissue (BAT)", "bat")   ~ "BAT",
							tissue_sum %in% c("whole animal", "whole body")          ~ "whole body",
							tissue_sum %in% c("whole blood", "blood")                ~ "blood",
							tissue_sum %in% c("skeletal muscle", "muscle")           ~ "muscle",
							tissue_sum == "adipose tissue"                           ~ "adipose",
							tissue_sum %in% c("serum", "plasma", "plasma/serum")     ~ "plasma/serum",
							tissue_sum == "erthrocyte"                               ~ "erythrocyte",
							TRUE                                                     ~ tissue_sum),
						# measurement category
						measurement_category = if_else(measurement_category == "gene expression",
													   "gene/protein expression", measurement_category),
						measurement_category = if_else(measurement_category == "oxidative stress",
													   "oxidative/nitrosative stress", measurement_category),
						# developmental stage
						stage = if_else(stage == "prenatal/postnatal", "both", stage))

	# `nutrition_sum` was reused to record a heat-wave design in Recio8; move that
	# information into a temperature-specific column so the nutrition columns can go.
		data <- data %>%
					mutate(temp_regime = if_else(str_detect(replace_na(nutrition_sum, ""), "heat wave"),
												 "heat wave", "constant"))

	# Temperature-specific derived moderators. t1 = treatment, t2 = control.
		data <- data %>%
					mutate(t1 = as.numeric(t1), t2 = as.numeric(t2),
						   temp_diff      = t1 - t2,
						   temp_direction = case_when(temp_diff  > 0 ~ "warmer",
													  temp_diff  < 0 ~ "cooler",
													  temp_diff == 0 ~ "equal",
													  TRUE           ~ NA_character_),
						   temp_magnitude = abs(temp_diff))

#### --------------------------------------------------  ####
# 1.4 Data checks  (writes audit reports; does not silently drop anything)
#### --------------------------------------------------  ####

		checks <- list()

	# (a) Missing or unusable summary statistics -------------------------------
		checks$missing_stats <- data %>%
			filter(if_any(c(mean_t1, sd_t1, n_t1, mean_t2, sd_t2, n_t2), ~ is.na(as.numeric(.x))) |
				   as.numeric(sd_t1) <= 0 | as.numeric(sd_t2) <= 0 |
				   as.numeric(n_t1)  <  2 | as.numeric(n_t2)  <  2) %>%
			select(row_key, study, tissue, descrp_measure, mean_t1, sd_t1, n_t1, mean_t2, sd_t2, n_t2)

	# (b) SD / SE internal consistency ----------------------------------------
	#     Convention is sd = se * sqrt(n). sd == se means an SEM was pasted into the
	#     SD column without conversion; any other mismatch implies a different n or
	#     a different error statistic than recorded.
		checks$sd_se <- data %>%
			mutate(across(c(sd_t1, se_t1, n_t1, sd_t2, se_t2, n_t2), as.numeric),
				   implied_n_t1 = (sd_t1 / se_t1)^2,
				   implied_n_t2 = (sd_t2 / se_t2)^2,
				   relerr_t1    = abs(sd_t1 - se_t1 * sqrt(n_t1)) / sd_t1,
				   relerr_t2    = abs(sd_t2 - se_t2 * sqrt(n_t2)) / sd_t2) %>%
			filter(replace_na(relerr_t1 > 0.02, FALSE) | replace_na(relerr_t2 > 0.02, FALSE)) %>%
			select(row_key, study, descrp_measure, sd_t1, se_t1, n_t1, implied_n_t1, relerr_t1,
											      sd_t2, se_t2, n_t2, implied_n_t2, relerr_t2)

	# (c) Implausible / invalid direction and category codes -------------------
		checks$dir_codes <- data %>%
			filter(!mito_efficiency_dir %in% c(0, 1, "0", "1")) %>%
			select(row_key, study, descrp_measure, measurement_category, mito_efficiency_dir)

	#     Same measure within a study coded inconsistently
		checks$dir_inconsistent <- data %>%
			mutate(m = str_to_lower(str_squish(descrp_measure))) %>%
			group_by(study, m) %>%
			filter(n_distinct(mito_efficiency_dir) > 1) %>%
			ungroup() %>%
			select(row_key, study, descrp_measure, mito_efficiency_dir, measurement_category)

		checks$cat_inconsistent <- data %>%
			mutate(m = str_to_lower(str_squish(descrp_measure))) %>%
			group_by(study, m) %>%
			filter(n_distinct(measurement_category) > 1) %>%
			ungroup() %>%
			select(row_key, study, descrp_measure, mito_efficiency_dir, measurement_category)

	# (d) Shared control groups -----------------------------------------------
	#     One control group compared against several treatment groups induces
	#     correlated effect sizes. Recorded, not removed: handled by the study /
	#     sample-dependency random effects.
		checks$shared_controls <- data %>%
			group_by(study, tissue, descrp_measure, mean_t2, sd_t2, n_t2) %>%
			filter(n() > 1) %>%
			summarise(n_reuse = n(), rows = paste(row_key, collapse = "; "), .groups = "drop")

	# (e) Suspicious value patterns -------------------------------------------
		checks$suspicious <- data %>%
			mutate(across(c(mean_t1, sd_t1, mean_t2, sd_t2), as.numeric)) %>%
			filter(mean_t1 == mean_t2 | (sd_t1 == sd_t2 & sd_t1 != 0)) %>%
			select(row_key, study, tissue, descrp_measure, mean_t1, sd_t1, mean_t2, sd_t2)

	# (f) Rows the extractors flagged ------------------------------------------
	#     exclude_missing_info: 0 = exclude (missing/ambiguous error or N),
	#     1 = keep, NA = not assessed (all new-file rows). NOTE the original script
	#     used `!exclude_missing_info == 0`, which drops every NA row through R's
	#     operator precedence; here the intent is written out explicitly.
		checks$flagged_exclude <- data %>%
			filter(!replace_na(as.character(exclude_missing_info), "1") %in% c("1", "NA")) %>%
			select(row_key, source_file, study, tissue, descrp_measure,
				   exclude_missing_info, mean_t1, sd_t1, n_t1, Notes)

		dir.create(here("output", "checks"), showWarnings = FALSE, recursive = TRUE)
		iwalk(checks, ~ write.csv(.x, here("output", "checks", paste0("check_", .y, ".csv")), row.names = FALSE))
		print(map_int(checks, nrow))

#### --------------------------------------------------  ####
# 1.5 Apply exclusions
#### --------------------------------------------------  ####

	# Flags rather than silent drops, so every decision is visible in the output and
	# can be toggled in the analysis.
		data <- data %>%
			mutate(
				# extractor-flagged missing/ambiguous information (0 = exclude).
				# "." appears once (WildNew3 liver COX) and is ambiguous -> query KHW;
				# treated as NOT excluded for now and reported in check_flagged_exclude.csv
				flag_missing_info = replace_na(as.character(exclude_missing_info), "1") == "0",
				# measures that are not mitochondrial (LDH, Na/K-ATPase, alkaline
				# phosphatase) or not a mitochondrial trait at all (plasma cort).
				# Kept as a flag rather than a hard drop because extraction round 2
				# deliberately recorded several such measures.
				flag_non_mito     = measurement_category %in% c("non-mitochondrial metabolic pathways",
																"glucocorticoids"),
				# unusable effect size
				flag_bad_stats    = row_key %in% checks$missing_stats$row_key)

		write.csv(data, here("output", "data", "temp_stress_full.csv"), row.names = FALSE)

	# Analysis dataset.
		data <- data %>% filter(!flag_missing_info, !flag_non_mito, !flag_bad_stats)

#### --------------------------------------------------  ####
# 2. Phylogeny
#### --------------------------------------------------  ####

	# Species binomials, with spelling corrections needed to match the Open Tree of
	# Life. Several are misspellings in the extraction sheet, and Dicentrarchus
	# labrax was entered under two different misspellings.
		data <- data %>%
			mutate(species_phylo = paste(genus, species, sep = "_"),
				   species_phylo = str_replace_all(species_phylo, " ", "_"),
				   species_phylo = case_when(
						species_phylo == "Meleagris_gallopavo_domesitcus"              ~ "Meleagris_gallopavo",
						species_phylo == "Cortunix_japonica"                           ~ "Coturnix_japonica",
						species_phylo == "Symphysodon_aequifasciatus"                   ~ "Symphysodon_aequifasciata",
						species_phylo == "Oncorhynchus_tschawyscha"                     ~ "Oncorhynchus_tshawytscha",
						species_phylo %in% c("Dicentrarachus_labrax",
											 "Dichentrarchus_labrax")                   ~ "Dicentrarchus_labrax",
						TRUE                                                            ~ species_phylo),
				   species_phylo2 = species_phylo)

	# Match names against the Open Tree of Life and build the tree
		tol_subtree <- rotl::tnrs_match_names(unique(data$species_phylo))
		print(tol_subtree[, c("search_string", "unique_name", "approximate_match")])

		tree <- rotl::tol_induced_subtree(tol_subtree$ott_id, label_format = "name")

	# Strip the "(genus in ...)" annotations OTL adds to ambiguous labels
		tree$tip.label <- gsub("_\\([^)]*\\)", "", tree$tip.label)

		dir.create(here("output", "phylo"), showWarnings = FALSE, recursive = TRUE)
		tree_checks(data.frame(data), tree, dataCol = "species_phylo")

		write.tree(tree, here("output", "phylo", "phylo.tre"))
		write.table(tree$tip.label, here("output", "phylo", "phylo_species.txt"),
					row.names = FALSE, col.names = FALSE)

	# Prune to the species actually present, then re-check
		tree <- tree_checks(data.frame(data), tree, dataCol = "species_phylo", type = "prune")
		tree_checks(data.frame(data), tree, dataCol = "species_phylo")

		write.tree(tree, here("output", "phylo", "phylo_pruned.tre"))
		write.table(gsub("_", " ", tree$tip.label), here("output", "phylo", "phylo_pruned_species.txt"),
					quote = FALSE, row.names = FALSE, col.names = FALSE)

	# The rotl tree above is topology only (no branch lengths). The time-calibrated
	# tree used in the models, output/phylo/phylo_pruned_species.nwk, was built by
	# uploading phylo_pruned_species.txt to TimeTree 5 (timetree.org, "Load a List of
	# Species") on 16 Sep 2026 and exporting to Newick. This is not re-run here --
	# redo it if the species list changes. TimeTree resolved all 29 names; it had no
	# data for Chiloscyllium plagiosum and used its congener C. punctatum instead
	# (the tip keeps the plagiosum label). The TimeTree and rotl topologies differ in
	# one split: TimeTree pairs Gallus with Meleagris, rotl pairs Gallus with Coturnix.
		stopifnot(setequal(read.tree(here("output", "phylo", "phylo_pruned_species.nwk"))$tip.label,
						   unique(data$species_phylo)))

	# Plot
		plot_tree <- ggtree(tree) +
			geom_tiplab(aes(label = gsub("_", " ", label)), size = 7.5, offset = 0.1, hjust = 0, align = FALSE) +
			scale_x_continuous(expand = expansion(mult = c(0, 0.8)))
		ggsave(here("output", "phylo", "phylo.png"), plot_tree, width = 22.888889, height = 8.604938)

#### --------------------------------------------------  ####
# 3. Effect size calculations
#### --------------------------------------------------  ####

	# SMDH (standardised mean difference allowing heteroscedastic population
	# variances; Bonett 2008, 2009). The data contain ratios (RCR, relative
	# expression), percentages, zeros and skewed measurement variables, all of which
	# complicate lnRR.
		data <- data %>% mutate(across(c(mean_t1, sd_t1, n_t1, mean_t2, sd_t2, n_t2), as.numeric))

		data <- metafor::escalc(measure = "SMDH",
								m1i = mean_t1, m2i = mean_t2,
								sd1i = sd_t1,  sd2i = sd_t2,
								n1i  = n_t1,   n2i  = n_t2,
								data = data, var.names = c("SMDH", "v_SMDH"))

	# Direction correction. Positive SMDH must always mean HIGHER mitochondrial
	# respiratory function. mito_efficiency_dir == 1: an increase in the measure
	# means higher function (no change). == 0: an increase means lower function
	# (flip). Codes are validated above, so anything else is an error, not a silent
	# pass-through.
		stopifnot(all(data$mito_efficiency_dir %in% c(0, 1, "0", "1")))

		data <- data %>%
			mutate(SMDH = if_else(as.character(mito_efficiency_dir) == "0", -1 * SMDH, SMDH))

	# Precision terms used by the publication-bias models
		data <- data %>%
			mutate(inv_n      = (n_t1 + n_t2) / (n_t1 * n_t2),
				   sqrt_inv_n = sqrt(inv_n),
				   ef_n       = (4 * n_t1 * n_t2) / (n_t1 + n_t2),   # effective sample size
				   depend     = interaction(sample_depend, study))

#### --------------------------------------------------  ####
# 4. Exploratory plots and the analysis dataset
#### --------------------------------------------------  ####

		dir.create(here("output", "fig_explore"), showWarnings = FALSE, recursive = TRUE)

	# Mean-variance relationship, by group
		control <- ggplot(data, aes(x = log(mean_t2), y = log(sd_t2))) +
			geom_point() + geom_smooth(method = "lm", se = TRUE) +
			geom_label_repel(aes(label = study), max.overlaps = 60, box.padding = 0.5,
							 point.padding = 0.5, segment.color = "grey50") +
			labs(title = "Mean-variance relationship, control", x = "log Mean", y = "log SD")

		trt <- ggplot(data, aes(x = log(mean_t1), y = log(sd_t1))) +
			geom_point() + geom_smooth(method = "lm", se = TRUE) +
			geom_label_repel(aes(label = study), max.overlaps = 60, box.padding = 0.5,
							 point.padding = 0.5, segment.color = "grey50") +
			labs(title = "Mean-variance relationship, treatment", x = "log Mean", y = "log SD")

		ggsave(here("output", "fig_explore", "fig_ex1.png"),
			   (control + my_theme() | trt + my_theme()) + plot_annotation(tag_levels = "A", tag_suffix = ")"),
			   width = 22.888889, height = 8.604938)

	# Effect size distribution and funnel
		ggsave(here("output", "fig_explore", "fig_ex2_SMDH.png"),
			   ggplot(data, aes(x = SMDH)) + geom_histogram(bins = 30) + my_theme() +
					labs(title = "Distribution of SMDH", x = "SMDH", y = "Frequency"),
			   width = 11, height = 8.6)

	# Extreme effects are kept but written out for sensitivity analysis
		write.csv(data %>% filter(abs(SMDH) > 5) %>%
					  dplyr::select(row_key, study, descrp_measure, units, mean_t1, sd_t1, n_t1,
									mean_t2, sd_t2, n_t2, SMDH, v_SMDH),
				  here("output", "checks", "check_extreme_effects.csv"), row.names = FALSE)

	# Drop extraction-workflow columns that carry no information for temperature
		temp_stress <- data %>%
			dplyr::select(!any_of(c("nutrition_sum", "nutrition_type", "fasting_period",
									"fasting_period_units", "CORT_values_available",
									"CI95_t1", "CI95_t2", "num_tissues_types")))

		dir.create(here("output", "data"), showWarnings = FALSE, recursive = TRUE)
		write.csv(temp_stress, here("output", "data", "temp_stress.csv"), row.names = FALSE)

	# Summary of the analysis dataset
		summary_table <- temp_stress %>%
			summarise(k = n(), spp = n_distinct(species_phylo), studies = n_distinct(study)) %>%
			mutate(envirn_type = "Temperature") %>%
			dplyr::select(envirn_type, everything())

		by_category <- temp_stress %>%
			group_by(measurement_category) %>%
			summarise(k = n(), studies = n_distinct(study), spp = n_distinct(species_phylo)) %>%
			arrange(desc(k))

		by_class <- temp_stress %>%
			group_by(class) %>%
			summarise(k = n(), studies = n_distinct(study), spp = n_distinct(species_phylo)) %>%
			arrange(desc(k))

		dir.create(here("output", "tables"), showWarnings = FALSE, recursive = TRUE)
		write.csv(summary_table, here("output", "tables", "data_sum_table.csv"), row.names = FALSE)
		write.csv(by_category,   here("output", "tables", "data_sum_by_category.csv"), row.names = FALSE)
		write.csv(by_class,      here("output", "tables", "data_sum_by_class.csv"), row.names = FALSE)

		print(summary_table); print(by_category); print(by_class)
