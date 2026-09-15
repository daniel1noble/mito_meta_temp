#### --------------------------------------------------  ####
# 1. Data Checking, processing and exploratory analysis
#### --------------------------------------------------  ####

	# Load the required libraries
		source("./R/func.R")
		check_and_install("pacman")
		install.packages("readxl"); library(readxl)
		devtools::install_github("daniel1noble/orchaRd")
		pacman::p_load(tidyverse, flextable, latex2exp, metafor, orchaRd, readxl, here, ggrepel, patchwork, rotl, ape, phytools, kutils, ggtree, janitor)

	# Load the data
		data <- read_excel(here("data", "mito_meta_data_merged_03072025.xlsx"))

	# Check the data
		str(data)

	# Add observation-level column
		data <- data %>%
			mutate(observation = row_number())

	# Check mean-variance relationship
		control <- ggplot(data, aes(x = log(mean_t1), y = log(sd_t1))) +
			geom_point() +
			geom_smooth(method = "lm", se = TRUE) +
			geom_label_repel(aes(label = study), max.overlaps = 60, box.padding = 0.5, point.padding = 0.5, segment.color = "grey50") +
			labs(title = "Mean-variance relationship of control", x = "log Mean", y = " log Standard Deviation")

		trt <- ggplot(data, aes(x = log(mean_t2), y = log(sd_t2))) +
			geom_point() +
			geom_smooth(method = "lm", se = TRUE) +
			geom_label_repel(aes(label = study), max.overlaps = 60, box.padding = 0.5, point.padding = 0.5, segment.color = "grey50") +
			labs(title = "Mean-variance relationship of treatment", x = "log Mean", y = " log Standard Deviation")

	plot <- (control + my_theme() | trt + my_theme()) + plot_annotation(tag_levels = "A", tag_suffix = ")")
	ggsave(here("output", "fig_explore", "fig_ex1.png"), plot, width = 22.888889, height = 8.604938)

	# Check levels of categorical variables that are of major interest to make sure now spelling errors etc
	unique(data$envirn_type)
	unique(data$measurement_category) # Some issues here with gene expression. Probably just call "gene/protein expression"
	unique(data$tissue_sum) # Some issues here "BAT" == "brown adipose tissue (BAT)" == "bat"; "whole body" == "whole animal"; "whole blood" == "blood"; "skeletal muscle" == "muscle"; 
	unique(data$stage) # re-cat pre/post to both

	# Check effects/studies in some measure categories
	tabyl(data, measurement_category) %>% arrange(desc(n))
	data  %>% group_by(measurement_category) %>% summarise(studies = n_distinct(study), effects = n()) %>% arrange(desc(effects))

	# Fix the issues
	data <- data %>%
		mutate(measurement_category = ifelse(measurement_category == "gene expression", "gene/protein expression", measurement_category),
				tissue_sum = ifelse(tissue_sum == "brown adipose tissue (BAT)", "BAT", tissue_sum),
				tissue_sum = ifelse(tissue_sum == "bat", "BAT", tissue_sum),
				tissue_sum = ifelse(tissue_sum == "whole animal", "whole body", tissue_sum),
				tissue_sum = ifelse(tissue_sum == "whole blood", "blood", tissue_sum),
				tissue_sum = ifelse(tissue_sum == "skeletal muscle", "muscle", tissue_sum),
				tissue_sum = ifelse(tissue_sum == "adipose tissue", "adipose", tissue_sum),
				tissue_sum = ifelse(tissue_sum == "serum", "plasma/serum", tissue_sum),
				tissue_sum = ifelse(tissue_sum == "erthrocyte", "erythrocyte", tissue_sum),
				     stage = ifelse(stage == "prenatal/postnatal", "both", stage))  %>% 
		filter(!measurement_category %in% c("non-mitochondrial metabolic pathways", "glucocorticoids"))  %>% filter(!exclude_missing_info == 0) %>% data.frame() # # Drop non-metabolic pathways and cort and Drop studies and or rows within studies that have missing (missing error, N) and/or ambiguous (uncertain errors) information

#### --------------------------------------------------  ####
# 2. Phylogeny
#### --------------------------------------------------  ####
	 
	# From these data we will create a phylogeny which we can trim based on data subsets. Create species first
	 data <- data %>%
		 mutate(species_phylo = paste(genus, species, sep = "_"),
		 		species_phylo = ifelse(species_phylo == "Meleagris_gallopavo domesitcus", "Meleagris_gallopavo", species_phylo))

	# Few more fixes in data to match species in data and tree
	 data <- data %>%
		 mutate(species_phylo = ifelse(species_phylo == "Trachemys_scripta elegans", "Trachemys_scripta_elegans", 
		 								if_else(species_phylo == "Oncorhynchus_tschawyscha", "Oncorhynchus_tshawytscha", 
										if_else(species_phylo == "Tamiasciurus_hudonicus", "Tamiasciurus_hudsonicus",
										if_else(species_phylo == "Cortunix_japonica", "Coturnix_japonica", 
										if_else(species_phylo == "Symphysodon_aequifasciatus", "Symphysodon_aequifasciata",
										if_else(species_phylo %in% c("Dicentrarachus_labrax", "Dichentrarchus_labrax"), "Dicentrarchus_labrax", species_phylo)))))),
		 species_phylo2 = species_phylo) 

	 # Create a phylogeny. Looks like all species are matched and in the tree
	 tol_subtree  <- rotl::tnrs_match_names(unique(data$species_phylo))

	 # Create a phylogeny
	 tree <- rotl::tol_induced_subtree(tol_subtree$ott_id, label_format = "name")

	 # Check the tree. Some warnings
	 length(unique(data$species_phylo)); length(tree$tip.label) # lost two species, which ones
	
	 # Fix labels so they match the data
	 tree$tip.label <- gsub("_\\([^)]*\\)", "", tree$tip.label)

	 # Check the tree again
		tree_checks(data, tree, dataCol = "species_phylo")
	 
	 # Write final tree
	 	write.tree(tree, here("output", "phylo", "phylo.tre"))
		write.table(tree$tip.label, here("output", "phylo", "phylo_species.txt"), row.names = FALSE, col.names = FALSE)

	 # Plot the tree
	 plot_tree <- ggtree(tree) + geom_tiplab(aes(label = gsub("_", " ", label)), size = 7.5, offset = 0.1, hjust = 0, align = FALSE) + scale_x_continuous(expand = expansion(mult = c(0, 0.8))) 
	 ggsave(here("output", "phylo", "phylo.png"), plot_tree, width = 22.888889, height = 8.604938)

	 # Check the tree
	 tree_checks(data, tree, dataCol = "species_phylo")

	 # Prune the tree
	 tree <- tree_checks(data, tree, dataCol = "species_phylo", type = "prune")

	 # Check the tree
	 tree_checks(data, tree, dataCol = "species_phylo")

	 # Write final tree. Can use the species names to try and get timetree as well. 
	 write.tree(tree, here("output", "phylo", "phylo_pruned.tre"))
	 write.table(gsub("_", " ", tree$tip.label), here("output", "phylo", "phylo_pruned_species.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

#### --------------------------------------------------  ####
# 3. Effect size calculations
#### --------------------------------------------------  ####

	## Let's calculate the effect size. We will use SMD with these data for a number of reasons. First, we have ratio data (RCR, gene expression which is relative etc) which can complicate lnRR interpretation. Second, we have lots of percentages. This is fine for lnRR and we can do a logit conversion but it can lead to some issues whereas SMD is a little more robust. We have also have some zero / negative values. Though these are rare, but are not defined with lnRR. 
		
		data <-  metafor::escalc(measure = "SMDH", m1i = mean_t1, m2i = mean_t2, 
											sd1i = sd_t1, sd2i = sd_t2, 
											n1i = n_t1 , n2i = n_t2, 
						data = data, var.names = c("SMDH", "v_SMDH"))

	# Have a look at effect size distributions. Some large effects we need to check, but pretty typical 
	plot_SMDH <- ggplot(data, aes(x = SMDH)) +
						geom_histogram(bins = 30) +
						labs(title = "Distribution of SMDH", x = "SMDH", y = "Frequency")

	# Lets have a look at funnel plots. No obvious issues.
	plot_funnel_SMDH <- metafor::funnel(data$SMDH, data$v_SMDH, yaxis="seinv")

	# Lest have a look at the extreme effect sizes to see if there are problems
		#write.csv(data  %>% filter(abs(SMDH) > 5)  %>% dplyr::select(study, descrp_measure, units...42, mean_t1, sd_t1, n_t1, mean_t2, sd_t2, n_t2, SMDH, v_SMDH), here("output", "checks", "extreme_effects.csv"), row.names = FALSE)

# We need to correct the direction for mito measurements as they have a different meaning depending on whether one groups mean is higher than the other. 
    head(data  %>%  filter(mito_efficiency_dir == 1)  %>% select(study, mito_efficiency_dir, SMDH, v_SMDH))
    head(data  %>%  filter(mito_efficiency_dir == 0)  %>% select(study, mito_efficiency_dir, SMDH, v_SMDH))
	head(data  %>%  filter(mito_efficiency_dir == "NA")  %>% select(study, mito_efficiency_dir, SMDH, v_SMDH))
	
	data  <- data %>%
	            mutate(SMDH = ifelse(mito_efficiency_dir %in% c(1, "NA"), SMDH, -1*SMDH))

#### --------------------------------------------------  ####
# 4. Data subsets
#### --------------------------------------------------  ####

	# Subset the nutrition stress data and remove columns that we don't seem to need
	nutri_stress <- data %>%
				filter(envirn_type == "nutrition")  %>% select(!c(Notes, CORT_values_available, CI95_t2, se_t2, CI95_t1, se_t1, num_tissues_types, admin, type, control_multiple_comparisons, measurement_methods))
	write.csv(nutri_stress, here("output","data",  "nutri_stress.csv"), row.names = FALSE)
	sum_nutri_stress <- nutri_stress  %>% summarise(k = n(), spp = n_distinct(species_phylo), studies = n_distinct(study))

	# Subset the CORT stress data and remove columns that we don't seem to need
	cort_stress <- data %>%
				filter(envirn_type == "cort")  %>% select(!c(Notes, CI95_t2, se_t2, CI95_t1, se_t1, num_tissues_types, admin, type, control_multiple_comparisons, measurement_methods))
	write.csv(cort_stress, here("output", "data", "cort_stress.csv"), row.names = FALSE)
	sum_cort_stress <- cort_stress  %>% summarise(k = n(), spp = n_distinct(species_phylo), studies = n_distinct(study))

	# Subset the deprivation stress data and remove columns that we don't seem to need
	deprive_stress <- data %>%
				filter(envirn_type == "care deprivation")  %>% select(!c(Notes, CI95_t2, se_t2, CI95_t1, se_t1, num_tissues_types, admin, type, control_multiple_comparisons, measurement_methods))
	write.csv(deprive_stress, here("output", "data", "deprive_stress.csv"), row.names = FALSE)
	sum_deprive_stress <- deprive_stress  %>% summarise(k = n(), spp = n_distinct(species_phylo), studies = n_distinct(study))

	# Subset the disturbance stress data and remove columns that we don't seem to need
	disturb_stress <- data %>%
				filter(envirn_type == "disturbance")  %>% select(!c(Notes, CI95_t2, se_t2, CI95_t1, se_t1, num_tissues_types, admin, type, control_multiple_comparisons, measurement_methods))
	write.csv(disturb_stress, here("output", "data", "disturb_stress.csv"), row.names = FALSE)
	sum_disturb_stress <- disturb_stress  %>% summarise(k = n(), spp = n_distinct(species_phylo), studies = n_distinct(study))

	# Summary table of all teh data subsets
	summary_table <- rbind(sum_nutri_stress, sum_cort_stress, sum_deprive_stress, sum_disturb_stress)
	summary_table <- summary_table %>% mutate(envirn_type = c("Nutrition", "CORT", "Care Deprivation", "Disturbance"))  %>% select(envirn_type, everything())
	write.csv(summary_table, here("output", "tables", "data_sum_table.csv"), row.names = FALSE)

# Check how many studies for temperature which we do not focus on because of challenges with interpretation
	temps <- data  %>% filter(envirn_type == "temp")
	unique(temps$study) # 33 studies
	dim(temps) # 459 effects

any(data$study == "s73") # Check if this study is in the data. It is not, so we can remove it from the list of studies to check





