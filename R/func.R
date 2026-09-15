#### -------------------------------------------- ####
#             Functions for analysis
#### -------------------------------------------- ####

#' @title check_and_install
#' @description Function to check and install a package
#' @param pkg Package name
check_and_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  } else {
    message(paste(pkg, "is already installed."))
  }
}

# Plot themes
my_theme  <- function() {
		list(theme_classic(),
			 theme(axis.text.x = element_text(size = 12),
				   axis.text.y = element_text(size = 12),
				  axis.title.x = element_text(size = 16),
				  axis.title.y = element_text(size = 16),
				      plot.tag = element_text(size = 20)))
	}

#' @title tree_checks
#' @description Checks and prunes trees
#' @param data Dataframe with species names
#' @param tree Phylogenetic tree
#' @param dataCol Column name of species in data
#' @param type Character. "checks" will return a list of species in the tree but not in the data and vice versa. "prune" will prune the tree of species not in the data.
 
  tree_checks <- function(data, tree, dataCol, type = c("checks", "prune")){
    type = match.arg(type)
    # How many unique species exist in data and tree
    Numbers <- matrix(nrow = 2, ncol = 1)
    Numbers[1,1] <- length(unique(data[,dataCol])) 
    Numbers[2,1] <- length(tree$tip.label) 
    rownames(Numbers)<- c("Species in data:", "Species in tree:")
    # Missing species or species not spelt correct      
    species_list1= setdiff(sort(tree$tip.label), sort(unique(data[,dataCol])))
    species_list2= setdiff(sort(unique(data[,dataCol])), sort(tree$tip.label) )
    if(type == "checks"){
      return(list(SpeciesNumbers = data.frame(Numbers), 
                  Species_InTree_But_NotData=species_list1, 
                  Species_InData_But_NotTree=species_list2))
    }
    if(type == "prune"){
      if(length(species_list2) >=1) stop("Sorry, you can only prune a tree when you have no taxa existing in the data that are not in the tree")
      return(ape::drop.tip(tree, species_list1))
    }
  }

#' @title round_df
#' @description Rounds all numeric columns in a dataframe
#' @param df Dataframe
#' @param digits Number of digits to round to
#' @return Dataframe with numeric columns rounded
round_df <- function(df, digits = 2) {
  df[] <- lapply(df, function(x) if (is.numeric(x)) round(x, digits) else x)
  return(df)
}


#' @title make_het_tables
#' @description Calculates heterogeneity statistics from a list of models
#' @param mod_list List of models
#' @param size Dataframe with size information (e.g., number of studies, species)
#' @param type Type of heterogeneity statistic to calculate. Options are "i2", "cv2", or "m2".
#' @details This function extracts heterogeneity statistics from a list of models and combines them into a single dataframe. It supports different types of heterogeneity statistics
#' (i2, cv2, m2) and can handle multiple models.
#' @return Dataframe with numeric columns rounded
make_het_tables  <- function(mod_list, size, type = c("i2", "cv2", "m2")){
  type = match.arg(type)

  het  <- lapply(mod_list[c(1,3,4)], function (x) switch(type, "i2" = orchaRd::i2_ml(x),
                                                              "cv2" = orchaRd::cvh2_ml(x),
                                                               "m2" = orchaRd::m2_ml(x)))
  het  <- data.frame(do.call("rbind", het))
  
  het.2 <- lapply(mod_list[2], function (x) switch(type, "i2" = orchaRd::i2_ml(x),
                                                        "cv2" = orchaRd::cvh2_ml(x),
                                                         "m2" = orchaRd::m2_ml(x)))
  het.2 <- data.frame(do.call("rbind", het.2))
  
  het[4,c(1,2,5)] <- het.2
  return(cbind(size, het))
}
