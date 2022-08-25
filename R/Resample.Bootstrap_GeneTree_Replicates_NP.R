#' Resample.Bootstrap_GeneTree_Replicates_NP: function to conduct non-parametric bootstrapping (i.e., random sampling with replacement) for a set of input gene trees
#'
#' This function returns a list containing multiple replicate sets of gene trees (each replicate is a bootstrapped dataset)
#' @param handle.GeneTrees MultiPhylo object of the input gene trees
#' @param numeric.NumberOfReps Numeric number of bootstrap replicates
#' @keywords Species tree, multispecies coalescent, phylogenetics, phylogenomics
#' @return list.BootStrap_GeneTree_ReplicateSets List containing sets of gene trees bootstrap replicates
#' @export
#' @examples
#'
#' #'
#' ################
#' # Load depends #
#' ################
#' library(SpeciesTopoTestR)
#' library(ape)
#' 
#' #################################
#' # Generate example species tree #
#' #################################
#' handle.SpeciesTree <- read.tree(text = "(A:1.0,(B:0.5,(C:1.0,D:1.0):1.5):0.5);")
#' 
#' #####################################################
#' # Simlate a set of gene trees for this species tree #
#' #####################################################
#' handle.Simulated_GeneTrees <- Simulate.GeneTrees_From_SpeciesTree3(handle.SpeciesTree = handle.SpeciesTree, 
#'                                                                   string.PathDir = '~/Desktop/', 
#'                                                                   numeric.NumberOfGeneTrees = 10)
#' 
#' #################################################
#' # Resample from example simulated gene tree set #
#' #################################################
#' list.Bootstrapped_GeneTreesSets <- Resample.Bootstrap_GeneTree_Replicates_NP(handle.GeneTrees = handle.Simulated_GeneTrees, numeric.NumberOfReps = 10)
#' 
#' 
#' 


############################################
# Resample.Bootstrap_GeneTree_Replicates_NP #
############################################
Resample.Bootstrap_GeneTree_Replicates_NP <- function(handle.GeneTrees, numeric.NumberOfReps){
  
  ###################
  # Summarize input #
  ###################
  string.CurrentDir <- getwd()
  numeric.NumberOfGenes <- length(handle.GeneTrees)
  list.BootStrap_GeneTree_ReplicateSets <- list()
  
  ####################################
  # Loop through replicate samplings #
  ####################################
  for (i in 1:numeric.NumberOfReps){
    
    ##############################
    # Randomly sample gene trees #
    ##############################
    list.BootStrap_GeneTree_ReplicateSets[[i]] <- sample(x = handle.GeneTrees, size = numeric.NumberOfGenes, replace = T)
    
  }
  
  ##########################################################################
  # Return to directory and return list of bootstrapped sets of gene trees #
  ##########################################################################
  setwd(string.CurrentDir)
  return(list.BootStrap_GeneTree_ReplicateSets)
}