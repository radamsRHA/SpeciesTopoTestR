#' Compute.GeneTree_Likelihoods: function to compute gene tree likelihoods using STELL
#'
#' This function returns a vector of gene tree likelihoods given a set of gene trees and a particular species tree
#' @param handle.GeneTrees Phylo object of gene trees
#' @param handle.SpeciesTree Phylo object of the species tree
#' @param string.PathDir String of the path to a parent directory used for simulating gene trees
#' @keywords Species tree, multispecies coalescent, phylogenetics, phylogenomics
#' @return vector.GeneTreeProbs Vector of gene tree likelihoods for each gene tree provided
#' @export
#' @examples
#'
#' ################
#' # Load depends #
#' ################
#' library(SpeciesTopoTestR)
#' library(ape)
#'
#' #################################
#' # Generate example species trees #
#' #################################
#' handle.SpeciesTree1 <- read.tree(text = "(A:1.0,(B:0.5,(C:1.0,D:1.0):1.5):0.5);")
#'
#' #################################################
#' # Simlate a set of gene trees for species tree1 #
#' #################################################
#' handle.Simulated_GeneTrees <- Simulate.GeneTrees_From_SpeciesTree3(handle.SpeciesTree = handle.SpeciesTree1,
#'                                                                   string.PathDir = '~/Desktop/',
#'                                                                   numeric.NumberOfGeneTrees = 100)
#'
#'
#' vector.GeneTree_Likelihoods <- Compute.GeneTree_Likelihoods(handle.SpeciesTree = handle.SpeciesTree1,
#'                                                             handle.GeneTrees = handle.Simulated_GeneTrees,
#'                                                             string.PathDir = '~/Desktop/')
#'

################################
# Compute.GeneTree_Likelihoods #
################################
Compute.GeneTree_Likelihoods <- function(handle.SpeciesTree, handle.GeneTrees, string.PathDir){

  #################################################################################
  # Define directory used for computing gene trees likelihoods given species tree #
  #################################################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_GeneTreeLikelihoods = paste(string.PathDir, '/GeneTreeLikelihoods_', Sys.Date(), sep = "")
  unlink(string.Path_Directory_GeneTreeLikelihoods, recursive = T)
  dir.create(string.Path_Directory_GeneTreeLikelihoods, showWarnings = T, recursive = T)
  setwd(string.Path_Directory_GeneTreeLikelihoods)

  ###################
  # Summarize input #
  ###################
  numeric.NumberOfGeneTrees <- length(handle.GeneTrees)
  string.SpeciesTree <- write.tree(phy = handle.SpeciesTree, file = "")
  vector.GeneTreeProbs <- rep(NA, numeric.NumberOfGeneTrees)

  ############################
  # Write trees to directory #
  ############################
  write.tree(phy = handle.SpeciesTree, file = paste0(string.Path_Directory_GeneTreeLikelihoods, "/SpeciesTree.tree"))
  write.tree(phy = handle.GeneTrees, file = paste0(string.Path_Directory_GeneTreeLikelihoods, "/GeneTrees.tree"))

  ##################################
  # Command to compute likelihoods #
  ##################################
  if (Sys.info()["sysname"] == "Darwin"){
    string.Command_STELLS <- "stells-v2-1-0-1-mac -g GeneTrees.tree -s SpeciesTree.tree -B"
  }
  if (Sys.info()["sysname"] == "Linux"){
    string.Command_STELLS <- "stells-v2-1-0-1-linux64 -g GeneTrees.tree -s SpeciesTree.tree -B"
  }

  #string.Command_STELLS <- "stells-v2-1-0-1-mac -g GeneTrees.tree -s SpeciesTree.tree -B"
  handle.STELLS_Results <- system(command = string.Command_STELLS, intern = T)
  numeric.GeneTree_Count <- 0
  vector.GeneTreeLikelihoods <- rep(NA, numeric.NumberOfGeneTrees)

  ##############################
  # Find gene tree likelihoods #
  ##############################
  vector.Grep_Lines_GeneTreeLiks <- grep(pattern = "Log-prob of Tree ", x = handle.STELLS_Results)

  for (i in 1:length(vector.Grep_Lines_GeneTreeLiks)){

    string.GeneTree_LnL_i <- handle.STELLS_Results[[vector.Grep_Lines_GeneTreeLiks[i]]]
    string.GeneTree_LnL_Split_i <- strsplit(x = string.GeneTree_LnL_i, paste0("Log-prob of Tree ", (i-1), ' '))[[1]][2]
    numeric.GeneTree_LnL <- as.numeric(string.GeneTree_LnL_Split_i)
    vector.GeneTreeLikelihoods[i] <- numeric.GeneTree_LnL
  }

  ##########################################
  # Return to directory and return results #
  ##########################################
  setwd(string.CurrentDir)
  return(vector.GeneTreeLikelihoods)

}
