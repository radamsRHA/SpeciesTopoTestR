#' Optimize.TopologySearch: function to find the ML topology using STELLS
#'
#' This function returns an optimized species topology and its maximized likelihood
#' @param handle.GeneTrees Phylo object containing a list of gene trees
#' @param numeric.Stells_algorthm Numeric defining with algorithm of STELLS to use: 0 or 1
#' @param string.PathDir String defining the path to a parent directory used for optimizing species tree
#' @keywords Species tree, multispecies coalescent, phylogenetics, phylogenomics
#' @return handle.Optimized_SpeciesTree Species tree with branch lengths (in coalescent units) that has been optimized with STELLS
#' @return numeric.MaximizedLnL Numeric containing the maximum likelihood given the optimized branch lengths of the species tree
#' @export
#' @examples
#'
#'

###########################
# Optimize.TopologySearch #
###########################
Optimize.TopologySearch <- function(handle.GeneTrees, string.PathDir){

  #################################################################################
  # Define directory used for computing gene trees likelihoods given species tree #
  #################################################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_OptimizeTopologySearch = paste(string.PathDir, '/Optimize_TopologySearch', Sys.Date(), sep = "")
  unlink(string.Path_Directory_OptimizeTopologySearch, recursive = T)
  dir.create(string.Path_Directory_OptimizeTopologySearch, showWarnings = T, recursive = T)
  setwd(string.Path_Directory_OptimizeTopologySearch)

  ############################
  # Write trees to directory #
  ############################
  numeric.NumberOfGeneTrees <- length(handle.GeneTrees)
  write.tree(phy = handle.GeneTrees, file = paste0(string.Path_Directory_OptimizeTopologySearch, "/GeneTrees.tree"))

  ##############
  # Run STELLS #
  ##############
  if (Sys.info()["sysname"] == "Darwin"){
    string.Command_STELLS <- "stells-v2-1-0-1-mac -g GeneTrees.tree"
  }
  if (Sys.info()["sysname"] == "Linux"){
    string.Command_STELLS <- "stells-v2-1-0-1-linux64 -g GeneTrees.tree"
  }

  #string.Command_STELLS <- "stells-v2-1-0-1-mac -g GeneTrees.tree"
  handle.STELLS_Results <- system(command = string.Command_STELLS, intern = T)
  print(handle.STELLS_Results)

  #####################################################################
  # Read results from STELLS and find optimum LnL and newick topology #
  #####################################################################
  numeric.Grep_Line_LnL <- grep(pattern = "Highest log likelihood of tree (found at tree space exploring) = ", x = handle.STELLS_Results, fixed = T)
  string.Line_LnL <- handle.STELLS_Results[[numeric.Grep_Line_LnL]]
  string.Split_Line_LnL <- strsplit(x = string.Line_LnL, split = "Highest log likelihood of tree (found at tree space exploring) = ", fixed = T)[[1]][2]
  numeric.Optimized_LnL <- as.numeric(string.Split_Line_LnL)

  numeric.Grep_Line_Newick <- grep(pattern = "The newick format of the inferred MLE species tree: ", x = handle.STELLS_Results)
  string.Line_Newick <- handle.STELLS_Results[[numeric.Grep_Line_Newick]]
  string.Split_Line_Newick <- strsplit(x = string.Line_Newick, split = "The newick format of the inferred MLE species tree: ")[[1]][2]
  string.Split_Line_Newick <- paste0(string.Split_Line_Newick, ';')
  handle.Optimized_SpeciesTree <- read.tree(text = string.Split_Line_Newick)

  #####################################################################
  # Return to directory and return the LnL and optimized species tree #
  #####################################################################
  setwd(string.CurrentDir)
  return(list(handle.Optimized_SpeciesTree = handle.Optimized_SpeciesTree,
              numeric.MaximizedLnL = numeric.Optimized_LnL))




}
