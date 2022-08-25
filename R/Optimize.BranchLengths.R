#' Optimize.BranchLengths: function to optimize branch lengths (in coalescent units) of an input species tree topology given a set of input gene tree topologies
#'
#' This function returns a species tree with branch lengths in coalescent units, which have been optimized using the MSC and the STELLS2 algorithm
#' @param handle.GeneTrees Phylo object containing a list of gene trees
#' @param handle.SpeciesTree Phylo object defining the input species tree
#' @param numeric.Stells_algorthm Numeric defining with algorithm of STELLS to use: 0 or 1
#' @param string.PathDir String defining the path to a parent directory used for optimizing species tree
#' @keywords Species tree, multispecies coalescent, phylogenetics, phylogenomics
#' @return handle.Optimized_SpeciesTree Species tree with branch lengths (in coalescent units) that has been optimized with STELLS
#' @return numeric.MaximizedLnL Numeric containing the maximum likelihood given the optimized branch lengths of the species tree
#' @export
#' @examples
#'
#'
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
#' #############################################################################
#' # Optimize branch lengths (in coalescent units) of given input species tree #
#' #############################################################################
#' handle.OptimizedBranchLengths <- Optimize.BranchLengths(handle.SpeciesTree = handle.SpeciesTree,
#'                                                         handle.GeneTrees = handle.Simulated_GeneTrees,
#'                                                         numeric.Stells_algorthm = 1,
#'                                                         string.PathDir = '~/Desktop/')

##########################
# Optimize.BranchLengths #
##########################
Optimize.BranchLengths <- function(handle.SpeciesTree, handle.GeneTrees, numeric.Stells_algorthm, string.PathDir){

  #################################################################################
  # Define directory used for computing gene trees likelihoods given species tree #
  #################################################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_GeneTreeLikelihoods = paste(string.PathDir, '/OptimizeBranches_SpeciesTree_', Sys.Date(), sep = "")
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
  if (numeric.Stells_algorthm == 0){
    if (Sys.info()["sysname"] == "Darwin"){
      string.Command_STELLS <- "stells-v2-1-0-1-mac -g GeneTrees.tree -s SpeciesTree.tree -O"
    }
    if (Sys.info()["sysname"] == "Linux"){
      string.Command_STELLS <- "stells-v2-1-0-1-linux64 -g GeneTrees.tree -s SpeciesTree.tree -O"
    }
  }
  if (numeric.Stells_algorthm == 1){
    if (Sys.info()["sysname"] == "Darwin"){
      string.Command_STELLS <- "stells-v2-1-0-1-mac -g GeneTrees.tree -s SpeciesTree.tree"
    }
    if (Sys.info()["sysname"] == "Linux"){
      string.Command_STELLS <- "stells-v2-1-0-1-linux64 -g GeneTrees.tree -s SpeciesTree.tree"
    }
  }

  ##############
  # Run STELLS #
  ##############
  handle.STELLS_Results <- system(command = string.Command_STELLS, intern = T)

  #####################################################################
  # Read results from STELLS and find optimum LnL and newick topology #
  #####################################################################
  numeric.Grep_Line_LnL <- grep(pattern = "Best branch length finder gives: ", x = handle.STELLS_Results)
  string.Line_LnL <- handle.STELLS_Results[[numeric.Grep_Line_LnL]]
  string.Split_Line_LnL <- strsplit(x = string.Line_LnL, split = "Best branch length finder gives: ")[[1]][2]
  numeric.Optimized_LnL <- as.numeric(string.Split_Line_LnL)

  numeric.Grep_Line_Newick <- grep(pattern = "The newick format of the inferred species tree with branch length", x = handle.STELLS_Results)
  string.Line_Newick <- handle.STELLS_Results[[numeric.Grep_Line_Newick]]
  string.Split_Line_Newick <- strsplit(x = string.Line_Newick, split = "The newick format of the inferred species tree with branch length: ")[[1]][2]
  string.Split_Line_Newick <- paste0(string.Split_Line_Newick, ';')
  handle.Optimized_SpeciesTree <- read.tree(text = string.Split_Line_Newick)

  #####################################################################
  # Return to directory and return the LnL and optimized species tree #
  #####################################################################
  setwd(string.CurrentDir)
  return(list(handle.Optimized_SpeciesTree = handle.Optimized_SpeciesTree, numeric.MaximizedLnL = numeric.Optimized_LnL))

}
