#' Simulate.GeneTrees_From_SpeciesTree3: function to simulate a set of gene trees given a species tree using HYBRID-LAMBDA under the MSC
#'
#' This function returns a list of simulated gene trees that have been generated for a given input species tree
#' @param handle.SpeciesTree Phylo object of the input species tree
#' @param string.PathDir String of the path to a parent directory that will be used for simulating gene trees
#' @param numeric.NumberOfGeneTrees Numeric number of gene trees to simulate
#' @keywords Species tree, multispecies coalescent, phylogenetics, phylogenomics
#' @return handle.SimulatedGeneTrees List of simulated gene trees
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
#' Simulate.GeneTrees_From_SpeciesTree(handle.SpeciesTree = handle.SpeciesTree,
#'                                     string.PathDir = '~/Desktop/',
#'                                     numeric.NumberOfGeneTrees = 10)


#######################################
# Simulate.GeneTrees_From_SpeciesTree3 #
#######################################
Simulate.GeneTrees_From_SpeciesTree3 <- function(handle.SpeciesTree, numeric.NumberOfGeneTrees, string.PathDir){

  #####################################################################
  # Define directory used for simulating gene trees from species tree #
  #####################################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_GeneTreeSimulation = paste(string.PathDir, '/GeneTreeSimulation_', Sys.Date(), sep = "")
  unlink(string.Path_Directory_GeneTreeSimulation, recursive = T)
  dir.create(string.Path_Directory_GeneTreeSimulation, showWarnings = T, recursive = T)
  setwd(string.Path_Directory_GeneTreeSimulation)

  ######################################################
  # Command for simulating gene trees with hybrid-coal #
  ######################################################
  string.Command_HybridLambda <- "hybrid-Lambda -spcu 'XXX' -num YYY -seed ZZZ -pop 1"
  string.Command_HybridLambda <- gsub(pattern = "XXX", replacement = write.tree(phy = handle.SpeciesTree), x = string.Command_HybridLambda)
  string.Command_HybridLambda <- gsub(pattern = "YYY", replacement = numeric.NumberOfGeneTrees, x = string.Command_HybridLambda)
  string.Command_HybridLambda <- gsub(pattern = "ZZZ", replacement = round(runif(1, 0, 10^10), digits = 0), x = string.Command_HybridLambda)
  #print(string.Command_HybridLambda)
  system(command = string.Command_HybridLambda, intern = T)

  #############################
  # Read simulated gene trees #
  #############################
  string.Path_To_SimulatedGeneTrees <- paste0(string.Path_Directory_GeneTreeSimulation, '/OUT_coal_unit')

  ################################
  # Remove underscores from name #
  ################################
  string.Command_SED <- "sed -i '' 's|_1||g' OUT_coal_unit"
  system(command = string.Command_SED, intern = T)

  #################################
  # Read simulated gene tree file #
  #################################
  handle.SimulatedGeneTrees <- read.tree(string.Path_To_SimulatedGeneTrees)

  #####################################################################
  # Return to directory and return results of simulated gene tree set #
  #####################################################################
  setwd(string.CurrentDir)
  return(handle.SimulatedGeneTrees)

}
