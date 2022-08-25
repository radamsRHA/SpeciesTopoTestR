#' Simulate.GeneTrees_From_SpeciesNetwork: function to simulate a set of gene trees given a species network using PhyloNet under the MSC
#'
#' This function returns a list of simulated gene trees that have been generated for a given input species network
#' @param string.SpeciesNetwork String of the species network in Rich newick format (can be read by dendroscope)
#' @param string.PathDir String of the path to a parent directory that will be used for simulating gene trees
#' @param numeric.NumberOfGeneTrees Numeric number of gene trees to simulate
#' @keywords Species tree, multispecies coalescent, phylogenetics, phylogenomics
#' @return handle.SimulatedGeneTrees List of simulated gene trees
#' @export
#' @examples
#'
#' ################
#' # Load depends #
#' ################
#' library(SpeciesTopoTestR)
#' library(ape)
#' 
#' ####################################
#' # Generate example species network #
#' ####################################
#' string.SpeciesNetwork <- "(((((C:1.0,D:1.0):1)#H1:0::0.25,A:1.0):2,B:1.0):2,#H1:0::0.75);"
#' 
#' ####################################################
#' # Simlate a set of gene trees for this species tree #
#' #####################################################
#' handle.SimulatedGeneTrees <- Simulate.GeneTrees_From_SpeciesNetwork(string.SpeciesNetwork = string.SpeciesNetwork,
#'                                                                     string.PathDir = '~/Desktop/',
#'                                                                     numeric.NumberOfGeneTrees = 10)

##########################################
# Simulate.GeneTrees_From_SpeciesNetwork #
##########################################
Simulate.GeneTrees_From_SpeciesNetwork <- function(string.SpeciesNetwork, numeric.NumberOfGeneTrees, string.PathDir){
  
  ########################################################################
  # Define directory used for simulating gene trees from species network #
  ########################################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_GeneTreeSimulation = paste(string.PathDir, '/GeneTreeSimulation_Networks', Sys.Date(), sep = "")
  unlink(string.Path_Directory_GeneTreeSimulation, recursive = T)
  dir.create(string.Path_Directory_GeneTreeSimulation, showWarnings = T, recursive = T)
  setwd(string.Path_Directory_GeneTreeSimulation)
  
  ###################################
  # Path to nexus file for Phylonet #
  ###################################
  string.NexusFile <- paste0(string.Path_Directory_GeneTreeSimulation, '/SimulatedGeneTrees_Network.nex')
  write(x = "#NEXUS\nBEGIN NETWORKS;", file = string.NexusFile)
  string.NetworkNexusCmd <- "Network net = XXX"
  string.NetworkNexusCmd <- gsub(pattern = "XXX", replacement = string.SpeciesNetwork, x = string.NetworkNexusCmd)
  write(x = string.NetworkNexusCmd, file = string.NexusFile, append = T)
  
  string.NumberGeneTreesCmd <- gsub(pattern = "XXX", replacement = numeric.NumberOfGeneTrees, x = "simGTinNetwork net XXX;")
  write(x = "END;\n\nBEGIN PHYLONET;", file = string.NexusFile, append = T)
  write(x = string.NumberGeneTreesCmd, file = string.NexusFile, append = T)
  write(x = "END;", file = string.NexusFile, append = T)
  
  ################
  # Run phylonet #
  ################
  string.Path_PhyloNet_jar <- system.file("extdata", "PhyloNet_3.8.0.jar", package="SpeciesTopoTestR")
  file.copy(from = string.Path_PhyloNet_jar, to = paste0(string.Path_Directory_GeneTreeSimulation, '/PhyloNet_3.8.0.jar'))
  setwd(dir = string.Path_Directory_GeneTreeSimulation)
  handle.PhyloNet_Results <- system(command = "java -jar PhyloNet_3.8.0.jar SimulatedGeneTrees_Network.nex > SimulatedGeneTrees.tree")
  setwd(dir = string.CurrentDir)
  
  ################################
  # Read in simulated gene trees #
  ################################
  handle.SimulatedGeneTrees <- read.tree(file = paste0(string.Path_Directory_GeneTreeSimulation, '/SimulatedGeneTrees.tree'))
  
  return(handle.SimulatedGeneTrees)
  
}