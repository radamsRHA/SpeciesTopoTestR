#' Optimize.NetworkSearch: function to find the ML network topology using PhyloNet
#'
#' This function returns an optimized network topology and its maximized likelihood
#' @param handle.GeneTrees Phylo object containing a list of gene trees
#' @param numeric.Stells_algorthm Numeric defining with algorithm of STELLS to use: 0 or 1
#' @param numeric.MaxReticulations Numeric defining the maximum number of reticulations for PhyloNet search
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
#' # Generate example species trees #
#' #################################
#' string.SpeciesNetwork <- "(((((C:1.0,D:1.0):1)#H1:0::0.5,A:1.0):2,B:1.0):2,#H1:0::0.5);"
#' string.SpeciesNetwork_2 <- "(((((C:1.0,A:1.0):1)#H1:0::0.5,D:1.0):2,B:1.0):2,#H1:0::0.5);"
#'
#' ##########################
#' # Simulate gene tree set #
#' ##########################
#' handle.Simulated_GeneTreeSet <- Simulate.GeneTrees_From_SpeciesNetwork(string.SpeciesNetwork = string.SpeciesNetwork,
#'                                                                        string.PathDir = '~/Desktop/',
#'                                                                        numeric.NumberOfGeneTrees = 5)
#'
#' ####################################
#' # Optimize network topology search #
#' ####################################
#' handle.OptimizedNetwork_Topology <- Optimize.NetworkSearch(handle.GeneTrees = handle.Simulated_GeneTreeSet,
#'                                                            numeric.MaxReticulations = 1,
#'                                                            string.PathDir = '~/Desktop/')
#'

##########################
# Optimize.NetworkSearch #
##########################
Optimize.NetworkSearch <- function(handle.GeneTrees, numeric.MaxReticulations, string.PathDir){

  ###################
  # Summarize input #
  ###################
  numeric.NumberOfGeneTrees <- length(handle.GeneTrees)

  #################################################################################
  # Define directory used for computing gene trees likelihoods given species tree #
  #################################################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_Optimize_SpeciesNetwork_Topology = paste(string.PathDir, '/Optimize_SpeciesNetwork_Topology_', Sys.Date(), sep = "")
  unlink(string.Path_Directory_Optimize_SpeciesNetwork_Topology, recursive = T)
  dir.create(string.Path_Directory_Optimize_SpeciesNetwork_Topology, showWarnings = T, recursive = T)
  setwd(string.Path_Directory_Optimize_SpeciesNetwork_Topology)

  ###################################
  # Path to nexus file for Phylonet #
  ###################################
  string.NexusFile <- paste0(string.Path_Directory_Optimize_SpeciesNetwork_Topology, '/Optimized_SpeciesNetwork_Topology.nex')
  write(x = "#NEXUS\n", file = string.NexusFile)
  write(x = "BEGIN TREES;\n", file = string.NexusFile, append = T)

  ############################
  # Write gene trees to file #
  ############################
  string.InferML <- "InferNetwork_ML ("
  for (i in 1:numeric.NumberOfGeneTrees){
    string.WriteTreeString <- "Tree geneTreeXXX = YYY"
    string.WriteTreeString <- gsub(pattern = "XXX", replacement = i, x = string.WriteTreeString)
    handle.GeneTree_i <- handle.GeneTrees[[i]]
    handle.GeneTree_i$edge.length <- NULL
    string.WriteTreeString <- gsub(pattern = "YYY", replacement = write.tree(phy = handle.GeneTree_i, file = ""), x = string.WriteTreeString)
    write(x = string.WriteTreeString, file = string.NexusFile, append = T)
    string.InferML <- paste0(string.InferML, 'geneTree', i, ',')
  }

  string.InferML <- paste0(string.InferML, 'X')
  string.InferML <- gsub(pattern = ",X", replacement = paste0(") ", numeric.MaxReticulations, " -o;"), x = string.InferML)
  write(x = "END;\n", file = string.NexusFile, append = T)
  write(x = "BEGIN PHYLONET;\n", file = string.NexusFile, append = T)
  write(string.InferML, file = string.NexusFile, append = T)
  write(x = "END;\n", file = string.NexusFile, append = T)

  ################
  # Run phylonet #
  ################
  string.Path_PhyloNet_jar <- system.file("extdata", "PhyloNet_3.8.0.jar", package="SpeciesTopoTestR")
  file.copy(from = string.Path_PhyloNet_jar, to = paste0(string.Path_Directory_Optimize_SpeciesNetwork_Topology, '/PhyloNet_3.8.0.jar'))
  setwd(dir = string.Path_Directory_Optimize_SpeciesNetwork_Topology)
  handle.PhyloNet_Results <- system(command = "java -jar PhyloNet_3.8.0.jar Optimized_SpeciesNetwork_Topology.nex > Optimized_NetworkTopology_Results.txt")

  #####################
  # Read results file #
  #####################
  handle.ResultsFile <- paste0(string.Path_Directory_Optimize_SpeciesNetwork_Topology, '/Optimized_NetworkTopology_Results.txt')
  handle.ResultsFile <- readLines(handle.ResultsFile)

  string.Optimized_SpeciesNetwork <- handle.ResultsFile[length(handle.ResultsFile)-1]
  string.ML <- handle.ResultsFile[length(handle.ResultsFile)]
  string.ML <- strsplit(x = string.ML, split = "-", fixed = T)[[1]][2]
  numeric.ML <- -1*as.numeric(string.ML)

  #####################################################################
  # Return to directory and return the LnL and optimized species tree #
  #####################################################################
  setwd(string.CurrentDir)
  return(list(string.Optimized_SpeciesNetwork = string.Optimized_SpeciesNetwork, numeric.MaximizedLnL = numeric.ML))


}
