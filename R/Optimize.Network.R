#' Optimize.Network: function to optimize a species network topology using Phylonet
#'
#' This function returns a species network (string) with branch lengths amd hybyridization edges that have been optimized using the MSC and the PhyloNet algorithm
#' @param handle.GeneTrees Phylo object containing a list of gene trees
#' @param string.SpeciesNetwork String of the species network in Rich newick format (can be read by dendroscope)
#' @param string.PathDir String defining the path to a parent directory used for optimizing species tree
#' @keywords Species tree, multispecies coalescent, phylogenetics, phylogenomics
#' @return string.Optimized_SpeciesNetwork Species network with branch lengths (in coalescent units) and hybridization edge that has been optimized with PhyloNet
#' @return numeric.MaximizedLnL Numeric containing the maximum likelihood given the optimized parameters of the species network
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
#'
#' ###############################################
#' # Optimize network using simulated gene trees #
#' ###############################################
#' handle.Optimized_Network <- Optimize.Network(string.SpeciesNetwork = string.SpeciesNetwork,
#'                                              handle.GeneTrees = handle.SimulatedGeneTrees,
#'                                              string.PathDir = '~/Desktop/')
#'

####################
# Optimize.Network #
####################
Optimize.Network <- function(string.SpeciesNetwork, handle.GeneTrees, string.PathDir){

  #################################################################################
  # Define directory used for computing gene trees likelihoods given species tree #
  #################################################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_GeneTreeLikelihoods = paste(string.PathDir, '/Optimize_SpeciesNetwork_', Sys.Date(), sep = "")
  unlink(string.Path_Directory_GeneTreeLikelihoods, recursive = T)
  dir.create(string.Path_Directory_GeneTreeLikelihoods, showWarnings = T, recursive = T)
  setwd(string.Path_Directory_GeneTreeLikelihoods)

  ###################
  # Summarize input #
  ###################
  numeric.NumberOfGeneTrees <- length(handle.GeneTrees)
  vector.GeneTreeProbs <- rep(NA, numeric.NumberOfGeneTrees)

  ###################################
  # Path to nexus file for Phylonet #
  ###################################
  string.NexusFile <- paste0(string.Path_Directory_GeneTreeLikelihoods, '/Optimized_SpeciesNetwork.nex')
  write(x = "#NEXUS\nBEGIN NETWORKS;", file = string.NexusFile)
  string.NetworkNexusCmd <- "Network net = XXX"
  string.NetworkNexusCmd <- gsub(pattern = "XXX", replacement = string.SpeciesNetwork, x = string.NetworkNexusCmd)
  write(x = string.NetworkNexusCmd, file = string.NexusFile, append = T)
  write(x = "END;\n", file = string.NexusFile, append = T)
  write(x = "BEGIN TREES;\n", file = string.NexusFile, append = T)

  ############################
  # Write gene trees to file #
  ############################
  string.CalcGTProb <- "CalGTProb net ("
  for (i in 1:numeric.NumberOfGeneTrees){
    string.WriteTreeString <- "Tree geneTreeXXX = YYY"
    string.WriteTreeString <- gsub(pattern = "XXX", replacement = i, x = string.WriteTreeString)
    handle.GeneTree_i <- handle.GeneTrees[[i]]
    handle.GeneTree_i$edge.length <- NULL
    string.WriteTreeString <- gsub(pattern = "YYY", replacement = write.tree(phy = handle.GeneTree_i, file = ""), x = string.WriteTreeString)
    write(x = string.WriteTreeString, file = string.NexusFile, append = T)
    string.CalcGTProb <- paste0(string.CalcGTProb, 'geneTree', i, ',')
  }

  string.CalcGTProb <- paste0(string.CalcGTProb, 'X')
  string.CalcGTProb <- gsub(pattern = ",X", replacement = ") -o;", x = string.CalcGTProb)
  write(x = "END;\n", file = string.NexusFile, append = T)
  write(x = "BEGIN PHYLONET;\n", file = string.NexusFile, append = T)
  write(string.CalcGTProb, file = string.NexusFile, append = T)
  write(x = "END;\n", file = string.NexusFile, append = T)

  ################
  # Run phylonet #
  ################
  string.Path_PhyloNet_jar <- system.file("extdata", "PhyloNet_3.8.0.jar", package="SpeciesTopoTestR")
  file.copy(from = string.Path_PhyloNet_jar, to = paste0(string.Path_Directory_GeneTreeLikelihoods, '/PhyloNet_3.8.0.jar'))
  setwd(dir = string.Path_Directory_GeneTreeLikelihoods)
  handle.PhyloNet_Results <- system(command = "java -jar PhyloNet_3.8.0.jar Optimized_SpeciesNetwork.nex > Optimized_Results.txt")

  #####################
  # Read results file #
  #####################
  handle.ResultsFile <- paste0(string.Path_Directory_GeneTreeLikelihoods, '/Optimized_Results.txt')
  handle.ResultsFile <- readLines(handle.ResultsFile)

  string.Optimized_SpeciesNetwork <- handle.ResultsFile[4]
  string.ML <- handle.ResultsFile[5]
  string.ML <- strsplit(x = string.ML, split = "-", fixed = T)[[1]][2]
  numeric.ML <- -1*as.numeric(string.ML)

  #####################################################################
  # Return to directory and return the LnL and optimized species tree #
  #####################################################################
  setwd(string.CurrentDir)
  return(list(string.Optimized_SpeciesNetwork = string.Optimized_SpeciesNetwork, numeric.MaximizedLnL = numeric.ML))

}
