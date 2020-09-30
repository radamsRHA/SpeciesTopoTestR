#' Compute.GeneTree_Likelihoods_SpeciesNetwork: function to compute gene tree likelihoods using PhyloNet for a species topology and/or network
#'
#' This function returns a vector of gene tree likelihoods given a set of gene trees and a particular species tree (networks included)
#' @param handle.GeneTrees Phylo object of gene trees
#' @param string.SpeciesNetwork String of the species network in Rich newick format (can be read by dendroscope)
#' @param string.PathDir String of the path to a parent directory used for simulating gene trees
#' @keywords Species tree, multispecies coalescent, phylogenetics, phylogenomics
#' @return vector.GeneTreeProbs Vector of gene tree likelihoods for each gene tree provided
#' @export
#' @examples
#'
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
#' string.SpeciesNetwork_2 <- "((A:1,B:1):1,(C:1,D:1):1);"
#'
#'
#' ####################################################
#' # Simlate a set of gene trees for this species tree #
#' #####################################################
#' handle.SimulatedGeneTrees <- Simulate.GeneTrees_From_SpeciesNetwork(string.SpeciesNetwork = string.SpeciesNetwork,
#'                                                                     string.PathDir = '~/Desktop/',
#'                                                                     numeric.NumberOfGeneTrees = 100)
#'
#' ############################
#' # Optimize species network #
#' ############################
#' handle.Optimized_Network1 <- Optimize.Network(string.SpeciesNetwork = string.SpeciesNetwork,
#'                                               handle.GeneTrees = handle.SimulatedGeneTrees,
#'                                               string.PathDir = '~/Desktop/')
#'
#'
#' #################################################################
#' # Compute gene tree likelihoods given optimized species network #
#' #################################################################
#' vector.GeneTreeLikelihoods_SpeciesNetwork <- Compute.GeneTree_Likelihoods_SpeciesNetwork(string.SpeciesNetwork = handle.Optimized_Network1$string.Optimized_SpeciesNetwork,
#'                                                                                          handle.GeneTrees = handle.SimulatedGeneTrees,
#'                                                                                          string.PathDir = '~/Desktop/')
#'

###############################################
# Compute.GeneTree_Likelihoods_SpeciesNetwork #
###############################################
Compute.GeneTree_Likelihoods_SpeciesNetwork <- function(string.SpeciesNetwork, handle.GeneTrees, string.PathDir){

  #################################################################################
  # Define directory used for computing gene trees likelihoods given species tree #
  #################################################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_GeneTreeLikelihoods_Network = paste(string.PathDir, '/GeneTreeLikelihoods_Network_', Sys.Date(), sep = "")
  unlink(string.Path_Directory_GeneTreeLikelihoods_Network, recursive = T)
  dir.create(string.Path_Directory_GeneTreeLikelihoods_Network, showWarnings = T, recursive = T)

  ###################
  # Summarize input #
  ###################
  numeric.NumberOfGeneTrees <- length(handle.GeneTrees)
  vector.GeneTreeProbs <- rep(NA, numeric.NumberOfGeneTrees)

  ###########################
  # Loop through gene trees #
  ###########################
  for (i in 1:numeric.NumberOfGeneTrees){

    ################################
    # Extract individual gene tree #
    ################################
    handle.GeneTree_i <- handle.GeneTrees[[i]]
    handle.GeneTree_i$edge.length <- NULL
    string.Path_Directory_GeneTree_Like_i = paste0(string.Path_Directory_GeneTreeLikelihoods_Network, '/GeneTree_', i)
    unlink(string.Path_Directory_GeneTree_Like_i, recursive = T)
    dir.create(string.Path_Directory_GeneTree_Like_i, showWarnings = T, recursive = T)
    setwd(string.Path_Directory_GeneTree_Like_i)

    ###################################
    # Path to nexus file for Phylonet #
    ###################################
    string.NexusFile <- paste0(string.Path_Directory_GeneTree_Like_i, '/ComputeGeneTreeLike_SpeciesNetwork.nex')
    write(x = "#NEXUS\nBEGIN NETWORKS;", file = string.NexusFile)
    string.NetworkNexusCmd <- "Network net = XXX"
    string.NetworkNexusCmd <- gsub(pattern = "XXX", replacement = string.SpeciesNetwork, x = string.NetworkNexusCmd)
    write(x = string.NetworkNexusCmd, file = string.NexusFile, append = T)
    write(x = "END;\n", file = string.NexusFile, append = T)
    write(x = "BEGIN TREES;\n", file = string.NexusFile, append = T)

    ###################
    # Write to string #
    ###################
    string.WriteTreeString <- "Tree geneTreeXXX = YYY"
    string.WriteTreeString <- gsub(pattern = "XXX", replacement = i, x = string.WriteTreeString)
    string.WriteTreeString <- gsub(pattern = "YYY", replacement = write.tree(phy = handle.GeneTree_i, file = ""), x = string.WriteTreeString)
    write(x = string.WriteTreeString, file = string.NexusFile, append = T)

    string.CalcGTProb <- "CalGTProb net (geneTreeXXX) -m ac;"
    string.CalcGTProb <- gsub(pattern = "XXX", replacement = i, x = string.CalcGTProb)
    write(x = "END;\n", file = string.NexusFile, append = T)
    write(x = "BEGIN PHYLONET;\n", file = string.NexusFile, append = T)
    write(string.CalcGTProb, file = string.NexusFile, append = T)
    write(x = "END;\n", file = string.NexusFile, append = T)

    ################
    # Run phylonet #
    ################
    string.Path_PhyloNet_jar <- system.file("extdata", "PhyloNet_3.8.0.jar", package="SpeciesTopoTestR")
    file.copy(from = string.Path_PhyloNet_jar, to = paste0(string.Path_Directory_GeneTree_Like_i, '/PhyloNet_3.8.0.jar'))
    setwd(dir = string.Path_Directory_GeneTree_Like_i)
    handle.PhyloNet_Results <- system(command = "java -jar PhyloNet_3.8.0.jar ComputeGeneTreeLike_SpeciesNetwork.nex > GeneTreeLikelihood.txt")

    #####################
    # Read results file #
    #####################
    handle.ResultsFile <- paste0(string.Path_Directory_GeneTree_Like_i, '/GeneTreeLikelihood.txt')
    handle.ResultsFile <- readLines(handle.ResultsFile)
    string.ML <- handle.ResultsFile[5]
    string.ML <- strsplit(x = string.ML, split = "-", fixed = T)[[1]][2]
    numeric.ML <- -1*as.numeric(string.ML)
    vector.GeneTreeProbs[i] <- numeric.ML


  }
  ##########################################
  # Return to directory and return results #
  ##########################################
  setwd(string.CurrentDir)
  return(vector.GeneTreeProbs)
}
