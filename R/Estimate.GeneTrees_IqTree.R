#' Estimate.GeneTrees_IqTree: function to estimate a set of gene trees for an input list of nexus alignments (simulated with Simulate.Gene.Alignments)
#'
#' This function returns (1) a list of ML gene tree estimates obtained via IQ-Tree
#' @param list.Simulated.Nexus.Alignments List of nexus alignments that have been simulated using Simulate.Gene.Alignments
#' @param path.PathParentDir Path to parent directory used for gene tree estimation
#' @param string.Model String to specify substitution model for IQ-tree
#' @keywords empirical bayesian estimation, species tree inference, multispecies coalescent
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
#' ########################################################
#' # Define species topologies for simulation and testing #
#' ########################################################
#' handle.SpeciesTree1 <- read.tree(text = "(A:3,(B:2,(C:1,D:1):1):1);")
#' handle.SpeciesTree2 <- read.tree(text = "(A:3,(D:2,(C:1,B:1):1):1);")
#'
#'
#' ##########################
#' # Simulate gene tree set #
#' ##########################
#' handle.Simulated_GeneTrees_T1 <- Simulate.GeneTrees_From_SpeciesTree(handle.SpeciesTree = handle.SpeciesTree1,
#'                                                                      string.PathDir = '~/Desktop/',
#'                                                                      numeric.NumberOfGeneTrees = 100)
#' ############################################
#' # Convert to mutation units branch lengths #
#' ############################################
#' handle.Converted_GeneTrees_T1 <- list()
#' class(handle.Converted_GeneTrees_T1) <- "multiPhylo"
#'
#' for (i in 1:length(handle.Simulated_GeneTrees_T1)){
#'
#'   handle.Simulated_GeneTrees_T1_i <- handle.Simulated_GeneTrees_T1[[i]]
#'   handle.Simulated_GeneTrees_T1_i$edge.length = handle.Simulated_GeneTrees_T1_i$edge.length* 0.00005
#'   handle.Converted_GeneTrees_T1[[i]] <- handle.Simulated_GeneTrees_T1_i
#' }
#'
#' ################################
#' # Simulate sequence alignments #
#' ################################
#' vector.pi <- c(0.3, 0.2, 0.2, 0.3)
#' names(vector.pi) <- c("A", "C", "G", "T")
#' list.SimulatedAlignments <- Simulate.GeneSequenceAlignments(list.Simulated.Gene.Tree.Set = handle.Converted_GeneTrees_T1,
#'                                                             path.PathParentDir = '~/Desktop/',
#'                                                             numeric.Ratio = 4.6,
#'                                                             numeric.LocusLength = 1000,
#'                                                             vector.BaseFrequencies = vector.pi)
#'
#' #######################
#' # Estimate gene trees #
#' #######################
#' handle.EstimatedGeneTrees <- Estimate.GeneTrees_IqTree(list.Simulated.Nexus.Alignments = list.SimulatedAlignments$list.Simulated.Nexus.Alignments,
#'                                                        path.PathParentDir = '~/Desktop/',
#'                                                        string.Model = "HKY")
#'
#' ################
#' # Conduct KH_2 #
#' ################
#' handle.Topologies <- list()
#' class(handle.Topologies) <- "multiPhylo"
#' handle.Topologies[[1]] <- handle.SpeciesTree1
#' handle.Topologies[[2]] <- handle.SpeciesTree2
#'
#' Conduct.SpeciesTopoTestR(handle.Topologies2Test = handle.Topologies,
#'                          handle.GeneTrees = handle.EstimatedGeneTrees$list.MLE.GeneTrees,
#'                          numeric.NumberOfReps = 1000,
#'                          string.Test = "KH",
#'                          numeric.Algorithm = 2,
#'                          boo.Networks = F,
#'                          string.PathDir = '~/Desktop/')

##############################################################################################################
# Estimate.GeneTrees_IqTree: function to estimate a set of gene trees given an input set of nexus alignments #
##############################################################################################################
Estimate.GeneTrees_IqTree <- function(list.Simulated.Nexus.Alignments, path.PathParentDir, string.Model){

  ##############################
  # Summarize input alignments #
  ##############################
  numeric.NumberOfNexusAlignments <- length(list.Simulated.Nexus.Alignments)
  list.MLE.GeneTrees <- list()
  class(list.MLE.GeneTrees) <- "multiPhylo"

  #########################################################
  # Set top-level parent dir for all gene tree estimation #
  #########################################################
  string.Estimation.Parent.Dir = paste(path.PathParentDir, '/Estimate_GeneTrees_IQTree', Sys.Date(), sep = "")
  unlink(string.Estimation.Parent.Dir, recursive = T)
  dir.create(string.Estimation.Parent.Dir, showWarnings = T, recursive = T)

  #########################################################
  # Loop through nexus alignments and estimate gene trees #
  #########################################################
  for (i in 1:numeric.NumberOfNexusAlignments){

    ##############################################
    # Make subdirectory for gene tree estimation #
    ##############################################
    string.Estimatation.GeneTree.Dir = paste(string.Estimation.Parent.Dir, '/GeneTree_', i, "_Estimation", sep = "")
    unlink(string.Estimatation.GeneTree.Dir, recursive = T)
    dir.create(string.Estimatation.GeneTree.Dir, showWarnings = T, recursive = T)

    ##########################################################
    # Extract nexus alignment and write to file in directory #
    ##########################################################
    handle.NexusAlignment <- list.Simulated.Nexus.Alignments[[i]]
    string.Path.To.Fasta.File <- paste(string.Estimatation.GeneTree.Dir, '/Alignment_GeneTree_', i, ".fasta", sep = "")
    write.FASTA(x = handle.NexusAlignment, file = string.Path.To.Fasta.File)

    string.IQ.Command <- "iqtree -s FASTA -m MMM -redo -t RANDOM -safe"
    string.IQ.Command <- gsub(pattern = "FASTA", replacement = paste('Alignment_GeneTree_', i, ".fasta", sep = ""), x = string.IQ.Command)
    string.IQ.Command <- gsub(pattern = "MMM", replacement = string.Model, x = string.IQ.Command)

    ###########################################
    # Change dir and estimate tree with RAXML #
    ###########################################
    setwd(string.Estimatation.GeneTree.Dir)
    system(command = string.IQ.Command, intern = T)

    ##################################
    # Save specific commands to file #
    ##################################
    write(x = string.IQ.Command, file = paste(string.Estimatation.GeneTree.Dir, '/IQCommand.txt', sep = ""))

    #############################
    # Read best ML tree from IQ #
    #############################
    string.Path.To.ML.GeneTree <- paste(string.Estimatation.GeneTree.Dir, '/Alignment_GeneTree_', i, '.fasta.treefile', sep = "")
    handle.MLE.GeneTree <- read.tree(file = string.Path.To.ML.GeneTree)
    list.MLE.GeneTrees[[i]] <- handle.MLE.GeneTree

  }

  ##################
  # Return results #
  ##################

  return(list(list.MLE.GeneTrees = list.MLE.GeneTrees))

}
