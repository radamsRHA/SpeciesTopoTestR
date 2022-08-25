#' Simulate.GeneSequenceAlignments: function to simulate a set of sequence alignments for a list of gene trees
#'
#' This function returns (1) a list of nexus alignments that have been simulated for the input gene trees
#' @param list.Simulated.Gene.Tree.Set Multiphylo list containing all the gene trees
#' @param path.PathParentDir Path to parent directory used for simulations
#' @param numeric.Ratio Kappa parameter for HKY substitution model
#' @param numeric.LocusLength Length of each locus
#' @param vector.BaseFrequencies Vector of base equibrium frequencies, with names equal to the nucleotide state
#' @keywords empirical bayesian estimation, species tree inference, multispecies coalescent
#' @export
#' @examples
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
#'                                                                      numeric.NumberOfGeneTrees = 3)
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
#'                                                             numeric.LocusLength = 100,
#'                                                             vector.BaseFrequencies = vector.pi)

############################################################################################################################
# Simulate.GeneSequenceAlignments: function to simulate a list of gene sequence alignments for a set of input gene trees #
############################################################################################################################
Simulate.GeneSequenceAlignments <- function(list.Simulated.Gene.Tree.Set, path.PathParentDir, numeric.Ratio, numeric.LocusLength, vector.BaseFrequencies){

  ##############################
  # Summarize input gene trees #
  ##############################
  numeric.NumberOfGeneTrees <- length(list.Simulated.Gene.Tree.Set)
  list.Simulated.Nexus.Alignments <- list()

  ##############################################################
  # Set top-level parent dir for all gene sequence simulations #
  ##############################################################
  string.Simulation.Parent.Dir = paste(path.PathParentDir, '/Simulate_GeneSequenceAlignments', Sys.Date(), sep = "")
  unlink(string.Simulation.Parent.Dir, recursive = T)
  dir.create(string.Simulation.Parent.Dir, showWarnings = T, recursive = T)

  ##################################################
  # Loop through gene trees and simulate alignment #
  ##################################################
  for (i in 1:numeric.NumberOfGeneTrees){

    #########################################################
    # Make subdirectory for gene tree (required by Seq-Gen) #
    #########################################################
    string.Simulation.Gene.Dir = paste(string.Simulation.Parent.Dir, '/Gene_', i, "_Simulation", sep = "")
    unlink(string.Simulation.Gene.Dir, recursive = T)
    dir.create(string.Simulation.Gene.Dir, showWarnings = T, recursive = T)

    ####################################################
    # Extract gene tree and write to file in directory #
    ####################################################
    handle.Gene.Tree.i <- list.Simulated.Gene.Tree.Set[[i]]
    string.Path.To.Gene.Tree.File <- paste(string.Simulation.Gene.Dir, '/GeneTree_', i, ".tree", sep = "")
    write.tree(phy = handle.Gene.Tree.i, file = string.Path.To.Gene.Tree.File)

    ######################################
    # Specify command string for seq-gen #
    ######################################
    string.SeqGen.Command <- "seq-gen -mHKY -l LENGTH -on -t TRATIO -f AAAA,CCCC,GGGG,TTTT GENETREE > AlignmentFile"
    string.SeqGen.Command <- gsub(pattern = "LENGTH", replacement = numeric.LocusLength, x = string.SeqGen.Command)
    string.SeqGen.Command <- gsub(pattern = "TRATIO", replacement = numeric.Ratio, x = string.SeqGen.Command)
    string.SeqGen.Command <- gsub(pattern = "GENETREE", replacement = paste('GeneTree_', i, ".tree", sep = ""), x = string.SeqGen.Command)
    string.SeqGen.Command <- gsub(pattern = "AlignmentFile", replacement = paste('GeneAlignment_', i, ".nex", sep = ""), x = string.SeqGen.Command)
    string.SeqGen.Command <- gsub(pattern = "AAAA", replacement = vector.BaseFrequencies[names(vector.BaseFrequencies) == "A"], x = string.SeqGen.Command)
    string.SeqGen.Command <- gsub(pattern = "CCCC", replacement = vector.BaseFrequencies[names(vector.BaseFrequencies) == "C"], x = string.SeqGen.Command)
    string.SeqGen.Command <- gsub(pattern = "GGGG", replacement = vector.BaseFrequencies[names(vector.BaseFrequencies) == "G"], x = string.SeqGen.Command)
    string.SeqGen.Command <- gsub(pattern = "TTTT", replacement = vector.BaseFrequencies[names(vector.BaseFrequencies) == "T"], x = string.SeqGen.Command)
    print(string.SeqGen.Command)

    #############################################
    # change directory and simulat with seq-gen #
    #############################################
    setwd(dir = string.Simulation.Gene.Dir)
    handle.SeqGen.Output <- system(command = string.SeqGen.Command, intern = F)

    ##################################
    # Save specific commands to file #
    ##################################
    write(x = string.SeqGen.Command, file = paste(string.Simulation.Gene.Dir, '/SeqGenCommand.txt', sep = ""))

    ###################
    # Read nexus file #
    ###################
    handle.Simulated.Nexus.File <- as.DNAbin(read.nexus.data(file = paste(string.Simulation.Gene.Dir, '/', 'GeneAlignment_', i, ".nex", sep = "")))
    list.Simulated.Nexus.Alignments[[i]] <- handle.Simulated.Nexus.File


  }
  ##################
  # Return results #
  ##################
  return(list(list.Simulated.Nexus.Alignments = list.Simulated.Nexus.Alignments))

}
