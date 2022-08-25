#' Conduct.SH_2: function to conduct the SH_2 STAR test given a set of input (all plausible, including the ML) species topologies using the RELL method
#'
#' This function returns a list containing p-values of the SH_2 STAR test for a set of input topologies
#' @param handle.SpeciesTopologies List of class multiPhylo containing the plausible topologies to be tested
#' @param handle.GeneTrees Phylo object containing a list of the input gene trees
#' @param numeric.NumberOfReps Number of bootstrap replicates to analyze
#' @param string.PathDir String defining the path to a parent directory used for conduct SH_2 STAR test
#' @keywords Species tree, multispecies coalescent, phylogenetics, phylogenomics
#' @return List Returns a list containing (1) matrix.Change_Delta_BS_Results, (2) vector.Pvalues, and (3) vector.Observed_Delta_ML.
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
#' #################################
#' # Generate example species trees #
#' #################################
#' handle.SpeciesTree1 <- read.tree(text = "(A:1.0,(B:0.5,(C:1.0,D:1.0):1.5):0.5);")
#' handle.SpeciesTree1$edge.length <- handle.SpeciesTree1$edge.length * 0.1
#' handle.SpeciesTree2 <- read.tree(text = "(C:1.0,(B:0.5,(A:1.0,D:1.0):1.5):0.5);")
#' handle.SpeciesTree3 <- read.tree(text = "(D:1.0,(B:0.5,(A:1.0,C:1.0):1.5):0.5);")
#' handle.SpeciesTopologies <- list()
#' class(handle.SpeciesTopologies) <- "multiPhylo"
#' handle.SpeciesTopologies[[1]] <- handle.SpeciesTree3
#' handle.SpeciesTopologies[[2]] <- handle.SpeciesTree2
#' handle.SpeciesTopologies[[3]] <- handle.SpeciesTree1
#'
#'
#'
#' #################################################
#' # Simlate a set of gene trees for species tree1 #
#' #################################################
#' handle.Simulated_GeneTrees <- Simulate.GeneTrees_From_SpeciesTree(handle.SpeciesTree = handle.SpeciesTree1,
#'                                                                   string.PathDir = '~/Desktop/',
#'                                                                   numeric.NumberOfGeneTrees = 100)
#'
#' ################################
#' # Conduct SH_2 STAR test for S #
#' ################################
#' Conduct.SH_2(handle.SpeciesTopologies = handle.SpeciesTopologies,
#'              handle.InputGeneTrees = handle.Simulated_GeneTrees,
#'              numeric.NumberOfReps = 1000,
#'              string.PathDir = '~/Desktop/')
#'
#'

################
# Conduct.SH_2 #
################
Conduct.SH_2 <- function(handle.SpeciesTopologies, handle.InputGeneTrees, numeric.NumberOfReps, string.PathDir){

  ###################
  # Summarize input #
  ###################
  numeric.NumberOfSpeciesTopologies <- length(handle.SpeciesTopologies)
  numeric.NumberOfInputGeneTrees <- length(handle.InputGeneTrees)

  ####################################################
  # Define directory used for conduct SH_1 STAR test #
  ####################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_SH2 = paste(string.PathDir, '/Conduct.SH_2_', Sys.Date(), sep = "")
  unlink(string.Path_Directory_SH2, recursive = T)
  dir.create(string.Path_Directory_SH2, showWarnings = T, recursive = T)

  #######################################################
  # Define directory used for optimizing all topologies #
  #######################################################
  string.Path_Directory_OptimizeTopologies = paste0(string.Path_Directory_SH2, '/OptimizeTopologies/')
  unlink(string.Path_Directory_OptimizeTopologies, recursive = T)
  dir.create(string.Path_Directory_OptimizeTopologies, showWarnings = T, recursive = T)

  #######################################################
  # Define directory used for optimizing all topologies #
  #######################################################
  string.Path_Directory_ComputeGeneTreeLikelihoods = paste0(string.Path_Directory_SH2, '/Compute_GeneTreeLikelihoods/')
  unlink(string.Path_Directory_ComputeGeneTreeLikelihoods, recursive = T)
  dir.create(string.Path_Directory_ComputeGeneTreeLikelihoods, showWarnings = T, recursive = T)

  #############################################################
  # Define subdirectory used for bootstrap replicate analyses #
  #############################################################
  string.Path_Directory_BootStrapDir = paste(string.Path_Directory_SH2, '/BootStrapping', sep = "")
  unlink(string.Path_Directory_BootStrapDir, recursive = T)
  dir.create(string.Path_Directory_BootStrapDir, showWarnings = T, recursive = T)


  ###########################################################
  # Loop through species topologies and optimize parameters #
  ###########################################################
  vector.SpeciesTopologies_Maximizedikelihoods <- rep(NA, numeric.NumberOfSpeciesTopologies)
  names(vector.SpeciesTopologies_Maximizedikelihoods) <- paste0("Topology_", 1:numeric.NumberOfSpeciesTopologies)
  handle.OptimizedTopologies <- list()
  class(handle.OptimizedTopologies) <- "multiPhylo"
  for (i in 1:numeric.NumberOfSpeciesTopologies){

    ################################
    # Make subdir for optimization #
    ################################
    handle.Topology_i <- handle.SpeciesTopologies[[i]]
    string.Path_Directory_OptimizeTopology_i = paste0(string.Path_Directory_OptimizeTopologies, '/SpeciesTopology_', i)
    unlink(string.Path_Directory_OptimizeTopology_i, recursive = T)
    dir.create(string.Path_Directory_OptimizeTopology_i, showWarnings = T, recursive = T)

    ###########################################################
    # Optimize Topology i and append to vector of likelihoods #
    ###########################################################
    handle.Optimized_SpeciesTree_i <- Optimize.BranchLengths(handle.SpeciesTree = handle.Topology_i,
                                                             handle.GeneTrees = handle.InputGeneTrees,
                                                             string.PathDir = string.Path_Directory_OptimizeTopology_i,
                                                             numeric.Stells_algorthm = 1)
    vector.SpeciesTopologies_Maximizedikelihoods[i] <- handle.Optimized_SpeciesTree_i$numeric.MaximizedLnL
    handle.OptimizedTopologies[[i]] <- handle.Optimized_SpeciesTree_i$handle.Optimized_SpeciesTree
  }

  ######################################
  # Find ML from vector of likelihoods #
  ######################################
  numeric.ML_SpeciesTopology <- max(vector.SpeciesTopologies_Maximizedikelihoods)
  numeric.Position_ML_Topology <- names(vector.SpeciesTopologies_Maximizedikelihoods)[vector.SpeciesTopologies_Maximizedikelihoods==numeric.ML_SpeciesTopology]

  #############################################################
  # Compute gene tree likelihoods for each maximized topology #
  #############################################################
  matrix.GeneTreeLikelihoods_Topologies <- matrix(nrow = numeric.NumberOfInputGeneTrees, ncol = numeric.NumberOfSpeciesTopologies)
  colnames(matrix.GeneTreeLikelihoods_Topologies) <- paste0("Topology_", 1:numeric.NumberOfSpeciesTopologies)
  rownames(matrix.GeneTreeLikelihoods_Topologies) <- paste0("GeneTree_", 1:numeric.NumberOfInputGeneTrees)
  for (i in 1:numeric.NumberOfSpeciesTopologies){

    ##########################
    # Get optimized topology #
    ##########################
    handle.OptimizedSpeciesTree_i <- handle.OptimizedTopologies[[i]]

    ###################################################
    # Define directory used for optimizing topology i #
    ###################################################
    string.Path_Directory_ComputeGeneTreeLikelihoods_Topology_i = paste0(string.Path_Directory_ComputeGeneTreeLikelihoods, '/Topology_', i)
    unlink(string.Path_Directory_ComputeGeneTreeLikelihoods_Topology_i, recursive = T)
    dir.create(string.Path_Directory_ComputeGeneTreeLikelihoods_Topology_i, showWarnings = T, recursive = T)

    ###########################################
    # Compute likelihoods of input gene trees #
    ###########################################
    vector.GeneTree_Likelihoods_SpeciesTree_i <- Compute.GeneTree_Likelihoods(handle.SpeciesTree = handle.OptimizedSpeciesTree_i,
                                                                             handle.GeneTrees = handle.InputGeneTrees,
                                                                             string.PathDir = string.Path_Directory_ComputeGeneTreeLikelihoods_Topology_i)

    matrix.GeneTreeLikelihoods_Topologies[,i] <- vector.GeneTree_Likelihoods_SpeciesTree_i
  }

  #########################################################
  # Step 1: compute observed deltas using the ML topology #
  #########################################################
  vector.Observed_Delta_ML <- numeric.ML_SpeciesTopology-vector.SpeciesTopologies_Maximizedikelihoods

  ########################################
  # Step 2: Conduct BootStrap Resampling #
  ########################################
  matrix.Delta_BS_Results <- matrix(nrow = numeric.NumberOfReps, ncol = numeric.NumberOfSpeciesTopologies)
  colnames(matrix.Delta_BS_Results) <- paste0("Topology_", 1:numeric.NumberOfSpeciesTopologies)
  rownames(matrix.Delta_BS_Results) <- paste0("Rep_", 1:numeric.NumberOfReps)
  list.BoostrapReplicates_GeneTrees <- Resample.Bootstrap_GeneTree_Replicates_NP(handle.GeneTrees = handle.InputGeneTrees,
                                                                                 numeric.NumberOfReps = numeric.NumberOfReps)

  #####################################
  # Step 2.1: Loop through replicates #
  #####################################
  for (i in 1:numeric.NumberOfReps){

    ###############################
    # Extract bootstrap replicate #
    ###############################
    #print(gsub("Conducting BS replicate XXX...", pattern = "XXX", replacement = i))
    handle.BootStrapped_GeneTrees_i <- list.BoostrapReplicates_GeneTrees[[i]]
    vector.SpeciesTopologies_Maximizedikelihoods_j <- rep(NA, numeric.NumberOfSpeciesTopologies)
    names(vector.SpeciesTopologies_Maximizedikelihoods_j) <- 1:numeric.NumberOfSpeciesTopologies

    #############################################################
    # Define subdirectory used for bootstrap replicate analyses #
    #############################################################
    string.Path_Directory_BootStrapRep_i = paste(string.Path_Directory_BootStrapDir, '/Rep_', i, sep = "")
    unlink(string.Path_Directory_BootStrapRep_i, recursive = T)
    dir.create(string.Path_Directory_BootStrapRep_i, showWarnings = T, recursive = T)

    ###########################################
    # sample positions to compute likelihoods #
    ###########################################
    vector.Bootstrap_Positions <- sample(x = 1:numeric.NumberOfInputGeneTrees, size = numeric.NumberOfInputGeneTrees, replace = T)

    #####################################
    # Loop through plausible topologies #
    #####################################
    for (j in 1:numeric.NumberOfSpeciesTopologies){

      #######################################
      # Get vector of gene tree likelihoods #
      #######################################
      vector.GeneTreeLikelihoods_Topology_j <- matrix.GeneTreeLikelihoods_Topologies[,j]
      #handle.BootStrapped_GeneTreesLikes_SpeciesTree_j <- sample(x = vector.GeneTreeLikelihoods_Topology_j, size = numeric.NumberOfInputGeneTrees, replace = T)
      handle.BootStrapped_GeneTreesLikes_SpeciesTree_j <- vector.GeneTreeLikelihoods_Topology_j[vector.Bootstrap_Positions]
      vector.SpeciesTopologies_Maximizedikelihoods_j[j] <- sum(handle.BootStrapped_GeneTreesLikes_SpeciesTree_j)
    }

    ##################################
    # Append to matrix of BS results #
    ##################################
    matrix.Delta_BS_Results[i,] <- vector.SpeciesTopologies_Maximizedikelihoods_j

  }

  ##################################
  # Step 3: Adjust the likelihoods #
  ##################################
  matrix.Adjusted_Likelihoods <- matrix.Delta_BS_Results
  for (j in 1:numeric.NumberOfSpeciesTopologies){
    matrix.Adjusted_Likelihoods[,j] <- matrix.Delta_BS_Results[,j] - mean(matrix.Delta_BS_Results[,j])
  }

  ###################################################
  # Step 4: Find MLEs of each bootstrap replicate i #
  ###################################################
  matrix.Change_Delta_BS_Results <- matrix(nrow = numeric.NumberOfReps, ncol = numeric.NumberOfSpeciesTopologies)
  colnames(matrix.Change_Delta_BS_Results) <- paste0("Topology_", 1:numeric.NumberOfSpeciesTopologies)
  rownames(matrix.Change_Delta_BS_Results) <- paste0("Rep_", 1:numeric.NumberOfReps)

  for (i in 1:numeric.NumberOfReps){

    ##########################################
    # Extract LnLs for bootstrap replicate i #
    ##########################################
    vector.LnLs_BootstrapReplicate_i <- matrix.Adjusted_Likelihoods[i,]
    names(vector.LnLs_BootstrapReplicate_i) <- colnames(matrix.Adjusted_Likelihoods)
    numeric.ML_SpeciesTopology_i <- max(vector.LnLs_BootstrapReplicate_i)
    numeric.Position_ML_Topology_i <- names(vector.LnLs_BootstrapReplicate_i)[vector.LnLs_BootstrapReplicate_i==numeric.ML_SpeciesTopology_i]
    vector.BS_Delta_ML_i <- numeric.ML_SpeciesTopology_i-vector.LnLs_BootstrapReplicate_i
    matrix.Change_Delta_BS_Results[i,] <- vector.BS_Delta_ML_i

  }

  ############################
  # Step 5: Compute p-values #
  ############################
  vector.Pvalues <- rep(NA, numeric.NumberOfSpeciesTopologies)
  for (j in 1:numeric.NumberOfSpeciesTopologies){

    numeric.ObservedDelta_j <- vector.Observed_Delta_ML[j]
    vetor.BS_Deltas_j <- matrix.Change_Delta_BS_Results[,j]
    numeric.pValue <- length(vetor.BS_Deltas_j[vetor.BS_Deltas_j>=numeric.ObservedDelta_j])/numeric.NumberOfReps
    vector.Pvalues[j] <- numeric.pValue
  }

  ##########################################
  # Return to directory and return results #
  ##########################################
  setwd(dir = string.CurrentDir)
  names(vector.Pvalues) <- names(vector.Observed_Delta_ML)
  return(list(matrix.Change_Delta_BS_Results = matrix.Change_Delta_BS_Results,
              vector.Pvalues = vector.Pvalues,
              vector.Observed_Delta_ML = vector.Observed_Delta_ML))


}
