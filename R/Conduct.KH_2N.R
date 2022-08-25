#' Conduct.KH_2N: function to conduct the KH_2N STAR test given two distinct species topologies (one or both are networks) and a set of input gene trees
#'
#' This function returns a list containing p-values of the KH_2N STAR test for two input species tree topologies
#' @param string.SpeciesNetwork1 String defining the first species topology (can be network or bifurcating)
#' @param string.SpeciesNetwork2 String defining the first species topology (can be network or bifurcating)
#' @param handle.GeneTrees Phylo object containing a list of the input gene trees
#' @param numeric.NumberOfReps Number of bootstrap replicates to analyze
#' @param string.PathDir String defining the path to a parent directory used for conduct KH_1 STAR test
#' @keywords Species tree, multispecies coalescent, phylogenetics, phylogenomics
#' @return List Returns a list containing (1) twosided pvalue, (2) upper p-values, (3) lower p-values, and (4) a vector of the bootstrapped test statistics delta
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
#' string.SpeciesNetwork_2 <- "((A:1,B:1):1,(C:1,D:1):1);"
#'
#'
#' ####################################################
#' # Simlate a set of gene trees for this species tree #
#' #####################################################
#' handle.SimulatedGeneTrees <- Simulate.GeneTrees_From_SpeciesNetwork(string.SpeciesNetwork = string.SpeciesNetwork,
#'                                                                     string.PathDir = '~/Desktop/',
#'                                                                     numeric.NumberOfGeneTrees = 5)
#'
#' Conduct.KH_2N(string.SpeciesNetwork1 = string.SpeciesNetwork,
#'               string.SpeciesNetwork2 = string.SpeciesNetwork_2,
#'               handle.InputGeneTrees = handle.SimulatedGeneTrees,
#'               numeric.NumberOfReps = 3,
#'               string.PathDir = '~/Desktop/')
#'

#################
# Conduct.KH_2N #
#################
Conduct.KH_2N <- function(string.SpeciesNetwork1, string.SpeciesNetwork2, handle.InputGeneTrees, numeric.NumberOfReps, string.PathDir){


  ###################
  # Summarize input #
  ###################
  numeric.NumberOfGeneTrees <- length(handle.InputGeneTrees)

  #####################################################
  # Define directory used for conduct KH_2N STAR test #
  #####################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_KH2N = paste(string.PathDir, '/Conduct.KH_2N_', Sys.Date(), sep = "")
  unlink(string.Path_Directory_KH2N, recursive = T)
  dir.create(string.Path_Directory_KH2N, showWarnings = T, recursive = T)

  ####################################################################################
  # Define subdirectory used for comput likelihoods given optimized species network1 #
  ####################################################################################
  string.Path_Directory_KH2_GeneTreeLikes_SpeciesNetwork1 = paste(string.Path_Directory_KH2N, '/GeneTreeLikes_Optimized_SpeciesNetwork_1', sep = "")
  unlink(string.Path_Directory_KH2_GeneTreeLikes_SpeciesNetwork1, recursive = T)
  dir.create(string.Path_Directory_KH2_GeneTreeLikes_SpeciesNetwork1, showWarnings = T, recursive = T)

  ####################################################################################
  # Define subdirectory used for comput likelihoods given optimized species network2 #
  ####################################################################################
  string.Path_Directory_KH2_GeneTreeLikes_SpeciesNetwork2 = paste(string.Path_Directory_KH2N, '/GeneTreeLikes_Optimized_SpeciesNetwork_2', sep = "")
  unlink(string.Path_Directory_KH2_GeneTreeLikes_SpeciesNetwork2, recursive = T)
  dir.create(string.Path_Directory_KH2_GeneTreeLikes_SpeciesNetwork2, showWarnings = T, recursive = T)

  #########################################################
  # Define subdirectory used for optimizing species tree1 #
  #########################################################
  string.Path_Directory_KH2_SpeciesNetwork1 = paste(string.Path_Directory_KH2N, '/Optimized_SpeciesNetwork_1', sep = "")
  unlink(string.Path_Directory_KH2_SpeciesNetwork1, recursive = T)
  dir.create(string.Path_Directory_KH2_SpeciesNetwork1, showWarnings = T, recursive = T)

  ######################################################################################
  # Compute observed test statistic for the differences in LnLs for the two topologies #
  ######################################################################################
  setwd(dir = string.Path_Directory_KH2_SpeciesNetwork1)
  handle.Optimized_SpeciesNetwork1_OBSERVED <- Optimize.Network(string.SpeciesNetwork = string.SpeciesNetwork1,
                                                                   handle.GeneTrees = handle.InputGeneTrees,
                                                                   string.PathDir = string.Path_Directory_KH2_SpeciesNetwork1)
  numeric.LnL_SpeciesNetwork1_OBSERVED <- handle.Optimized_SpeciesNetwork1_OBSERVED$numeric.MaximizedLnL
  string.Optimized_Network_1 <- handle.Optimized_SpeciesNetwork1_OBSERVED$string.Optimized_SpeciesNetwork

  #########################################################
  # Define subdirectory used for optimizing species tree2 #
  #########################################################
  string.Path_Directory_KH2_SpeciesNetwork2 = paste(string.Path_Directory_KH2N, '/Optimized_SpeciesNetwork_2', sep = "")
  unlink(string.Path_Directory_KH2_SpeciesNetwork2, recursive = T)
  dir.create(string.Path_Directory_KH2_SpeciesNetwork2, showWarnings = T, recursive = T)

  ######################################################################################
  # Compute observed test statistic for the differences in LnLs for the two topologies #
  ######################################################################################
  setwd(dir = string.Path_Directory_KH2_SpeciesNetwork2)
  handle.Optimized_SpeciesNetwork2_OBSERVED <- Optimize.Network(string.SpeciesNetwork = string.SpeciesNetwork2,
                                                                   handle.GeneTrees = handle.InputGeneTrees,
                                                                   string.PathDir = string.Path_Directory_KH2_SpeciesNetwork2)

  numeric.LnL_SpeciesNetwork2_OBSERVED <- handle.Optimized_SpeciesNetwork2_OBSERVED$numeric.MaximizedLnL
  string.Optimized_Network_2 <- handle.Optimized_SpeciesNetwork2_OBSERVED$string.Optimized_SpeciesNetwork

  ##################################
  # Step 1: Compute observed delta #
  ##################################
  numeric.Delta_Observed <- numeric.LnL_SpeciesNetwork1_OBSERVED - numeric.LnL_SpeciesNetwork2_OBSERVED

  ###########################################################################
  # Step 2.1: Compute gene tree likelihoods for each optimized species tree #
  ###########################################################################
  vector.GeneTree_Likelihoods_SpeciesNetwork1 <- Compute.GeneTree_Likelihoods_SpeciesNetwork(string.SpeciesNetwork = string.Optimized_Network_1,
                                                                           handle.GeneTrees = handle.InputGeneTrees,
                                                                           string.PathDir = string.Path_Directory_KH2_GeneTreeLikes_SpeciesNetwork1)

  vector.GeneTree_Likelihoods_SpeciesNetwork2 <- Compute.GeneTree_Likelihoods_SpeciesNetwork(string.SpeciesNetwork = string.Optimized_Network_2,
                                                                                             handle.GeneTrees = handle.InputGeneTrees,
                                                                                             string.PathDir = string.Path_Directory_KH2_GeneTreeLikes_SpeciesNetwork2)

  ##########################
  # Compute observed delta #
  ##########################
  vector.ObservedDelta <- vector.GeneTree_Likelihoods_SpeciesNetwork1 - vector.GeneTree_Likelihoods_SpeciesNetwork2

  ###############################################################################
  # Step 2.2: Conducting RELL bootstrap resampling of the gene tree likelihoods #
  ###############################################################################
  vector.Delta_BootstrapReplicates_RELL <- rep(NA, numeric.NumberOfReps)

  ###########################
  # Loop through replicates #
  ###########################
  for (i in 1:numeric.NumberOfReps){

    ###############################
    # Extract bootstrap replicate #
    ###############################
    print(gsub("Conducting BS replicate XXX...", pattern = "XXX", replacement = i))

    vector.DeltaBoostrapReplication_i <- sample(x = vector.ObservedDelta, size = numeric.NumberOfGeneTrees, replace = T)
    vector.Delta_BootstrapReplicates_RELL[i] <- sum(vector.DeltaBoostrapReplication_i)

  }

  ###############################################
  # Step 3: Center the delta LnLs by their mean #
  ###############################################
  vector.Delta_BootstrapReplicates <- vector.Delta_BootstrapReplicates_RELL[!is.na(vector.Delta_BootstrapReplicates_RELL)]
  vector.Centered_Delta_BootstrapReplicates <- vector.Delta_BootstrapReplicates - mean(vector.Delta_BootstrapReplicates)
  ###############################
  # Step 5: Compute signficance #
  ###############################
  numeric.Upper_Pvalue <- length(vector.Centered_Delta_BootstrapReplicates[vector.Centered_Delta_BootstrapReplicates>=numeric.Delta_Observed])/length(vector.Centered_Delta_BootstrapReplicates)
  numeric.Lower_Pvalue <- length(vector.Centered_Delta_BootstrapReplicates[vector.Centered_Delta_BootstrapReplicates<=numeric.Delta_Observed])/length(vector.Centered_Delta_BootstrapReplicates)

  numeric.MinPvalue <- min(c(numeric.Upper_Pvalue, numeric.Lower_Pvalue))
  numeric.TwoSided_Pvalue <- numeric.MinPvalue*2

  ##########################################
  # Return to directory and return results #
  ##########################################
  setwd(dir = string.CurrentDir)
  return(list(TwoSided_Pvalue = numeric.TwoSided_Pvalue,
              Upper_Pvalue = numeric.Upper_Pvalue,
              Lower_Pvalue = numeric.Lower_Pvalue,
              Observed_Delta = numeric.Delta_Observed,
              BS_Delta = vector.Delta_BootstrapReplicates))


}
