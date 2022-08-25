#' Conduct.KH_3N: function to conduct the KH_2N STAR test given two distinct species topologies (one or both are networks) and a set of input gene trees
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

#################
# Conduct.KH_3N #
#################
Conduct.KH_3N <- function(string.SpeciesNetwork1, string.SpeciesNetwork2, handle.InputGeneTrees, numeric.NumberOfReps, string.PathDir){


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
    vector.BootstrapDelta <- sample(x = vector.ObservedDelta, size = length(vector.ObservedDelta), replace = T)
    vector.Delta_BootstrapReplicates_RELL[i] <- sum(vector.BootstrapDelta)

  }

  ###############################################
  # Step 3: Center the delta LnLs by their mean #
  ###############################################
  vector.Delta_BootstrapReplicates <- vector.Delta_BootstrapReplicates_RELL[!is.na(vector.Delta_BootstrapReplicates_RELL)]
  vector.Centered_Delta_BootstrapReplicates <- vector.Delta_BootstrapReplicates - mean(vector.Delta_BootstrapReplicates)

  ############################
  # Step 5: Compute variance #
  ############################
  numeric.Variance_BootstrapReplicates <- var(vector.Centered_Delta_BootstrapReplicates)
  numeric.Pvalue <- 2*(1-pnorm(q = numeric.Delta_Observed, mean = 0, sd = sqrt(numeric.Variance_BootstrapReplicates)))
  pvalue2sided=2*pnorm(-abs(numeric.Delta_Observed),mean = 0, sd = sqrt(numeric.Variance_BootstrapReplicates))

  print(c(numeric.Delta_Observed, "Observed delta"))
  #print(numeric.Variance_BootstrapReplicates)
  print(c(sqrt(numeric.Variance_BootstrapReplicates), "estimated variance"))
  print(c(pvalue2sided, "two-sided P-value"))
}
