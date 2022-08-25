#' Conduct.KH_3: function to conduct the KH_2 STAR test given two distinct species topologies and a set of input gene trees
#'
#' This function returns a list containing p-values of the KH_2 STAR test for two input species tree topologies
#' @param handle.SpeciesTree1 Phylogenetic tree defining the first species topology
#' @param handle.SpeciesTree2 Phylogenetic tree defining the second species topology
#' @param handle.GeneTrees Phylo object containing a list of the input gene trees
#' @param numeric.NumberOfReps Number of bootstrap replicates to analyze
#' @param string.PathDir String defining the path to a parent directory used for conduct KH_1 STAR test
#' @keywords Species tree, multispecies coalescent, phylogenetics, phylogenomics
#' @return List Returns a list containing (1) twosided pvalue, (2) upper p-values, (3) lower p-values, and (4) a vector of the bootstrapped test statistics delta
#' @export
#' @examples
#'
#'
#'

################
# Conduct.KH_3 #
################
Conduct.KH_3 <- function(handle.SpeciesTree1, handle.SpeciesTree2, handle.InputGeneTrees, numeric.NumberOfReps, string.PathDir){

  ###################
  # Summarize input #
  ###################
  numeric.NumberOfGeneTrees <- length(handle.InputGeneTrees)

  ####################################################
  # Define directory used for conduct KH_1 STAR test #
  ####################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_KH2 = paste(string.PathDir, '/Conduct.KH_2_', Sys.Date(), sep = "")
  unlink(string.Path_Directory_KH2, recursive = T)
  dir.create(string.Path_Directory_KH2, showWarnings = T, recursive = T)

  #################################################################################
  # Define subdirectory used for comput likelihoods given optimized species tree1 #
  #################################################################################
  string.Path_Directory_KH2_GeneTreeLikes_SpeciesTree1 = paste(string.Path_Directory_KH2, '/GeneTreeLikes_Optimized_SpeciesTree1', sep = "")
  unlink(string.Path_Directory_KH2_GeneTreeLikes_SpeciesTree1, recursive = T)
  dir.create(string.Path_Directory_KH2_GeneTreeLikes_SpeciesTree1, showWarnings = T, recursive = T)

  #################################################################################
  # Define subdirectory used for comput likelihoods given optimized species tree2 #
  #################################################################################
  string.Path_Directory_KH2_GeneTreeLikes_SpeciesTree2 = paste(string.Path_Directory_KH2, '/GeneTreeLikes_Optimized_SpeciesTree2', sep = "")
  unlink(string.Path_Directory_KH2_GeneTreeLikes_SpeciesTree2, recursive = T)
  dir.create(string.Path_Directory_KH2_GeneTreeLikes_SpeciesTree2, showWarnings = T, recursive = T)

  #########################################################
  # Define subdirectory used for optimizing species tree1 #
  #########################################################
  string.Path_Directory_KH2_SpeciesTree1 = paste(string.Path_Directory_KH2, '/Optimized_SpeciesTree1', sep = "")
  unlink(string.Path_Directory_KH2_SpeciesTree1, recursive = T)
  dir.create(string.Path_Directory_KH2_SpeciesTree1, showWarnings = T, recursive = T)

  ######################################################################################
  # Compute observed test statistic for the differences in LnLs for the two topologies #
  ######################################################################################
  setwd(dir = string.Path_Directory_KH2_SpeciesTree1)
  handle.Optimized_SpeciesTree1_OBSERVED <- Optimize.BranchLengths(handle.SpeciesTree = handle.SpeciesTree1,
                                                                   handle.GeneTrees = handle.InputGeneTrees,
                                                                   string.PathDir = string.Path_Directory_KH2_SpeciesTree1,
                                                                   numeric.Stells_algorthm = 1)
  numeric.LnL_SpeciesTree1_OBSERVED <- handle.Optimized_SpeciesTree1_OBSERVED$numeric.MaximizedLnL


  #########################################################
  # Define subdirectory used for optimizing species tree2 #
  #########################################################
  string.Path_Directory_KH2_SpeciesTree2 = paste(string.Path_Directory_KH2, '/Optimized_SpeciesTree2', sep = "")
  unlink(string.Path_Directory_KH2_SpeciesTree2, recursive = T)
  dir.create(string.Path_Directory_KH2_SpeciesTree2, showWarnings = T, recursive = T)

  ######################################################################################
  # Compute observed test statistic for the differences in LnLs for the two topologies #
  ######################################################################################
  setwd(dir = string.Path_Directory_KH2_SpeciesTree1)
  handle.Optimized_SpeciesTree2_OBSERVED <- Optimize.BranchLengths(handle.SpeciesTree = handle.SpeciesTree2,
                                                                   handle.GeneTrees = handle.InputGeneTrees,
                                                                   string.PathDir = string.Path_Directory_KH2_SpeciesTree2,
                                                                   numeric.Stells_algorthm = 1)
  numeric.LnL_SpeciesTree2_OBSERVED <- handle.Optimized_SpeciesTree2_OBSERVED$numeric.MaximizedLnL

  ##################################
  # Step 1: Compute observed delta #
  ##################################
  numeric.Delta_Observed <- numeric.LnL_SpeciesTree1_OBSERVED - numeric.LnL_SpeciesTree2_OBSERVED

  ###########################################################################
  # Step 2.1: Compute gene tree likelihoods for each optimized species tree #
  ###########################################################################
  vector.GeneTree_Likelihoods_SpeciesTree1 <- Compute.GeneTree_Likelihoods(handle.SpeciesTree = handle.Optimized_SpeciesTree1_OBSERVED$handle.Optimized_SpeciesTree,
                                                                           handle.GeneTrees = handle.InputGeneTrees,
                                                                           string.PathDir = string.Path_Directory_KH2_GeneTreeLikes_SpeciesTree1)
  vector.GeneTree_Likelihoods_SpeciesTree2 <- Compute.GeneTree_Likelihoods(handle.SpeciesTree = handle.Optimized_SpeciesTree2_OBSERVED$handle.Optimized_SpeciesTree,
                                                                           handle.GeneTrees = handle.InputGeneTrees,
                                                                           string.PathDir = string.Path_Directory_KH2_GeneTreeLikes_SpeciesTree2)

  ##########################
  # Compute observed delta #
  ##########################
  vector.ObservedDelta <- vector.GeneTree_Likelihoods_SpeciesTree1 - vector.GeneTree_Likelihoods_SpeciesTree2

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
    vector.BootStrapDelta_i <- sample(x = vector.ObservedDelta, size = length(vector.ObservedDelta), replace = T)
    vector.Delta_BootstrapReplicates_RELL[i] <- sum(vector.BootStrapDelta_i)

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

  print(c(numeric.Delta_Observed, "Observed Delta"))
  #print(numeric.Variance_BootstrapReplicates)
  print(c(sqrt(numeric.Variance_BootstrapReplicates), "Estimated Variance"))
  print(c(pvalue2sided, "two-sided p-value"))

}