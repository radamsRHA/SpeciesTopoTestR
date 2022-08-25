#' Conduct.SOWH_2v2: function to conduct the SOWH_1 STAR test given two distinct species topologies and a set of input gene trees
#'
#' This function returns a list containing p-values of the SOWH_1 STAR test
#' @param handle.SpeciesTree_1 Phylogenetic tree defining the first species topology
#' @param handle.GeneTrees Phylo object containing a list of the input gene trees
#' @param numeric.NumberOfReps Number of bootstrap replicates to analyze
#' @param string.PathDir String defining the path to a parent directory used for conduct SOWH_1 STAR test
#' @keywords Species tree, multispecies coalescent, phylogenetics, phylogenomics
#' @return List Returns a list containing (1) vector.Null_BS_Delta_Stat, (2) numeric.pValue, and (3) numeric.Delta_Observed
#' @export
#' @examples
#' #'
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
#' handle.SpeciesTreeX1 <- read.tree(text = "(A:1.0,(B:0.5,(C:1.0,D:1.0):1.5):0.5);")
#' handle.SpeciesTreeX1$edge.length <- handle.SpeciesTreeX1$edge.length * 1
#' handle.SpeciesTreeX2 <- read.tree(text = "(C:1.0,(B:0.5,(A:1.0,D:1.0):1.5):0.5);")
#'
#'
#'
#' #################################################
#' # Simlate a set of gene trees for species tree1 #
#' #################################################
#' handle.Simulated_GeneTrees <- Simulate.GeneTrees_From_SpeciesTree(handle.SpeciesTree = handle.SpeciesTreeX1,
#'                                                                   string.PathDir = '~/Desktop/',
#'                                                                   numeric.NumberOfGeneTrees = 10)
#'
#'
#'################################
#'# Conduct SOWH_2 STAR test for S #
#'################################
#'Conduct.SOWH_2v2(handle.SpeciesTree_1 = handle.SpeciesTreeX2,
#'               handle.InputGeneTrees = handle.Simulated_GeneTrees,
#'               numeric.NumberOfReps = 100,
#'               string.PathDir = '~/Desktop/')
#'
#'
#'

##################
# Conduct.SOWH_2v2 #
#################
Conduct.SOWH_2v2 <- function(handle.SpeciesTree_1, handle.InputGeneTrees, numeric.NumberOfReps, string.PathDir){
  
  ###################
  # Summarize input #
  ###################
  numeric.NumberOfGeneTrees <- length(handle.InputGeneTrees)
  
  ####################################################
  # Define directory used for conduct SOWH_1 STAR test #
  ####################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_SOWH_2 = paste(string.PathDir, '/Conduct.SOWH_2_', Sys.Date(), sep = "")
  unlink(string.Path_Directory_SOWH_2, recursive = T)
  dir.create(string.Path_Directory_SOWH_2, showWarnings = T, recursive = T)
  
  #############################################################
  # Define subdirectory used for bootstrap replicate analyses #
  #############################################################
  string.Path_Directory_BootStrapDir = paste(string.Path_Directory_SOWH_2, '/Parametric_BootStrapping', sep = "")
  unlink(string.Path_Directory_BootStrapDir, recursive = T)
  dir.create(string.Path_Directory_BootStrapDir, showWarnings = T, recursive = T)
  
  #########################################################
  # Define subdirectory used for optimizing species tree1 #
  #########################################################
  string.Path_Directory_SOWH_SpeciesTree1 = paste(string.Path_Directory_SOWH_2, '/Optimized_SpeciesTree1', sep = "")
  unlink(string.Path_Directory_SOWH_SpeciesTree1, recursive = T)
  dir.create(string.Path_Directory_SOWH_SpeciesTree1, showWarnings = T, recursive = T)
  
  ######################################################################################
  # Compute observed test statistic for the differences in LnLs for the two topologies #
  ######################################################################################
  setwd(dir = string.Path_Directory_SOWH_SpeciesTree1)
  handle.Optimized_SpeciesTree1_OBSERVED <- Optimize.BranchLengths(handle.SpeciesTree = handle.SpeciesTree_1,
                                                                   handle.GeneTrees = handle.InputGeneTrees,
                                                                   string.PathDir = string.Path_Directory_SOWH_SpeciesTree1,
                                                                   numeric.Stells_algorthm = 1)
  numeric.LnL_SpeciesTree1_OBSERVED <- handle.Optimized_SpeciesTree1_OBSERVED$numeric.MaximizedLnL
  handle.Optimized_SpeciesTopology_1 <- handle.Optimized_SpeciesTree1_OBSERVED$handle.Optimized_SpeciesTree
  
  ###########################################################
  # Define subdirectory used for optimizing species tree ML #
  ###########################################################
  string.Path_Directory_SOWH_Find_ML = paste(string.Path_Directory_SOWH_2, '/Optimize_Find_ML', sep = "")
  unlink(string.Path_Directory_SOWH_Find_ML, recursive = T)
  dir.create(string.Path_Directory_SOWH_Find_ML, showWarnings = T, recursive = T)
  
  ######################################################################################
  # Compute observed test statistic for the differences in LnLs for the two topologies #
  ######################################################################################
  setwd(dir = string.Path_Directory_SOWH_Find_ML)
  handle.Optimized_SpeciesTree_OBSERVED_ML <- Optimize.TopologySearch(handle.GeneTrees = handle.InputGeneTrees,
                                                                      string.PathDir = string.Path_Directory_SOWH_Find_ML)
  numeric.LnL_SpeciesTree_OBSERVED_ML <- handle.Optimized_SpeciesTree_OBSERVED_ML$numeric.MaximizedLnL
  handle.Optimized_SpeciesTopology_ML <- handle.Optimized_SpeciesTree_OBSERVED_ML$handle.Optimized_SpeciesTree
  
  ##################################
  # Step 1: Compute observed delta #
  ##################################
  numeric.Delta_Observed <- numeric.LnL_SpeciesTree_OBSERVED_ML - numeric.LnL_SpeciesTree1_OBSERVED
  
  ############################################
  # Step 2: Parametric bootstrap simulations #
  ############################################
  vector.Null_BS_Delta_Stat <- rep(NA, numeric.NumberOfReps)
  
  for (i in 1:numeric.NumberOfReps){
    
    ###############################
    # Extract bootstrap replicate #
    ###############################
    print(gsub("Conducting BS replicate XXX...", pattern = "XXX", replacement = i))
    
    #############################################################
    # Define subdirectory used for bootstrap replicate analyses #
    #############################################################
    string.Path_Directory_BootStrapRep_i = paste(string.Path_Directory_BootStrapDir, '/Rep_', i, sep = "")
    unlink(string.Path_Directory_BootStrapRep_i, recursive = T)
    dir.create(string.Path_Directory_BootStrapRep_i, showWarnings = T, recursive = T)
    
    ######################################
    # Parametric BS simulations using T1 #
    ######################################
    setwd(dir = string.Path_Directory_BootStrapRep_i)
    handle.Parametric_Simulated_GeneTrees <- Simulate.GeneTrees_From_SpeciesTree2(handle.SpeciesTree = handle.Optimized_SpeciesTopology_1,
                                                                                 numeric.NumberOfGeneTrees = numeric.NumberOfGeneTrees,
                                                                                 string.PathDir = string.Path_Directory_BootStrapRep_i)
    
    ####################
    # Find ML topology #
    ####################
    handle.Optimized_SpeciesTree_OBSERVED_ML_i <- Optimize.TopologySearch(handle.GeneTrees = handle.Parametric_Simulated_GeneTrees,
                                                                          string.PathDir = string.Path_Directory_BootStrapRep_i)
    numeric.LnL_SpeciesTree_OBSERVED_ML_i <- handle.Optimized_SpeciesTree_OBSERVED_ML_i$numeric.MaximizedLnL
    
    ##########################################
    # Compute gene tree likelihoods given T1 #
    ##########################################
    vector.GeneTreeLikelihoods_T1_i <- Compute.GeneTree_Likelihoods(handle.SpeciesTree = handle.Optimized_SpeciesTopology_1,
                                                                    handle.GeneTrees = handle.Parametric_Simulated_GeneTrees,
                                                                    string.PathDir = string.Path_Directory_BootStrapRep_i)
    numeric.LnL_SpeciesTree_OBSERVED_T1_i <- sum(vector.GeneTreeLikelihoods_T1_i)
    
    #####################################
    # Compute test statistics for rep i #
    #####################################
    numeric.Delta_Null_i <- numeric.LnL_SpeciesTree_OBSERVED_ML_i - numeric.LnL_SpeciesTree_OBSERVED_T1_i
    vector.Null_BS_Delta_Stat[i] <- numeric.Delta_Null_i
  }
  
  ############################
  # Step 5: Compute p-values #
  ############################
  numeric.pValue <- length(vector.Null_BS_Delta_Stat[vector.Null_BS_Delta_Stat>=numeric.Delta_Observed])/numeric.NumberOfReps
  
  ##########################################
  # Return to directory and return results #
  ##########################################
  setwd(dir = string.CurrentDir)
  return(list(vector.Null_BS_Delta_Stat = vector.Null_BS_Delta_Stat,
              numeric.pValue = numeric.pValue,
              numeric.Delta_Observed = numeric.Delta_Observed))
  
}
