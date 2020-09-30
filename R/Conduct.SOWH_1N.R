#' Conduct.SOWH_1N: function to conduct the SOWH_1N STAR test given a species network topology
#'
#' This function returns a list containing p-values of the SOWH_1 STAR test
#' @param string.SpeciesNetwork_1 String defining the first species topology (can be network or bifurcating)
#' @param handle.GeneTrees Phylo object containing a list of the input gene trees
#' @param numeric.NumberOfReps Number of bootstrap replicates to analyze
#' @param numeric.MaxReticulations Numeric defining the maximum number of reticulations for PhyloNet search
#' @param string.PathDir String defining the path to a parent directory used for conduct SOWH_1 STAR test
#' @keywords Species tree, multispecies coalescent, SWOH, Network
#' @return List Returns a list containing (1) vector.Null_BS_Delta_Stat, (2) numeric.pValue, and (3) numeric.Delta_Observed
#' @export
#' @examples
#' #'
#' #'
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
#'                                                                        numeric.NumberOfGeneTrees = 3)
#' 
#' 
#' 
#' ######################################
#' # Conduct SOWH_1N STAR using network #
#' ######################################
#' Conduct.SOWH_1N(string.SpeciesNetwork_1 = string.SpeciesNetwork_2, 
#'                 handle.InputGeneTrees = handle.Simulated_GeneTreeSet, 
#'                 numeric.NumberOfReps = 2, 
#'                 numeric.MaxReticulations = 1,
#'                 string.PathDir = '~/Desktop/')

###################
# Conduct.SOWH_1N #
###################
Conduct.SOWH_1N <- function(string.SpeciesNetwork_1, handle.InputGeneTrees, numeric.NumberOfReps, numeric.MaxReticulations, string.PathDir){
  
  ###################
  # Summarize input #
  ###################
  numeric.NumberOfGeneTrees <- length(handle.InputGeneTrees)
  
  ####################################################
  # Define directory used for conduct SOWH_1 STAR test #
  ####################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_SOWH_1N = paste(string.PathDir, '/Conduct.SOWH_1N_', Sys.Date(), sep = "")
  unlink(string.Path_Directory_SOWH_1N, recursive = T)
  dir.create(string.Path_Directory_SOWH_1N, showWarnings = T, recursive = T)
  
  #############################################################
  # Define subdirectory used for bootstrap replicate analyses #
  #############################################################
  string.Path_Directory_BootStrapDir = paste(string.Path_Directory_SOWH_1N, '/Parametric_BootStrapping', sep = "")
  unlink(string.Path_Directory_BootStrapDir, recursive = T)
  dir.create(string.Path_Directory_BootStrapDir, showWarnings = T, recursive = T)
  
  ############################################################
  # Define subdirectory used for optimizing species network1 #
  ############################################################
  string.Path_Directory_SOWH_SpeciesNetwork1 = paste(string.Path_Directory_SOWH_1N, '/Optimized_SpeciesNetwork1', sep = "")
  unlink(string.Path_Directory_SOWH_SpeciesNetwork1, recursive = T)
  dir.create(string.Path_Directory_SOWH_SpeciesNetwork1, showWarnings = T, recursive = T)
  
  ######################################################################################
  # Compute observed test statistic for the differences in LnLs for the two topologies #
  ######################################################################################
  setwd(dir = string.Path_Directory_SOWH_SpeciesNetwork1)
  handle.Optimized_SpeciesNetwork1_OBSERVED <- Optimize.Network(string.SpeciesNetwork = string.SpeciesNetwork_1,
                                                                handle.GeneTrees = handle.InputGeneTrees,
                                                                string.PathDir = string.Path_Directory_SOWH_SpeciesNetwork1)
  numeric.LnL_SpeciesNetwork1_OBSERVED <- handle.Optimized_SpeciesNetwork1_OBSERVED$numeric.MaximizedLnL
  string.Optimized_Network_1 <- handle.Optimized_SpeciesNetwork1_OBSERVED$string.Optimized_SpeciesNetwork
  
  ##############################################################
  # Define subdirectory used for optimizing species network ML #
  ##############################################################
  string.Path_Directory_SOWH_Find_ML_N = paste(string.Path_Directory_SOWH_1N, '/Optimize_Find_ML_Network', sep = "")
  unlink(string.Path_Directory_SOWH_Find_ML_N, recursive = T)
  dir.create(string.Path_Directory_SOWH_Find_ML_N, showWarnings = T, recursive = T)
  
  ######################################################################################
  # Compute observed test statistic for the differences in LnLs for the two topologies #
  ######################################################################################
  setwd(dir = string.Path_Directory_SOWH_Find_ML_N)
  handle.Optimized_SpeciesNetwork_OBSERVED_ML <- Optimize.NetworkSearch(handle.GeneTrees = handle.InputGeneTrees,
                                                                        numeric.MaxReticulations = numeric.MaxReticulations,
                                                                        string.PathDir = string.Path_Directory_SOWH_Find_ML_N)
  numeric.LnL_SpeciesNetwork_OBSERVED_ML <- handle.Optimized_SpeciesNetwork_OBSERVED_ML$numeric.MaximizedLnL
  string.Optimized_SpeciesNetwork_ML <- handle.Optimized_SpeciesNetwork_OBSERVED_ML$string.Optimized_SpeciesNetwork
  
  ##################################
  # Step 1: Compute observed delta #
  ##################################
  numeric.Delta_Observed <- numeric.LnL_SpeciesNetwork_OBSERVED_ML - numeric.LnL_SpeciesNetwork1_OBSERVED
  
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
    handle.Parametric_Simulated_GeneTrees <- Simulate.GeneTrees_From_SpeciesNetwork(string.SpeciesNetwork = string.Optimized_Network_1,
                                                                                    numeric.NumberOfGeneTrees =numeric.NumberOfGeneTrees,
                                                                                    string.PathDir = string.Path_Directory_BootStrapRep_i)
    ####################
    # Find ML topology #
    ####################
    handle.Optimized_SpeciesNetwork_OBSERVED_ML_i <- Optimize.NetworkSearch(handle.GeneTrees = handle.Parametric_Simulated_GeneTrees,
                                                                            numeric.MaxReticulations = numeric.MaxReticulations,
                                                                            string.PathDir = string.Path_Directory_BootStrapRep_i)
    numeric.LnL_SpeciesNetwork_OBSERVED_ML_i <- handle.Optimized_SpeciesNetwork_OBSERVED_ML_i$numeric.MaximizedLnL
    
    ###############
    # Optimize T1 #
    ###############
    handle.Optimized_SpeciesNetwork1_OBSERVED_i <- Optimize.Network(string.SpeciesNetwork = string.Optimized_Network_1,
                                                                       handle.GeneTrees = handle.Parametric_Simulated_GeneTrees,
                                                                       string.PathDir = string.Path_Directory_BootStrapRep_i)
    numeric.LnL_SpeciesNetwork1_OBSERVED_i <- handle.Optimized_SpeciesNetwork1_OBSERVED_i$numeric.MaximizedLnL
    
    #####################################
    # Compute test statistics for rep i #
    #####################################
    numeric.Delta_Null_i <- numeric.LnL_SpeciesNetwork_OBSERVED_ML_i - numeric.LnL_SpeciesNetwork1_OBSERVED_i
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