#' Conduct.KH_1N: function to conduct the KH_1N STAR test given two distinct species topologies (includes one or more network topologies) and a set of input gene trees
#'
#' This function returns a list containing p-values of the KH_1 STAR test for two input species tree topologies
#' @param string.SpeciesNetwork1 Strint defining the first species topology (can be network or bifurcating)
#' @param string.SpeciesNetwork2 Strint defining the first species topology (can be network or bifurcating)
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
#'                                                                     numeric.NumberOfGeneTrees = 10)
#' 
#' 
#' handle.KH_1N_Results <- Conduct.KH_1N(string.SpeciesNetwork1 = string.SpeciesNetwork, 
#'                                       string.SpeciesNetwork2 = string.SpeciesNetwork_2, 
#'                                       handle.InputGeneTrees = handle.SimulatedGeneTrees, 
#'                                       numeric.NumberOfReps = 3, 
#'                                       string.PathDir = '~/Desktop/')
#' 
#' 

################
# Conduct.KH_1N #
################
Conduct.KH_1N <- function(string.SpeciesNetwork1, string.SpeciesNetwork2, handle.InputGeneTrees, numeric.NumberOfReps, string.PathDir){
  
  ####################################################
  # Define directory used for conduct KH_1 STAR test #
  ####################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_KH1 = paste(string.PathDir, '/Conduct.KH_1N_', Sys.Date(), sep = "")
  unlink(string.Path_Directory_KH1, recursive = T)
  dir.create(string.Path_Directory_KH1, showWarnings = T, recursive = T)
  
  #############################################################
  # Define subdirectory used for bootstrap replicate analyses #
  #############################################################
  string.Path_Directory_BootStrapDir = paste(string.Path_Directory_KH1, '/BootStrapping', sep = "")
  unlink(string.Path_Directory_BootStrapDir, recursive = T)
  dir.create(string.Path_Directory_BootStrapDir, showWarnings = T, recursive = T)
  
  #########################################################
  # Define subdirectory used for optimizing species tree1 #
  #########################################################
  string.Path_Directory_KH1_SpeciesNetwork1 = paste(string.Path_Directory_KH1, '/Optimized_SpeciesNetwork1', sep = "")
  unlink(string.Path_Directory_KH1_SpeciesNetwork1, recursive = T)
  dir.create(string.Path_Directory_KH1_SpeciesNetwork1, showWarnings = T, recursive = T)
  
  ######################################################################################
  # Compute observed test statistic for the differences in LnLs for the two topologies #
  ######################################################################################
  setwd(dir = string.Path_Directory_KH1_SpeciesNetwork1)
  handle.Optimized_SpeciesNetwork1_OBSERVED <- Optimize.Network(string.SpeciesNetwork = string.SpeciesNetwork1,
                                                                   handle.GeneTrees = handle.InputGeneTrees,
                                                                   string.PathDir = string.Path_Directory_KH1_SpeciesNetwork1)
  numeric.LnL_SpeciesNetwork1_OBSERVED <- handle.Optimized_SpeciesNetwork1_OBSERVED$numeric.MaximizedLnL
  
  #########################################################
  # Define subdirectory used for optimizing species tree2 #
  #########################################################
  string.Path_Directory_KH1_SpeciesNetwork2 = paste(string.Path_Directory_KH1, '/Optimized_SpeciesNetwork2', sep = "")
  unlink(string.Path_Directory_KH1_SpeciesNetwork2, recursive = T)
  dir.create(string.Path_Directory_KH1_SpeciesNetwork2, showWarnings = T, recursive = T)
  
  ######################################################################################
  # Compute observed test statistic for the differences in LnLs for the two topologies #
  ######################################################################################
  setwd(dir = string.Path_Directory_KH1_SpeciesNetwork2)
  handle.Optimized_SpeciesNetwork2_OBSERVED <- Optimize.Network(string.SpeciesNetwork = string.SpeciesNetwork2,
                                                                handle.GeneTrees = handle.InputGeneTrees,
                                                                string.PathDir = string.Path_Directory_KH1_SpeciesNetwork2)
  
  numeric.LnL_SpeciesNetwork2_OBSERVED <- handle.Optimized_SpeciesNetwork2_OBSERVED$numeric.MaximizedLnL
  
  ##################################
  # Step 1: Compute observed delta #
  ##################################
  numeric.Delta_Observed <- numeric.LnL_SpeciesNetwork1_OBSERVED - numeric.LnL_SpeciesNetwork2_OBSERVED
  
  ########################################
  # Step 2: Conduct BootStrap Resampling #
  ########################################
  list.BoostrapReplicates_GeneTrees <- Resample.Bootstrap_GeneTree_Replicates_NP(handle.GeneTrees = handle.InputGeneTrees,
                                                                                 numeric.NumberOfReps = numeric.NumberOfReps)
  
  ###############################################################
  # Step 3: Optimize species trees for each bootstrap replicate #
  ###############################################################
  vector.Delta_BootstrapReplicates <- rep(NA, numeric.NumberOfReps)
  
  ###########################
  # Loop through replicates #
  ###########################
  for (i in 1:numeric.NumberOfReps){
    
    ###############################
    # Extract bootstrap replicate #
    ###############################
    print(gsub("Conducting BS replicate XXX...", pattern = "XXX", replacement = i))
    handle.BootStrapped_GeneTrees_i <- list.BoostrapReplicates_GeneTrees[[i]]
    
    #############################################################
    # Define subdirectory used for bootstrap replicate analyses #
    #############################################################
    string.Path_Directory_BootStrapRep_i = paste(string.Path_Directory_BootStrapDir, '/Rep_', i, sep = "")
    unlink(string.Path_Directory_BootStrapRep_i, recursive = T)
    dir.create(string.Path_Directory_BootStrapRep_i, showWarnings = T, recursive = T)
    
    ##################################################
    # Define subdirectory used for optimizing tree 1 #
    ##################################################
    string.Path_Directory_BootStrapRep_SpeciesNetwork1_i = paste(string.Path_Directory_BootStrapRep_i, '/Optimized_SpeciesNetwork1_', i, sep = "")
    unlink(string.Path_Directory_BootStrapRep_SpeciesNetwork1_i, recursive = T)
    dir.create(string.Path_Directory_BootStrapRep_SpeciesNetwork1_i, showWarnings = T, recursive = T)
    
    ############################################
    # Optimize species tree 1 for BS replicate #
    ############################################
    handle.Optimized_SpeciesNetwork_1_i <- Optimize.Network(string.SpeciesNetwork = string.SpeciesNetwork1,
                                                               handle.GeneTrees = handle.BootStrapped_GeneTrees_i,
                                                               string.PathDir = string.Path_Directory_BootStrapRep_SpeciesNetwork1_i)
    numeric.LnL_SpeciesNetwork1_i <- handle.Optimized_SpeciesNetwork_1_i$numeric.MaximizedLnL
    
    ##################################################
    # Define subdirectory used for optimizing tree 2 #
    ##################################################
    string.Path_Directory_BootStrapRep_SpeciesNetwork2_i = paste(string.Path_Directory_BootStrapRep_i, '/Optimized_SpeciesNetwork2_', i, sep = "")
    unlink(string.Path_Directory_BootStrapRep_SpeciesNetwork2_i, recursive = T)
    dir.create(string.Path_Directory_BootStrapRep_SpeciesNetwork2_i, showWarnings = T, recursive = T)
    
    ############################################
    # Optimize species tree 2 for BS replicate #
    ############################################
    handle.Optimized_SpeciesNetwork_2_i <- Optimize.Network(string.SpeciesNetwork = string.SpeciesNetwork2,
                                                               handle.GeneTrees = handle.BootStrapped_GeneTrees_i,
                                                               string.PathDir = string.Path_Directory_BootStrapRep_SpeciesNetwork2_i)
    numeric.LnL_SpeciesNetwork2_i <- handle.Optimized_SpeciesNetwork_2_i$numeric.MaximizedLnL
    
    ################################################
    # Compute test statistic delta for replicate i #
    ################################################
    numeric.Delta_LnL_BootStrapReplicate_i <- numeric.LnL_SpeciesNetwork1_i - numeric.LnL_SpeciesNetwork2_i
    vector.Delta_BootstrapReplicates[i] <- numeric.Delta_LnL_BootStrapReplicate_i
    
    ###############################
    # Compress subdirectory files #
    ###############################
    setwd(dir = string.Path_Directory_BootStrapRep_i)
    system(command = "find . -type f -exec bzip2 -9 {} +")
    
  }
  
  ###############################################
  # Step 4: Center the delta LnLs by their mean #
  ###############################################
  vector.Delta_BootstrapReplicates <- vector.Delta_BootstrapReplicates[!is.na(vector.Delta_BootstrapReplicates)]
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