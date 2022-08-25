#' Conduct.KH_1: function to conduct the KH_1 STAR test given two distinct species topologies and a set of input gene trees
#'
#' This function returns a list containing p-values of the KH_1 STAR test for two input species tree topologies
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
#' handle.SpeciesTree1$edge.length <- handle.SpeciesTree1$edge.length*0.2
#' handle.SpeciesTree2 <- read.tree(text = "(C:1.0,(B:0.5,(A:1.0,D:1.0):1.5):0.5);")
#'
#' #################################################
#' # Simlate a set of gene trees for species tree1 #
#' #################################################
#' handle.Simulated_GeneTrees <-  system.file("extdata", "ExampleGeneTreeSet.tree", package="SpeciesTopoTestR") # load previous simulations
#' handle.Simulated_GeneTrees <- read.tree(file = handle.Simulated_GeneTrees) # load previous simulations
#' 
#' handle.Simulated_GeneTrees <- Simulate.GeneTrees_From_SpeciesTree3(handle.SpeciesTree = handle.SpeciesTree1,
#'                                                                   string.PathDir = '~/Desktop/',
#'                                                                   numeric.NumberOfGeneTrees = 10)
#'
#' ##############################################################
#' # Conduct KH_1 STAR test for the two species tree topologies #
#' ##############################################################
#' Conduct.KH_1(handle.SpeciesTree1 = handle.SpeciesTree2,
#'              handle.SpeciesTree2 = handle.SpeciesTree1,
#'              handle.InputGeneTrees = handle.Simulated_GeneTrees,
#'              numeric.NumberOfReps = 100,
#'              string.PathDir = '~/Desktop/')

################
# Conduct.KH_1 #
################
Conduct.KH_1 <- function(handle.SpeciesTree1, handle.SpeciesTree2, handle.InputGeneTrees, numeric.NumberOfReps, string.PathDir){


  ####################################################
  # Define directory used for conduct KH_1 STAR test #
  ####################################################
  string.CurrentDir <- getwd()
  string.Path_Directory_KH1 = paste(string.PathDir, '/Conduct.KH_1_', Sys.Date(), sep = "")
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
  string.Path_Directory_KH1_SpeciesTree1 = paste(string.Path_Directory_KH1, '/Optimized_SpeciesTree1', sep = "")
  unlink(string.Path_Directory_KH1_SpeciesTree1, recursive = T)
  dir.create(string.Path_Directory_KH1_SpeciesTree1, showWarnings = T, recursive = T)

  ######################################################################################
  # Compute observed test statistic for the differences in LnLs for the two topologies #
  ######################################################################################
  setwd(dir = string.Path_Directory_KH1_SpeciesTree1)
  handle.Optimized_SpeciesTree1_OBSERVED <- Optimize.BranchLengths(handle.SpeciesTree = handle.SpeciesTree1,
                                                                    handle.GeneTrees = handle.InputGeneTrees,
                                                                    string.PathDir = string.Path_Directory_KH1_SpeciesTree1,
                                                                    numeric.Stells_algorthm = 1)
  numeric.LnL_SpeciesTree1_OBSERVED <- handle.Optimized_SpeciesTree1_OBSERVED$numeric.MaximizedLnL

  #########################################################
  # Define subdirectory used for optimizing species tree2 #
  #########################################################
  string.Path_Directory_KH1_SpeciesTree2 = paste(string.Path_Directory_KH1, '/Optimized_SpeciesTree2', sep = "")
  unlink(string.Path_Directory_KH1_SpeciesTree2, recursive = T)
  dir.create(string.Path_Directory_KH1_SpeciesTree2, showWarnings = T, recursive = T)

  ######################################################################################
  # Compute observed test statistic for the differences in LnLs for the two topologies #
  ######################################################################################
  setwd(dir = string.Path_Directory_KH1_SpeciesTree2)
  handle.Optimized_SpeciesTree2_OBSERVED <- Optimize.BranchLengths(handle.SpeciesTree = handle.SpeciesTree2,
                                                                   handle.GeneTrees = handle.InputGeneTrees,
                                                                   string.PathDir = string.Path_Directory_KH1_SpeciesTree2,
                                                                   numeric.Stells_algorthm = 1)
  numeric.LnL_SpeciesTree2_OBSERVED <- handle.Optimized_SpeciesTree2_OBSERVED$numeric.MaximizedLnL

  ##################################
  # Step 1: Compute observed delta #
  ##################################
  numeric.Delta_Observed <- numeric.LnL_SpeciesTree1_OBSERVED - numeric.LnL_SpeciesTree2_OBSERVED

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
    string.Path_Directory_BootStrapRep_SpeciesTree1_i = paste(string.Path_Directory_BootStrapRep_i, '/Optimized_SpeciesTree1_', i, sep = "")
    unlink(string.Path_Directory_BootStrapRep_SpeciesTree1_i, recursive = T)
    dir.create(string.Path_Directory_BootStrapRep_SpeciesTree1_i, showWarnings = T, recursive = T)

    ############################################
    # Optimize species tree 1 for BS replicate #
    ############################################
    handle.Optimized_SpeciesTree_1_i <- Optimize.BranchLengths(handle.SpeciesTree = handle.SpeciesTree1,
                                                               handle.GeneTrees = handle.BootStrapped_GeneTrees_i,
                                                               string.PathDir = string.Path_Directory_BootStrapRep_SpeciesTree1_i,
                                                               numeric.Stells_algorthm = 1)
    numeric.LnL_SpeciesTree1_i <- handle.Optimized_SpeciesTree_1_i$numeric.MaximizedLnL

    ##################################################
    # Define subdirectory used for optimizing tree 2 #
    ##################################################
    string.Path_Directory_BootStrapRep_SpeciesTree2_i = paste(string.Path_Directory_BootStrapRep_i, '/Optimized_SpeciesTree2_', i, sep = "")
    unlink(string.Path_Directory_BootStrapRep_SpeciesTree2_i, recursive = T)
    dir.create(string.Path_Directory_BootStrapRep_SpeciesTree2_i, showWarnings = T, recursive = T)

    ############################################
    # Optimize species tree 1 for BS replicate #
    ############################################
    handle.Optimized_SpeciesTree_2_i <- Optimize.BranchLengths(handle.SpeciesTree = handle.SpeciesTree2,
                                                               handle.GeneTrees = handle.BootStrapped_GeneTrees_i,
                                                               string.PathDir = string.Path_Directory_BootStrapRep_SpeciesTree2_i,
                                                               numeric.Stells_algorthm = 1)
    numeric.LnL_SpeciesTree2_i <- handle.Optimized_SpeciesTree_2_i$numeric.MaximizedLnL


    ################################################
    # Compute test statistic delta for replicate i #
    ################################################
    numeric.Delta_LnL_BootStrapReplicate_i <- numeric.LnL_SpeciesTree1_i - numeric.LnL_SpeciesTree2_i
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
