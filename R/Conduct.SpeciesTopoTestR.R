#' Conduct.SpeciesTopoTestR: function to conduct an array of different likelihood-based tests of species topologies
#'
#' This function returns a list containing the results of a topology test
#' @param handle.Topologies2Test List of topologies to be tested. Use format "multiPhylo" for bifurcating topologies, and a standard list for networks (each network is just a string)
#' @param handle.GeneTrees Phylo object containing a list of the input gene trees
#' @param string.Test String defining the test to run, can be "KH", "SH", "SOWH"
#' @param numeric.Algorithm Numeric specifiying the algorithm KH (1, 2 or 3), SH (1 or 2), SOWH (1 or 2)
#' @param boo.Networks Boolien specifiying whether the input topologies include networks (True) or not (False)
#' @param handle.GeneTrees Phylo object containing a list of the input gene trees
#' @param numeric.NumberOfReps Number of bootstrap replicates to analyze
#' @param numeric.MaxReticulations Number of maximum reticulating edges. Only used for the SOWH test with network topologies
#' @param string.PathDir String defining the path to a parent directory used for conduct KH_1 STAR test
#' @keywords Species tree, multispecies coalescent, phylogenetics
#' @return List Returns a list containing p-values and details for the assumed topology test
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
#'
#' ########################
#' # Run SpeciesTopoTestR #
#' ########################
#' handle.SpeciesTopologies <- list()
#' class(handle.SpeciesTopologies) <- "multiPhylo"
#' handle.SpeciesTopologies[[1]] <- handle.SpeciesTree1
#' handle.SpeciesTopologies[[2]] <- handle.SpeciesTree2
#'
#' list.SpeciesTopologies <- list()
#' list.SpeciesTopologies[[1]] <- write.tree(phy = handle.SpeciesTree1, file = "")
#' list.SpeciesTopologies[[2]] <- write.tree(phy = handle.SpeciesTree2, file = "")
#'
#'
#' Conduct.SpeciesTopoTestR(handle.Topologies2Test = list.SpeciesTopologies,
#'                          handle.GeneTrees = handle.Simulated_GeneTrees_T1,
#'                          numeric.NumberOfReps = 3,
#'                          string.Test = "SOWH",
#'                          numeric.Algorithm = 2,
#'                          boo.Networks = T,
#'                          string.PathDir = '~/Desktop/',
#'                          numeric.MaxReticulations = 1)
#'

############################
# Conduct.SpeciesTopoTestR #
############################
Conduct.SpeciesTopoTestR <- function(handle.Topologies2Test, handle.GeneTrees, numeric.NumberOfReps, string.Test, numeric.Algorithm, boo.Networks, string.PathDir, numeric.MaxReticulations){

  #####################################################
  # Define directory used for conduct SH_2N STAR test #
  #####################################################
  string.Dir_SpeciesTopoTestR <- paste0(string.Test, "_", numeric.Algorithm)

  if (boo.Networks == T){
    string.Dir_SpeciesTopoTestR <- paste0(string.Dir_SpeciesTopoTestR, "N")
  }

  ##############################################
  # Define directory used for SpeciesTopoTestR #
  ##############################################
  string.CurrentDir <- getwd()
  string.Path_Directory_SpeciesTopoTestR = paste0(string.PathDir, '/SpeciesTopoTestR_', string.Dir_SpeciesTopoTestR, "_", Sys.Date())
  unlink(string.Path_Directory_SpeciesTopoTestR, recursive = T)
  dir.create(string.Path_Directory_SpeciesTopoTestR, showWarnings = T, recursive = T)

  #######################
  # Define tests to run #
  #######################
  if (string.Dir_SpeciesTopoTestR == "KH_1"){

    #######
    # KH1 #
    #######
    handle.Results_SpeciesTopoTestR_KH_1 <- Conduct.KH_1(handle.SpeciesTree1 = handle.Topologies2Test[[1]],
                                                         handle.SpeciesTree2 = handle.Topologies2Test[[2]],
                                                         handle.InputGeneTrees = handle.GeneTrees,
                                                         numeric.NumberOfReps = numeric.NumberOfReps,
                                                         string.PathDir = string.Path_Directory_SpeciesTopoTestR)

    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_KH_1

  }
  if (string.Dir_SpeciesTopoTestR == "KH_1N"){

    #########
    # KH1_N #
    #########
    handle.Results_SpeciesTopoTestR_KH_1N <- Conduct.KH_1N(string.SpeciesNetwork1 =  handle.Topologies2Test[[1]],
                                                         string.SpeciesNetwork2 = handle.Topologies2Test[[2]],
                                                         handle.InputGeneTrees = handle.GeneTrees,
                                                         numeric.NumberOfReps = numeric.NumberOfReps,
                                                         string.PathDir = string.Path_Directory_SpeciesTopoTestR)
    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_KH_1N

  }
  if (string.Dir_SpeciesTopoTestR == "KH_2"){

    #######
    # KH2 #
    #######
    handle.Results_SpeciesTopoTestR_KH_2 <- Conduct.KH_2(handle.SpeciesTree1 = handle.Topologies2Test[[1]],
                                                         handle.SpeciesTree2 = handle.Topologies2Test[[2]],
                                                         handle.InputGeneTrees = handle.GeneTrees,
                                                         numeric.NumberOfReps = numeric.NumberOfReps,
                                                         string.PathDir = string.Path_Directory_SpeciesTopoTestR)
    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_KH_2

  }
  if (string.Dir_SpeciesTopoTestR == "KH_2N"){

    #########
    # KH2_N #
    #########
    handle.Results_SpeciesTopoTestR_KH_2N <- Conduct.KH_2N(string.SpeciesNetwork1 =  handle.Topologies2Test[[1]],
                                                           string.SpeciesNetwork2 = handle.Topologies2Test[[2]],
                                                           handle.InputGeneTrees = handle.GeneTrees,
                                                           numeric.NumberOfReps = numeric.NumberOfReps,
                                                           string.PathDir = string.Path_Directory_SpeciesTopoTestR)
    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_KH_2N

  }
  if (string.Dir_SpeciesTopoTestR == "KH_3"){

    #######
    # KH3 #
    #######
    handle.Results_SpeciesTopoTestR_KH_3 <- Conduct.KH_3(handle.SpeciesTree1 = handle.Topologies2Test[[1]],
                                                         handle.SpeciesTree2 = handle.Topologies2Test[[2]],
                                                         handle.InputGeneTrees = handle.GeneTrees,
                                                         numeric.NumberOfReps = numeric.NumberOfReps,
                                                         string.PathDir = string.Path_Directory_SpeciesTopoTestR)
    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_KH_3
  }
  if (string.Dir_SpeciesTopoTestR == "KH_3N"){

    #########
    # KH3_N #
    #########
    handle.Results_SpeciesTopoTestR_KH_3N <- Conduct.KH_3N(string.SpeciesNetwork1 =  handle.Topologies2Test[[1]],
                                                           string.SpeciesNetwork2 = handle.Topologies2Test[[2]],
                                                           handle.InputGeneTrees = handle.GeneTrees,
                                                           numeric.NumberOfReps = numeric.NumberOfReps,
                                                           string.PathDir = string.Path_Directory_SpeciesTopoTestR)
    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_KH_3N
  }


  if (string.Dir_SpeciesTopoTestR == "SH_1"){

    #######
    # SH1 #
    #######
    handle.Results_SpeciesTopoTestR_SH_1 <- Conduct.SH_1(handle.SpeciesTopologies = handle.Topologies2Test,
                                                         handle.InputGeneTrees = handle.GeneTrees,
                                                         numeric.NumberOfReps = numeric.NumberOfReps,
                                                         string.PathDir = string.Path_Directory_SpeciesTopoTestR)

    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_SH_1

  }
  if (string.Dir_SpeciesTopoTestR == "SH_1N"){

    #######
    # SH1 #
    #######
    handle.Results_SpeciesTopoTestR_SH_1N <- Conduct.SH_1N(list.SpeciesTopologiesNetworks = handle.Topologies2Test,
                                                         handle.InputGeneTrees = handle.GeneTrees,
                                                         numeric.NumberOfReps = numeric.NumberOfReps,
                                                         string.PathDir = string.Path_Directory_SpeciesTopoTestR)

    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_SH_1N

  }

  if (string.Dir_SpeciesTopoTestR == "SH_2"){

    #######
    # SH2 #
    #######
    handle.Results_SpeciesTopoTestR_SH_2 <- Conduct.SH_2(handle.SpeciesTopologies = handle.Topologies2Test,
                                                         handle.InputGeneTrees = handle.GeneTrees,
                                                         numeric.NumberOfReps = numeric.NumberOfReps,
                                                         string.PathDir = string.Path_Directory_SpeciesTopoTestR)

    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_SH_2

  }
  if (string.Dir_SpeciesTopoTestR == "SH_2N"){

    #######
    # SH2 #
    #######
    handle.Results_SpeciesTopoTestR_SH_2N <- Conduct.SH_2N(list.SpeciesTopologiesNetworks = handle.Topologies2Test,
                                                           handle.InputGeneTrees = handle.GeneTrees,
                                                           numeric.NumberOfReps = numeric.NumberOfReps,
                                                           string.PathDir = string.Path_Directory_SpeciesTopoTestR)

    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_SH_2N

  }


  if (string.Dir_SpeciesTopoTestR == "SOWH_1"){

    #########
    # SOWH1 #
    #########
    handle.Results_SpeciesTopoTestR_SOWH_1 <- Conduct.SOWH_1(handle.SpeciesTree_1 = handle.Topologies2Test[[1]],
                                                         handle.InputGeneTrees = handle.GeneTrees,
                                                         numeric.NumberOfReps = numeric.NumberOfReps,
                                                         string.PathDir = string.Path_Directory_SpeciesTopoTestR)

    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_SOWH_1

  }
  if (string.Dir_SpeciesTopoTestR == "SOWH_1N"){

    #########
    # SOWH1N #
    #########
    handle.Results_SpeciesTopoTestR_SOWH_1N <- Conduct.SOWH_1N(string.SpeciesNetwork_1 = handle.Topologies2Test[[1]],
                                                             handle.InputGeneTrees = handle.GeneTrees,
                                                             numeric.NumberOfReps = numeric.NumberOfReps,
                                                             numeric.MaxReticulations = numeric.MaxReticulations,
                                                             string.PathDir = string.Path_Directory_SpeciesTopoTestR)

    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_SOWH_1N

  }

  if (string.Dir_SpeciesTopoTestR == "SOWH_2"){

    #########
    # SOWH2 #
    #########
    handle.Results_SpeciesTopoTestR_SOWH_2 <- Conduct.SOWH_2v3(handle.SpeciesTree_1 = handle.Topologies2Test[[1]],
                                                             handle.InputGeneTrees = handle.GeneTrees,
                                                             numeric.NumberOfReps = numeric.NumberOfReps,
                                                             string.PathDir = string.Path_Directory_SpeciesTopoTestR)

    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_SOWH_2

  }
  if (string.Dir_SpeciesTopoTestR == "SOWH_2N"){

    #########
    # SOWH2N #
    #########
    handle.Results_SpeciesTopoTestR_SOWH_2N <- Conduct.SOWH_2N(string.SpeciesNetwork_1 = handle.Topologies2Test[[1]],
                                                               handle.InputGeneTrees = handle.GeneTrees,
                                                               numeric.NumberOfReps = numeric.NumberOfReps,
                                                               numeric.MaxReticulations = numeric.MaxReticulations,
                                                               string.PathDir = string.Path_Directory_SpeciesTopoTestR)

    handle.SpeciesTopoTestR_Results <- handle.Results_SpeciesTopoTestR_SOWH_2N

  }





 return(handle.SpeciesTopoTestR_Results)
}
