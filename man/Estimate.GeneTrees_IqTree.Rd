% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Estimate.GeneTrees_IqTree.R
\name{Estimate.GeneTrees_IqTree}
\alias{Estimate.GeneTrees_IqTree}
\title{Estimate.GeneTrees_IqTree: function to estimate a set of gene trees for an input list of nexus alignments (simulated with Simulate.Gene.Alignments)}
\usage{
Estimate.GeneTrees_IqTree(
  list.Simulated.Nexus.Alignments,
  path.PathParentDir,
  string.Model
)
}
\arguments{
\item{list.Simulated.Nexus.Alignments}{List of nexus alignments that have been simulated using Simulate.Gene.Alignments}

\item{path.PathParentDir}{Path to parent directory used for gene tree estimation}

\item{string.Model}{String to specify substitution model for IQ-tree}
}
\description{
This function returns (1) a list of ML gene tree estimates obtained via IQ-Tree
}
\examples{



################
# Load depends #
################
library(SpeciesTopoTestR)
library(ape)

########################################################
# Define species topologies for simulation and testing #
########################################################
handle.SpeciesTree1 <- read.tree(text = "(A:3,(B:2,(C:1,D:1):1):1);")
handle.SpeciesTree2 <- read.tree(text = "(A:3,(D:2,(C:1,B:1):1):1);")


##########################
# Simulate gene tree set #
##########################
handle.Simulated_GeneTrees_T1 <- Simulate.GeneTrees_From_SpeciesTree(handle.SpeciesTree = handle.SpeciesTree1,
                                                                     string.PathDir = '~/Desktop/',
                                                                     numeric.NumberOfGeneTrees = 100)
############################################
# Convert to mutation units branch lengths #
############################################
handle.Converted_GeneTrees_T1 <- list()
class(handle.Converted_GeneTrees_T1) <- "multiPhylo"

for (i in 1:length(handle.Simulated_GeneTrees_T1)){

  handle.Simulated_GeneTrees_T1_i <- handle.Simulated_GeneTrees_T1[[i]]
  handle.Simulated_GeneTrees_T1_i$edge.length = handle.Simulated_GeneTrees_T1_i$edge.length* 0.00005
  handle.Converted_GeneTrees_T1[[i]] <- handle.Simulated_GeneTrees_T1_i
}

################################
# Simulate sequence alignments #
################################
vector.pi <- c(0.3, 0.2, 0.2, 0.3)
names(vector.pi) <- c("A", "C", "G", "T")
list.SimulatedAlignments <- Simulate.GeneSequenceAlignments(list.Simulated.Gene.Tree.Set = handle.Converted_GeneTrees_T1,
                                                            path.PathParentDir = '~/Desktop/',
                                                            numeric.Ratio = 4.6,
                                                            numeric.LocusLength = 1000,
                                                            vector.BaseFrequencies = vector.pi)

#######################
# Estimate gene trees #
#######################
handle.EstimatedGeneTrees <- Estimate.GeneTrees_IqTree(list.Simulated.Nexus.Alignments = list.SimulatedAlignments$list.Simulated.Nexus.Alignments,
                                                       path.PathParentDir = '~/Desktop/',
                                                       string.Model = "HKY")

################
# Conduct KH_2 #
################
handle.Topologies <- list()
class(handle.Topologies) <- "multiPhylo"
handle.Topologies[[1]] <- handle.SpeciesTree1
handle.Topologies[[2]] <- handle.SpeciesTree2

Conduct.SpeciesTopoTestR(handle.Topologies2Test = handle.Topologies,
                         handle.GeneTrees = handle.EstimatedGeneTrees$list.MLE.GeneTrees,
                         numeric.NumberOfReps = 1000,
                         string.Test = "KH",
                         numeric.Algorithm = 2,
                         boo.Networks = F,
                         string.PathDir = '~/Desktop/')
}
\keyword{bayesian}
\keyword{coalescent}
\keyword{empirical}
\keyword{estimation,}
\keyword{inference,}
\keyword{multispecies}
\keyword{species}
\keyword{tree}
