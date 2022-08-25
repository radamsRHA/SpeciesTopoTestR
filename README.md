
---
# SpeciesTopoTestR: Likelihood-based hypothesis tests of species topologies in R
**NOTE See the file https://github.com/radamsRHA/SpeciesTopoTestR/blob/master/SpeciesTopoTestR.pdf for detailed instructions**

## Installing R package SpeciesTopoTestR from github
The R package SpeciesTopoTestR is freely available to download and distribute from github <https://github.com/radamsRHA/SpeciesTopoTestR/>. To install and load SpeciesTopoTestR, you must first install the R package `devtools`, 

```
install.packages("devtools")
```
Now using devtools we can install `SpeciesTopoTestR` from github:

```
library(devtools)
install_github("radamsRHA/SpeciesTopoTestR")
library(SpeciesTopoTestR) # Load package 
```
`SpeciesTopoTestR` also requires the following dependencies to be installed:

```
install.packages('ape')  
```

And the following external software is required to be install and in your $PATH:  
* STELLS2: https://github.com/yufengwudcs/STELLS2  
* PHYLONET: https://bioinfocs.rice.edu/phylonet  
* HYBRID-LAMBDA: https://github.com/hybridLambda/hybrid-Lambda  

These two software packages can also be found in the https://github.com/radamsRHA/SpeciesTopoTestR/inst/extdata directory.   
For OSX (mac) this means that the executable `stells-v2-1-0-1-mac` must be installed and in your path.  
For linux, the executable `stells-v2-1-0-1-linux64` must be installed and in your path.   
Additionally, the executable `hybrid-Lambda` must also be in your path if you want to test network hybridizations. 

This can also be added manually to the `$PATH` using the following example command lines (opening the terminal and typing the following commands):  

`export PATH=$PATH:/home/adamsr/PROGRAMS/STELLS2-master/`  
`export PATH=$PATH:/home/adamsr/PROGRAMS/hybrid-Lambda-dev/`  

where `/home/adamsr/PROGRAMS/` is replaced with the system path to your install directories. 

**IMPORTANT: Both the executables for STELLS2 must be found in your $PATH command line argument**  
**IMPORTANT: For network testing, the executables for HYBRID-LAMBDA must be found in your $PATH command line argument**  

To begin using `SpeciesTopoTestR` try using the examples associated with each function. Importantly, the three tests implemented by SpeciesTopoTest R (KH\* SH\*, and SOWH\*) can be conducted using the primary function `Conduct.SpeciesTopoTestR` with specific options. See examples below:

## Example: KH* Test

We can use `Conduct.SpeciesTopoTestR` with option "KH" to run the KH* test for a given datasets. First, let's load the R package `SpeciesTopoTestR` and its dependancies:

```
################
# Load depends #
################
library(ape)
library(phytools)
library(SpeciesTopoTestR)
```

Now, let's specify a target topology S1 and load some example gene trees (also simulated from S11) for the KH* test:


```
#################################
# Generate example species trees #
#################################
handle.SpeciesTree1 <- read.tree(text = "(A:1.0,(B:0.5,(C:1.0,D:1.0):1.5):0.5);")

handle.Simulated_GeneTrees <-  system.file("extdata", "ExampleGeneTreeSet.tree", package="SpeciesTopoTestR")
handle.Simulated_GeneTrees <- read.tree(file = handle.Simulated_GeneTrees)
```


Next, let's conduct the KH* test using an alternative topology S2

```
###########################
# Alternative topology S2 #
###########################
handle.SpeciesTree2 <- read.tree(text = "(C:1.0,(B:0.5,(A:1.0,D:1.0):1.5):0.5);")

#########################################
# Combined S1 and S2 into a single list #
#########################################
handle.SpeciesTopologies <- list()
class(handle.SpeciesTopologies) <- "multiPhylo"
handle.SpeciesTopologies[[1]] <- handle.SpeciesTree1
handle.SpeciesTopologies[[2]] <- handle.SpeciesTree2

################
# RUN KH* test #
################

RESULTS <- Conduct.SpeciesTopoTestR(handle.Topologies2Test = handle.SpeciesTopologies,
                         handle.GeneTrees = handle.Simulated_GeneTrees,
                         numeric.NumberOfReps = 100,
                         string.Test = "KH", # RUN KH* test
                         numeric.Algorithm = 1,  # RUN KH* test with Algorithm 2
                         boo.Networks = F, # No network topologies
                         string.PathDir = '~/Desktop/')
                         
RESULTS$TwoSided_Pvalue # The two-tailed p-value is the final p-value

```

## Example: SH* Test

We can use `Conduct.SpeciesTopoTestR` with option "SH" to run the SH* test for a given datasets. 


We will define three total topologies that will be placed in a `multiPhylo` list object `handle.SpeciesTopologies` alongside the target topology `handle.SpeciesTree1`:

```
#################################
# Generate example species trees #
#################################
handle.SpeciesTree2 <- read.tree(text = "(C:1.0,(B:0.5,(A:1.0,D:1.0):1.5):0.5);")
handle.SpeciesTree3 <- read.tree(text = "(D:1.0,(B:0.5,(A:1.0,C:1.0):1.5):0.5);")
handle.SpeciesTopologies <- list()
class(handle.SpeciesTopologies) <- "multiPhylo"

handle.SpeciesTopologies[[1]] <- handle.SpeciesTree3
handle.SpeciesTopologies[[2]] <- handle.SpeciesTree2
handle.SpeciesTopologies[[3]] <- handle.SpeciesTree1
```

Let's conduct the SH* test:

```
RESULTS <- Conduct.SpeciesTopoTestR(handle.Topologies2Test = handle.SpeciesTopologies,
                         handle.GeneTrees = handle.Simulated_GeneTrees,
                         numeric.NumberOfReps = 100,
                         string.Test = "SH", # RUN SH* test
                         numeric.Algorithm = 1,  # RUN SH* test with Algorithm 1
                         boo.Networks = F, # No network topologies
                         string.PathDir = '~/Desktop/')
                         
RESULTS$vector.Pvalues 

```

## Example: SWOH* Test

We can use `Conduct.SpeciesTopoTestR` with option "SOWH" to run the SWOH* test for a given datasets. 


Next, let's run the SWOH* test using this focal topology:

```
#####################################################
# Build multiPhylo list to hold the single topology #
#####################################################
handle.SpeciesTopologies <- list()
class(handle.SpeciesTopologies) <- "multiPhylo"
handle.SpeciesTopologies[[1]] <- handle.SpeciesTree1

RESULTS <- Conduct.SpeciesTopoTestR(handle.Topologies2Test = handle.SpeciesTopologies,
                         handle.GeneTrees = handle.Simulated_GeneTrees,
                         numeric.NumberOfReps = 100,
                         string.Test = "SOWH", # RUN SOWH* test
                         numeric.Algorithm = 1,  # RUN SOWH* test with Algorithm 1
                         boo.Networks = F, # No network topologies
                         string.PathDir = '~/Desktop/')
RESULTS$numeric.pValue

```

We can also run the SOWH* test using a different topology to compare the results:

```
#####################################################
# Build multiPhylo list to hold the single topology #
#####################################################
handle.SpeciesTopologies <- list()
class(handle.SpeciesTopologies) <- "multiPhylo"
handle.SpeciesTopologies[[1]] <- handle.SpeciesTree2

RESULTS <- Conduct.SpeciesTopoTestR(handle.Topologies2Test = handle.SpeciesTopologies,
                                    handle.GeneTrees = handle.Simulated_GeneTrees,
                                    numeric.NumberOfReps = 100,
                                    string.Test = "SOWH", # RUN SOWH* test
                                    numeric.Algorithm = 1,  # RUN SOWH* test with Algorithm 1
                                    boo.Networks = F, # No network topologies
                                    string.PathDir = '~/Desktop/')
RESULTS$numeric.pValue
```
## NOTE: formatting networks for SpeciesTopoTestR

SpeciesTopoTestR uses the same format as the program PHYLONET (https://bioinfocs.rice.edu/phylonet) for network-based species topologies. Instead of placing network topologies in a `multiPhylo` list object, we simply uses the strings describing the topologies. 

We can use the following example by placing two topology strings `string.SpeciesNetwork` and `string.SpeciesNetwork_2` into a single list `list.SpeciesTopologies`:

```
####################################
# Generate example species network #
####################################
string.SpeciesNetwork <- "(((((C:1.0,D:1.0):1)#H1:0::0.25,A:1.0):2,B:1.0):2,#H1:0::0.75);"
string.SpeciesNetwork_2 <- "((A:1,B:1):1,(C:1,D:1):1);"

list.SpeciesTopologies <- list()
list.SpeciesTopologies[[1]] <- string.SpeciesNetwork
list.SpeciesTopologies[[2]] <- string.SpeciesNetwork_2

Conduct.SpeciesTopoTestR(handle.Topologies2Test = list.SpeciesTopologies,
                         handle.GeneTrees = handle.Simulated_GeneTrees,
                         numeric.NumberOfReps = 3,
                         string.Test = "KH",
                         numeric.Algorithm = 3,
                         boo.Networks = T,
                         string.PathDir = '~/Desktop/',
                         numeric.MaxReticulations = 1)

```




