MoCo: **Mo**tion-**Co**ntrolled Brain Phenotype Differences Between Groups <img src="fig/MoCo.png" width="120" align="right"/>
===================================================

MoCo is an R package designed to remove motion artifacts in brain phenotype analysis. Please note that MoCo is still under development.

## Installation

Ensure you have the following R packages installed:

-   SuperLearner
-   haldensify
-   devtools

You can install them by running the following code:

``` r
if(!require(c("SuperLearner","haldensify", "devtools"))){
    install.packages(c("SuperLearner","haldensify", "devtools"))
}
```

Then, you can install MoCo from GitHub using the following code:

```{r}
library(devtools)
install_github("thebrisklab/MoCo")

library(MoCo)
```

## Tutorial

The `moco()` is the main function of MoCo. The input and output of the function is in the figure below.

<img src="fig/input_output.png" width="300" align="center"/>

