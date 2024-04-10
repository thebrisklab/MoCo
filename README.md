# MoCo

MoCo is an R package designed for motion control in MRI studies. Please note that MoCo is now a private repository. To install MoCo, follow the steps below:

-   Download the Repository: Obtain the zip file for the MoCo repository.

-   Unzip the File: Extract the contents of the zip file.

-   Navigate to MoCo Folder: Access the MoCo folder in your file system.

-   Install MoCo:

```{r}
# Load necessary libraries
library(devtools)
library(roxygen2)

# Set working directory to MoCo folder
# Please change the path based on your computer
setwd("./MoCo")

# Generate documentation
devtools::document()

# Load the package
devtools::load_all()

# Install the package
devtools::install()

library(MoCo)
```
