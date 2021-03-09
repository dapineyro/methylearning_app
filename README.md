
# Methylearning Shiny App
 
* Author: David Piñeyro
* License: GPL-3
* Date: 2018-06-01

This is the companion Shiny app for the methylearning package. 

## Installation

In order to run this app, you need several R packages 
previously installed and working in your R session. To install
all of them use:

```{r}
# install CRAN packages

install.packages(shiny)
install.packages(data.table)
install.packages(parallel)
install.packages(xtable)
install.packages(e1071)
install.packages(shinythemes)

# install Bioconductor and Bioconductor packages
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite()
biocLite("GEOquery")
```

You will also need to install `methylearning` package and all its dependencies.

```{r}
# In your R session, type:
# Set one level up methylearning folder.
setwd("path/to/folder/containing/methylearning/folder")
devtools::install("methylearning")
```

Dependencies should be installed automatically. If you experience any problem, it would be probably related to some system packages lacking or wrongly configured. For instance, in Ubuntu 17.10, apart from `r-base` and `r-base-dev` system packages and the `devtools` `R` package, you will also need some extra system packages. To install them, in the command shell type (root access required):

```
sudo apt-get update
sudo apt-get upgrade
sudo apt-get install libcurl4-openssl-dev libssl-dev libxml2-dev
# install and configure Java.
sudo apt-get install default-jre default-jdk
sudo R CMD javareconf
# install MySQL server.
sudo apt-get install libmariadbclient-dev
```

## Run the app.
Start a **new** R session and write:

```{r}
# Change your working directory to the methylearning_app/ directory.
setwd("/path/to/methylearning_app")
shiny::runApp("app.R")
```

**IMPORTANT NOTE:** `methylearning` package make use of many other packages. It is possible that you reach your maximum number of allowed DLLs in your `R` session (i.e. more than 100). If you experience any error realated to a failed loaded package please, start `methylearning` Shiny APP from a new `R` session. If you still experience problems type the code shown below in your command shell (Linux or similar OS) to increase your DLL limit.

```
# Create an .Renviron in your home directory, which contains the 
# new value for the R_MAX_NUM_DLLS variable (it will be increased to 150).
echo "R_MAX_NUM_DLLS=150" > ~/.Renviron
```
