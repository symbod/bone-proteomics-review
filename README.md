# Meta-analysis of Proteomics Data from Osteoblasts, Bone and Blood

## Overview

This repository contains the code and data for the paper titled
**"Meta-analysis of Proteomics Data from Osteoblasts, Bone and Blood:
Insights into Druggable Targets, Active Factors and Potential Biomarkers
for Bone Biomaterial Design"** by Schmidt, Adamowicz & Arend et al.

## Setup Instructions

### Prerequisites

Make sure you have the following packages installed in R:

-   `rmarkdown`
-   `gprofiler2`
-   `data.table`
-   `ggplot2`
-   `igraph`
-   `reticulate`
-   `plotly`
-   `tidyr`
-   `dplyr`
-   `GGally`
-   `ComplexHeatmap`
-   `circlize`
-   `readxl`
-   `ggpattern`
-   `stringr`
-   `openxlsx`
-   `ggtext`
-   `ggpubr`
-   `scatterpie`
-   `ggsankey` (install from GitHub)

### Installation

``` r
# Install required R packages
required_packages <- c("rmarkdown", "gprofiler2", "data.table", "ggplot2", "igraph", "reticulate", "plotly", "tidyr", "dplyr", "GGally", "ComplexHeatmap", "circlize", "readxl", "ggpattern", "stringr", "openxlsx", "ggtext", "ggpubr", "scatterpie")
for(package in required_packages){
  if(!require(package,character.only = TRUE, quietly = TRUE)) install.packages(package, dependencies = TRUE, quietly = TRUE)
  library(package, character.only = TRUE, quietly = TRUE)
}

# Install ggsankey from GitHub
remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
```

### Virtual Environment Setup

``` r
# Set up Python virtual environment
reticulate::virtualenv_create("bbr")
reticulate::py_install("drugstone", pip = TRUE, ignore_installed=TRUE, envname = "bbr")
reticulate::use_virtualenv("bbr")

# Import Python package
ds <- reticulate::import("drugstone")
ds$print_license()
ds$accept_license()
```

### Running the Analysis

To run the analysis, follow the workflow described in the `Workflow.Rmd`
file. The workflow includes data harmonization, intersection analysis,
network enrichment, and visualization steps.

## Data

The data required for the analysis is stored in the `data` directory.
Ensure all input files are placed in the correct directories before
running the scripts.

## Visualizations

The repository includes various scripts to generate visualizations such
as Venn diagrams, Sankey plots, heatmaps, and network graphs. Refer to
the `Workflow.Rmd` file for detailed instructions on generating these
visualizations.

## Contributing

Contributions are welcome! Please create a pull request or open an issue
to discuss any changes or improvements.

## License

This project is licensed under the  GPL-3.0 License. See the LICENSE file for
details.

## Contact

For any questions or inquiries, please contact Klaudia Adamowicz or Lis
Arend.
