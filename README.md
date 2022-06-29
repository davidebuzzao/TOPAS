# TOPAS (TOP-down Attachment of Seeds) ###

TOPAS is an iterative approach which aims at connecting the largest number of seed nodes in a top-down fashion through connectors which guarantee the highest flow of a Random Walk with Restart in a network of functional associations.

* Version 1.0

### Step 1: Installation

TOPAS is implemented in R (≥4.1), and it was tested in a version-controlled conda environment. 

Download this repository and either execute A or B:

A.  follow the instructions to install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) and then execute the following to generate a conda environment with all you need to run TOPAS:
```
conda env create -f environment.yml
```

B.  manually install the following libraries inside a R (≥4.1) environment:

  * r-tidyverse=1.3.1
  * r-igraph=1.2.11
  * r-dnet=1.1.7
  * bioconductor-biocparallel=1.28.0
  * r-r.utils=2.11.0
  * r-dosnow=1.0.20
  * r-foreach=1.5.1

### Step 2: Program Execution

You can copy this R snippet into your local computer to produce a disease module as an example. You may want to change the path to your local repository.

```{r}
#!/usr/bin/Rscript

setwd('~/topas') 
source('src/TOPAS.R')

## Load network
network_example = 
  readr::read_tsv(paste0('example/FunCoup5_pfc80.tsv.gz'), 
                  col_names = c('V1','V2'), 
                  col_types = 'cc')
## Seed gene set
seeds_example = 
  readr::read_table('example/adrenal_gland_diseases.txt', 
                    col_names = 'V',
                    col_types = 'c')

module_example = 
  TOPAS(
    network = network_example,
    seeds = seeds_example$V,
    expansion_steps = 2,
    cores = 4
    )
```

*NOTE*: TOPAS is only designed to process undirected networks!

### Contacts ###

* Davide Buzzao (davide.buzzao@scilifelab.se)
* Erik L.L. Sonnhammer (erik.sonnhammer@scilifelab.se)