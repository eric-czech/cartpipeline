# Mock CAR T-Cell Target Pipeline Project

This project prototypes a (naive) pipeline for generating CAR T-Cell therapy targets for a variety of different types of cancers.  

[CAR T-Cell Therapy](https://www.cancer.gov/about-cancer/treatment/research/car-t-cells) targets are typically proteins expressed on the surface of tumorous cells (independent of MHC presentation) with certain properties that make them appropriate for artificially induced immune responses.  More specifically, this pipeline would ideally identify proteins meeting the following criteria ([Atanackovic et al 2016](http://onlinelibrary.wiley.com/doi/10.1111/bjh.13889/full)):

1. Must be expressed on cell surface (must be accessible by extracellular chimeric receptors)
2. Should be expressed in as few non-cancerous cells as possible
3. Should be expressed by as much of the afflicted population as possible
4. Should be homogeneously expressed by cancerous cells of any one patient
5. Should be critical to the function of the cell so that it is less likely to be downregulated by selective pressure induced by immunotherapy

Currently, this project at least tries to identify targets according to the first 3 criteria using [TCGA](https://cancergenome.nih.gov/) expression data (via [cBioPortal](http://www.cbioportal.org/study?id=brca_tcga#summary)) and the [Human Protein Atlas](http://www.proteinatlas.org/) (HPA).  HPA data is used to identify/prioritize candidates meeting criteria \#1 and \#2 since it provides common locations of proteins for genes (e.g. in cellular membranes, intracellularly, within mitochondria, etc) as well as tissue specificity scores denoting how specific to a certain tissue type any one protein is.  TCGA RNA-seq data then provides further evidence for \#2 since RNA-seq data is differential based on tumor/normal pairs as well as for \#3 by allowing for examination of over-expression within the patient population.

## Overview

This pipeline has 3 major parts:

1. Python Scripts - Used to collect and process data
2. [Ketrew](https://github.com/hammerlab/ketrew)/[Coclobas](https://github.com/hammerlab/coclobas) - Used for pipeline task scheduling and management
3. Docker - Used by Ketrew to containerize individual tasks

In no particular order, some noteworthy details on parts above are:

- 
