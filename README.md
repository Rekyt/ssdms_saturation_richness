
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Is prediction of species richness from Stacked Species Distribution Models biased by habitat saturation?

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/224002794.svg)](https://zenodo.org/badge/latestdoi/224002794)
[![Launch Rstudio
Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Rekyt/ssdms_saturation_richness/master?urlpath=rstudio)
<!-- badges: end -->

This repository contains the data and code for our paper:

> Grenié M., Violle C, Munoz F. * Is prediction of species richness from
> Stacked Species Distribution Models biased by habitat saturation?*.
> accepted in *Ecological Indicators*.

## How to cite

Please cite this compendium as:

> Grenié M., Violle C, Munoz F., (2019). *Compendium of R code and data
> for Is prediction of species richness from Stacked Species
> Distribution Models biased by habitat saturation?*. Accessed 02 déc.
> 2019. Online at <https://doi.org/10.5281/zenodo.3552836>

## 🔧 How to download or install

You can download the compendium as a zip from from this URL:
</archive/master.zip>

Or you can install this compendium as an R package,
\`cssdms.saturation.richness, from GitHub with:

``` r
# install.packages("devtools")
remotes::install_github("Rekyt/ssdms_saturation_richness")
```

## 💻 How to run the analyses

This compendium uses [`drake`](https://docs.ropensci.org/drake/) to make
analyses reproducible. To redo the analyses and rebuild the manuscript
run the following lines (from the `comsat` folder):

``` r
# install.packages("devtools")
pkgload::load_all()  # Load all functions included in the package
make(saturation_workflow())  # Run Analyses
```

Beware that some code make time a long time to run, and it may be useful
to run analyses in parallel.

\#\#You can run the analyses by clicking on the `Binder` badge:
[![Launch Rstudio
Binder](http://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/Rekyt/ssdms_saturation_richness/master?urlpath=rstudio)

Dependencies

As noted in the `DESCRPTION` files this project depends on:

  - [`virtualspecies`](https://cran.r-project.org/package=virtualspecies),
    to simulate species;
  - [`drake`](https://cran.r-project.org/package=drake), to execute a
    reproducible workflow;
  - the `tidyverse` (`dplyr`, `ggplot2`, `purrr`, and `tidyr`) for data
    wrangling;
  - [`ggpubr`](https://cran.r-project.org/package=ggpubr) to customize
    plot
