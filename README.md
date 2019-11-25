
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Is prediction of species richness from Stacked Species Distribution Models biased by habitat saturation?

This repository contains the data and code for our paper:

> Grenié M., Violle C, Munoz F. * Is prediction of species richness from
> Stacked Species Distribution Models biased by habitat saturation?*.
> accepted in *Ecological Indicators*.

## How to cite

Please cite this compendium as:

> Grenié M., Violle C, Munoz F., (2019). *Compendium of R code and data
> for Is prediction of species richness from Stacked Species
> Distribution Models biased by habitat saturation?*. Accessed 25 nov.
> 2019. Online at <https://doi.org/xxx/xxx>

## 🔧 How to download or install

You can download the compendium as a zip from from this URL:
</archive/master.zip>

Or you can install this compendium as an R package, `comsat`, from
GitHub with:

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

## Dependencies

As noted in the `DESCRPTION` files this project depends on:

  - [`virtualspecies`](https://cran.r-project.org/package=virtualspecies),
    to simulate species;
  - [`drake`](https://cran.r-project.org/package=drake), to execute a
    reproducible workflow;
  - the `tidyverse` (`dplyr`, `ggplot2`, `purrr`, and `tidyr`) for data
    wrangling;
  - [`ggpubr`](https://cran.r-project.org/package=ggpubr) to customize
    plot
