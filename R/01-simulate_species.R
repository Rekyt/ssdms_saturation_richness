# This script contains function to simulate species

#' Simulating environmental layer
#'
#' This function simulates a single linear gradient from a raster
#' @param length length of the linear gradient (default = 2000)
#'
#' @importFrom raster raster
#' @export
create_env_raster = function(length = 2000) {
    art_raster = raster(matrix(1:2000))
    names(art_raster) = "unidim"

    return(art_raster)
}

#' Simulate Random Species with gaussian responses
#'
#' Wrapper around the function `generateRandomSp()` of `virtualspecies`. This
#' function generates species with gaussian environmental response, with random
#' response curve parameters
#'
#' @param n number of species to simulate
#' @param env_raster environmental raster
#'
#' @importFrom virtualspecies generateRandomSp
#' @export
simulate_species = function(n = 100, env_raster, given_seed = 27022018) {
    set.seed(given_seed)

    virt_species = lapply(seq_len(n),
                         function(x) {
                             generateRandomSp(env_raster,
                                              approach = "response",
                                              relations = "quadratic",
                                              realistic.sp = TRUE,
                                              plot = FALSE)
                         })

    names(virt_species) = paste0("virt_", seq_along(virt_species))

    return(virt_species)
}

#' Randomly sample species presence-absence
#'
#' This function sample presences absences from a virtual species list on a
#' given number of sites. It outputs a data.frame containing presences and
#' absences of species all in the same sites.
#'
#' @param species_list a list of virtual species produced by the
#'                     `virtualspecies` package
#' @param n_cell number of cells to sample
#'
#' @import raster
#' @importFrom sp coordinates
#' @import tidyverse
#' @export
sample_species_pa = function(species_list, virt_raster, n_cell = 300) {
    species_list %>%
        # Stack all species PA raster
        purrr::map(~.x$pa.raster) %>%
        stack() %>%
        # Sample the same sites for all species
        sampleRandom(n_cell, cell = TRUE) %>%
        dplyr::as_data_frame() %>%
        # Add cell coordinates as two separate columns
        dplyr::mutate(coords = cell %>%
                          purrr::map(~raster::xyFromCell(virt_raster, .x)) %>%
                          purrr::map(as.data.frame)
        ) %>%
        tidyr::unnest(coords)
}

# Affect Habitat saturation
change_alt_hab_sat = function(species_list, hab_sat_level) {
    species_prob = species_list %>%
        purrr::map_dfr(~.x$details$parameters$unidim$args %>%
                           {do.call(virtualspecies::quadraticFun,
                                    c(list(1:2000), .))} %>%
                           scales::rescale(to = c(0, 1)) %>%
                           data.frame(cell = 1:2000, env_suitab = .,
                                      hab_sat = hab_sat_level) %>%
                           mutate(env_sat = env_suitab * hab_sat,
                                  presence = sapply(env_sat, function(x) {
                                      prob = ifelse(x > 1, 1, x)
                                      rbinom(1, 1, prob)
                                  })),
                       .id = "species") %>%
        mutate(hab_type = "alt") %>%
        nest(-hab_sat, -hab_type, .key = "obs")
}
