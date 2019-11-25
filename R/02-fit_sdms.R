# This script contains the functions to fit SDMs based on data

# Functions for regular sdm models ---------------------------------------------

compute_threshold = function(df, first_col = "presence", second_col = "pred") {
    r = SDMTools::optim.thresh(df %>%
                      dplyr::select(!!first_col) %>%
                      unlist(),
                  df %>%
                      dplyr::select(!!second_col) %>%
                      unlist())

    max(r$`max.sensitivity+specificity`)
}

add_prediction = function(df, model, prediction_col = "pred",
                          ...) {
    df %>%
        dplyr::mutate(!!prediction_col := predict(model, df, ...))
}

#' Standard way of building SDMs
#' @import tidyverse
#' @export
build_sdm = function(given_hab_sat) {
    given_hab_sat %>%
        tidyr::unnest(obs) %>%
        tidyr::nest(-hab_sat, -hab_type, -species, .key = "obs") %>%
        # SDM + variable extraction
        dplyr::mutate(sp_glm = obs %>%
                   purrr::map(~glm(presence ~ cell + I(cell^2), data = .x,
                            control = list(maxit = 60),
                            family = binomial(link = "logit"))),
               sp_gam = obs %>%
                   purrr::map(~mgcv::gam(presence ~ s(cell, bs = "cr"),
                                         data = .x)),
               # Compute species prevalence in sample data
               samp_prev = obs %>%
                   purrr::map_dbl(~sum(.x$presence)/nrow(.x))) %>%
        # Filter out species totally absent from dataset
        dplyr::filter(samp_prev > 0) %>%
        dplyr::mutate(
            # Add predicted probability column into data
            obs = obs %>%
                purrr::map2(sp_glm,
                     ~add_prediction(.x, .y, "pred_glm", type = "response")) %>%
                purrr::map2(sp_gam,
                     ~add_prediction(.x, .y, "pred_gam", type = "response")),
            # Compute probability threshold
            thresh_glm = purrr::map_dbl(obs, compute_threshold,
                                 second_col = "pred_glm"),
            thresh_gam = purrr::map_dbl(obs, compute_threshold,
                                 second_col = "pred_gam"),
            # Maximum predicted probability
            max_prob = purrr::map_dbl(obs, ~max(.x$pred_glm)))
}
