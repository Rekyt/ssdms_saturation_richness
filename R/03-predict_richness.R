# This file contains function to predict species richness at each site

#' Predict richness from sdm table
#'
#' @import dplyr tidyr
#' @importFrom stringr str_extract
#' @export
predict_richness = function(sdm_df) {

    # Initial Grouped Data Frame
    by_cell_df = sdm_df %>%
        unnest(obs) %>%
        group_by(cell)

    # Observed Species Richness data.frame
    obs_rich_df = by_cell_df %>%
        summarise(hab_sat = mean(hab_sat),
                  obs_rich = sum(presence))

    # data.frame containing predicted species richness from summed probabilities
    pred_sum_df = by_cell_df %>%
        dplyr::select(cell, contains("pred")) %>%
        gather_by("pred") %>%
        mutate(pred_type = pred_type %>%
                   stringr::str_extract("[:alpha:]+$") %>%
                   paste0("pred_sum_", .)) %>%
        group_by(cell, pred_type) %>%
        summarise(pred_sum = sum(pred_value)) %>%
        spread(pred_type, pred_sum)

    # data.frame containing predicted species richness from thresholded prob.
    pred_thresh_df = by_cell_df %>%
        dplyr::select(hab_sat, cell, contains("thresh"),
                      contains("pred")) %>%
        gather_by("pred") %>%
        gather_by("thresh") %>%
        mutate(pred_type = pred_type %>%
                   stringr::str_extract("[:alpha:]+$"),
               thresh_type = thresh_type %>%
                   stringr::str_extract("[:alpha:]+$")) %>%
        filter(pred_type == thresh_type) %>%
        group_by(cell, pred_type) %>%
        summarise(pred_thresh = sum(ifelse(pred_value >= thresh_value,
                                                1, 0))) %>%
        mutate(pred_type = pred_type %>%
                   paste0("pred_thresh_", .)) %>%
        spread(pred_type, pred_thresh)

    list(obs_rich_df, pred_sum_df, pred_thresh_df) %>%
        Reduce(function(x, y) inner_join(x, y, by = c("cell")), .)
}


#' Function to extract deviation of species richness
#'
#' @import tidyr dplyr
#' @export
extract_rich_component = function(rich_df, component) {

    component = enquo(component)

    comp_name = quo_name(component)

    comp_type  = paste0(comp_name, "_type")
    comp_value = paste0(comp_name, "_value")

    rich_df %>%
        ungroup() %>%
        dplyr::select(hab_sat, cell, obs_rich, contains(comp_name)) %>%
        tidyr::gather(!!comp_type, !!comp_value, contains(comp_name)) %>%
        extract(!!comp_type, c(comp_type, "model_type"),
                regex = "([:alnum:]+_[:alnum:]+)_([:alnum:]+)")
}

gather_by = function(df, prefix) {

    name_col  = paste0(prefix, "_type")
    value_col = paste0(prefix, "_value")

    gather(df, !!name_col, !!value_col, contains(prefix))
}

#' Composition accuracy from threshold predictions
#'
#' This function computes the accuracy (specificity and sensitivity) of the
#' prediction of composition of assemblages. It quantifies the amount to which
#' the prediction using thresholds predict well the real presence and the real
#' absence of species at the assemblage level
#'
#' @import tidyr dplyr
#' @export
predict_thresh_composition_accuracy = function(sdm_df) {
    sdm_df %>%
        unnest(obs) %>%
        dplyr::select(hab_sat, cell, presence, contains("thresh"),
                      contains("pred")) %>%
        gather_by("pred") %>%
        gather_by("thresh") %>%
        mutate(pred_type = pred_type %>%
                   stringr::str_extract("[:alpha:]+$"),
               thresh_type = thresh_type %>%
                   stringr::str_extract("[:alpha:]+$")) %>%
        filter(pred_type == thresh_type) %>%
        dplyr::select(-thresh_type) %>%
        mutate(bin_pred = ifelse(pred_value >= thresh_value, 1, 0)) %>%
        group_by(hab_sat, cell, pred_type) %>%
        summarise(n_pp = sum(ifelse(presence == 1 & bin_pred == 1, 1, 0)),
                  n_pa = sum(ifelse(presence == 1 & bin_pred == 0, 1, 0)),
                  n_ap = sum(ifelse(presence == 0 & bin_pred == 1, 1, 0)),
                  n_aa = sum(ifelse(presence == 0 & bin_pred == 0, 1, 0)),
                  assemblage_specificity = n_aa / (n_aa + n_pa),
                  assemblage_sensitivity = n_pp / (n_pp + n_ap))
}

#' Composition accuracy from predicted probabilities
#'
#' This function computes the accuracy (specificity and sensitivity) of the
#' prediction of composition of assemblages. It quantifies the amount to which
#' the predicted probabilities match well the real presence and the real
#' absence of species at the assemblage level. To do so it needs to draw
#' distribution from the probabilities and computes the average specificity and
#' sensitivity at the assemblage level
#'
#' @import tidyr dplyr
#' @export
predict_prob_composition_accuracy = function(sdm_df, n_draws = 10) {

    # Number of draws from predicted presence probabilities
    n_draws = n_draws

    sdm_df %>%
        unnest(obs) %>%
        dplyr::select(hab_sat, cell, presence, contains("pred")) %>%
        gather_by("pred") %>%
        mutate(pred_type = pred_type %>%  # Simplify model name
                   stringr::str_extract("[:alpha:]+$"),
               # Make sure predicted probabilities have values in [0;1]
               pred_value = case_when(pred_value < 0 ~ 0,
                                      pred_value > 1 ~ 1,
                                      TRUE ~ pred_value),
               # Random binomial draw proportional to predicted probability
               pred_draw = purrr::map(pred_value, ~rbinom(n_draws, 1, .x))) %>%
        unnest(pred_draw) %>%
        # Add to each line the draw number
        mutate(draw_number = rep_len(1:n_draws, length.out = nrow(.))) %>%
        # Compute sensitivity and specificity
        group_by(hab_sat, cell, pred_type, draw_number) %>%
        summarise(n_pp = sum(ifelse(presence == 1 & pred_draw == 1, 1, 0)),
                  n_pa = sum(ifelse(presence == 1 & pred_draw == 0, 1, 0)),
                  n_ap = sum(ifelse(presence == 0 & pred_draw == 1, 1, 0)),
                  n_aa = sum(ifelse(presence == 0 & pred_draw == 0, 1, 0)),
                  assemblage_specificity = n_aa / (n_aa + n_pa),
                  assemblage_sensitivity = n_pp / (n_pp + n_ap)) %>%
        summarise_at(vars(contains("assemblage")), mean)
}

get_pred_obs_slope_signif = function(pred_rich) {

    pred_rich %>%
        filter(model_type == "glm") %>%
        nest(-hab_sat, -pred_type) %>%
        mutate(obs_pred_mod = purrr::map(data, ~lm(obs_rich ~ pred_value,
                                                   data = .x)),
               mod_coef     = purrr::map(obs_pred_mod, broom::tidy)) %>%
        unnest(mod_coef) %>%
        filter(term == "pred_value") %>%
        mutate(slope_t_value = (estimate - 1)/std.error,
               pval_slope = 2*pt(abs(slope_t_value), 1998, lower.tail = FALSE))
}
