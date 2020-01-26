# Script that contains functions to make figures


# Helper functions -------------------------------------------------------------

#' @importFrom ggplot2 as_labeller
get_beta_hab_labs = function(given_df) {
    unique(given_df$hab_sat) %>%
        purrr::set_names(., .) %>%
        purrr::map_chr(~paste0("Habitat Saturation = ", .x * 100, "%")) %>%
        as_labeller()
}

# Figure 2: Predictions on a single species ------------------------------------
#' @importFrom virtualspecies quadraticFun
compare_theor_obtained_responses = function(full_sdm_df, virtual_sp,
                                            given_sp = "virt_9") {

    # Get given virtual species
    given_virtual = virtual_sp[[given_sp]]

    # Get parameters and function
    virtual_params = given_virtual[["details"]][["parameters"]][["unidim"]]

    pa_conv = given_virtual[["PA.conversion"]]

    # Compute theoretical suitability from parameters estimated from virtual
    # species' characteristics
    theoretical_curve = data_frame(
        cell = seq(virtual_params$min, virtual_params$max)) %>%
        mutate(pred_glm = do.call(virtual_params$fun,
                                  list(list(x = cell),
                                       as.list(virtual_params$args)) %>%
                                      purrr::flatten())) %>%
        mutate(pred_glm = scales::rescale(pred_glm, to = c(0, 1)),
               env_suitab = virtualspecies::logisticFun(pred_glm,
                                                        pa_conv[["alpha"]] %>%
                                                            as.numeric(),
                                                        pa_conv[["beta"]] %>%
                                                            as.numeric()) %>%
                   scales::rescale(to = c(0, 1)),
               total_prob = (pred_glm * env_suitab) %>%
                   scales::rescale(to = c(0, 1)))

    # Get range of when thresholds predict presences or absences
    pres_abs_segments = full_sdm_df %>%
        filter(species == given_sp) %>%
        unnest(obs) %>%
        mutate(pred_thresh = ifelse(pred_glm >= thresh_glm, 1, 0)) %>%
        group_by(hab_sat, pred_thresh) %>%
        summarise(min = first(cell), max = last(cell)) %>%
        mutate(bla = lst(range(min)),
               two = lst(range(max))) %>%
        unnest(bla, two) %>%
        ungroup() %>%
        rowwise() %>%
        mutate(real_range = case_when(min == bla & two < max ~ c(two, max),
                                      two == max & min < bla ~ c(min, bla),
                                      two == max & bla < min ~ c(min, max)) %>%
                   lst()) %>%
        filter(!is.na(real_range[[1]])) %>%
        dplyr::select(-min, -max, -bla, -two) %>%
        mutate(min_val = real_range[[1]],
               max_val = real_range[[2]],
               # This defines the relative height of threshold predictions
               # that are going to be above and beneath the rest of the plot
               height = case_when(pred_thresh == 1 ~ pred_thresh +
                                      hab_sat * 0.09,
                                  pred_thresh == 0 ~ - hab_sat * 0.09))

    # Plot
    full_sdm_df %>%
        filter(species == given_sp) %>%
        unnest(obs) %>%
        ggplot(aes(cell, pred_glm, color = as.factor(hab_sat))) +
        geom_line(size = 1.2) +
        geom_line(data = theoretical_curve, aes(y = total_prob, color = NULL,
                                                linetype = "theory"),
                  show.legend = c(color = FALSE, linetype = TRUE), size = 0.6) +
        geom_segment(data = pres_abs_segments, size = 1,
                     aes(x = min_val, xend = max_val, y = pred_thresh,
                         yend = pred_thresh), alpha = 1/7) +
        scale_color_viridis_d(direction = -1, guide = guide_legend(order = 1),
                              labels = function(x) x %>%
                                  as.character() %>%
                                  as.numeric() %>%
                                  scales::percent()) +
        scale_linetype_manual(values = c(theory = 3),
                              labels = c(theory = "Expected\nRelationship"),
                              guide = guide_legend(order = 2)) +
        labs(x        = "Environment",
             y        = "Predicted Presence Probability",
             color    = "Habitat Saturation",
             linetype = NULL) +
        theme_bw()
}

# Figure 3: Observed vs. Predicted Richness ------------------------------------
plot_rich_pred_obs = function(hab_sat_df) {

    facet_labs = get_beta_hab_labs(hab_sat_df)

    par(lheight = 0.5)

    ggplot(hab_sat_df, aes(pred_value, obs_rich, color = pred_type)) +
        geom_abline(slope = 1, intercept = 0, linetype = 2) +
        geom_point(alpha = 1/20, size = 1/2) +
        stat_smooth(se = FALSE) +
        ggpubr::stat_cor(aes(label = tolower(..r.label..)),
                         method = "spearman") +
        facet_wrap(~hab_sat, labeller = facet_labs, ncol = 4) +
        scale_color_brewer(labels = c(pred_sum = "Probabilities",
                                      pred_thresh = "Thresholds"),
                           type = "qual") +
        labs(x        = "Predicted Species Richness",
             y        = "Observed Species Richness",
             color    = "Methods") +
        theme_bw() +
        theme(legend.position  = "top",
              panel.grid       = element_blank(),
              strip.background = element_blank(),
              aspect.ratio     = 1)
}

# Figure 4: Models RMSE, Bias and Variance -------------------------------------

#' @import ggplot2
plot_rich_rmse_bias_var = function(pred_rich_rmse) {
    pred_rich_rmse %>%
        filter(model_type == "glm", measure_type != "rich_cor") %>%
        ggplot(aes(as.factor(hab_sat), measure, color = pred_type)) +
        geom_point(size = 2) +
        geom_line(aes(group = pred_type), size = 1) +
        facet_wrap(
            vars(measure_type), scales = "free_y",
            labeller = labeller(
                measure_type = c(
                    rmse = "Root Mean Squared Error",
                    variance = "Variance",
                    bias     = "Bias")), ncol = 3) +
        labs(x = "Habitat Saturation", color = "Method", y = "") +
        scale_color_brewer(labels = c(pred_sum = "Probabilities",
                                      pred_thresh = "Thresholds"),
                           type = "qual") +
        scale_x_discrete(
            labels = function(x) scales::percent(as.numeric(x))
        ) +
        theme_bw(15) +
        theme(aspect.ratio = 1,
              strip.background = element_blank(),
              axis.text = element_text(color = "black", size = rel(0.8)),
              axis.text.x = element_text(angle = 80, hjust = 1),
              legend.position = "top")
}
