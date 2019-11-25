# Script that contains function for drake setup

#' Full analysis plan to generate the different stesp
#'
#' @import drake
#' @importFrom dplyr bind_rows
#' @export
saturation_workflow = function() {

    # Constants ----------------------------------------------------------------

    # Habitat saturation values to be used in simulations
    hab_sat_values = round(seq(0.1, 1.7, length.out = 8), 1)

    # Generate environmental gradient, species and assemblages
    master_plan = drake_plan(
        art_raster    = create_env_raster(2000),
        art_species   = simulate_species(env_raster = art_raster,
                                         given_seed = 27022018),
        sampled_sites = sample_species_pa(art_species, art_raster, 2000))

    # Modulate species habitat saturations
    hab_sat_alt = drake_plan(
        hab_sat_alternative = target(
            change_alt_hab_sat(art_species, hab_sat_val),
            transform = map(hab_sat_val = !!hab_sat_values)),
        # Fit SDMS -------------------------------------------------------------
        sdm_df_alt = target(
            build_sdm(hab_sat_alternative),
            transform = map(hab_sat_alternative)
        ),
        all_sdm = target(
            dplyr::bind_rows(sdm_df_alt),
            transform = combine(sdm_df_alt)
        ),
        # Estimate species richness --------------------------------------------
        rich_df_alt = target(
            predict_richness(sdm_df_alt),
            transform = map(sdm_df_alt)
        ),
        all_rich = target(
            dplyr::bind_rows(rich_df_alt),
            transform = combine(rich_df_alt)
        ),
        pred_rich = extract_rich_component(all_rich, pred),
        pred_rich_rmse = pred_rich %>%
            group_by(hab_sat, pred_type, model_type) %>%
            summarise(rmse = sqrt(mean((pred_value - obs_rich)^2)),
                      rich_cor = cor(obs_rich, pred_value,
                                     method = "spearman"),
                      variance = mean((pred_value - mean(pred_value))^2),
                      bias     = mean(pred_value - obs_rich)) %>%
            tidyr::gather("measure_type", "measure", rmse, rich_cor,
                          variance, bias) %>%
            mutate(measure_type = factor(measure_type,
                                         levels = c("rich_cor", "rmse", "bias",
                                                    "variance"))),
        # Test significance of slope against one -------------------------------
        rich_slopes = all_rich %>%
            tidyr::nest(richness_df = c(cell, obs_rich,
                                        starts_with("pred"))) %>%
            mutate(
                richness_lm = purrr::map(
                    richness_df, ~lm(obs_rich ~ pred_sum_glm, data = .x)),
                richness_coef = purrr::map(richness_lm, broom::tidy)
            ) %>%
            unnest(richness_coef) %>%
            filter(term == "pred_sum_glm") %>%
            # This test to what extent the slope of probability-based richness
            # against observed richness is different than one
            # see https://stats.stackexchange.com/a/137122/223164
            mutate(p_value_against_slope_one = pt(q = (estimate - 1)/std.error,
                                                  df = 2000 - 2,
                                                  lower.tail = FALSE))
    )

    # Draw figures -------------------------------------------------------------
    figure_plan = drake_plan(
        fig_theo_real_one_sp = compare_theor_obtained_responses(all_sdm,
                                                                art_species,
                                                                "virt_9"),
        fig_rich_pred_obs    = plot_rich_pred_obs(pred_rich),
        fig_rich_rmse_cor    = plot_rich_rmse_bias_var(pred_rich_rmse)
    )

    # Manuscript ---------------------------------------------------------------

    manuscript_plan = drake_plan(
        manuscript = target({
            file_out("manuscript/estimating_richness_sdm.docx")
            file_out("manuscript/estimating_richness_sdm.pdf")
            file_in("manuscript/saturation_ms.bib")
            rmarkdown::render(
                knitr_in("manuscript/estimating_richness_sdm.Rmd"),
                "all",
                encoding = "UTF-8")
        })
    )


    # Final plan ---------------------------------------------------------------
    bind_plans(master_plan,
               hab_sat_alt,
               figure_plan,
               manuscript_plan)
}
