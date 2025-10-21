library(tidyverse)
source('scripts/simulations/generate_data.R') 

### Baseline Set of Parameters 
base_params <- 
  list('n_subjects' = 10000,
       'n_trials' = 36,
       'treatment_model' = 
         list('(Intercept)' = -2.48522,
              'baseline_bmi' = 0.06729,
              'gender' = -1.00825,
              'race' = 0.40550,
              'site[NC]' = -0.38618,
              'site[SC]' = 0.34849,
              'baseline_age' = -0.04932,
              't2dm' = -0.03380,
              'insulin' = -0.01976,
              'hypertension' = 0.70059,
              'hypertension_rx' = -0.34257,
              'dyslipidemia' = 0.42442,
              'antilipemic_rx' = -0.21276,
              'smoking_status[former]' = 2.16270,
              'smoking_status[never]' = 1.65722),
       
       'sigma_rss' = 0.025
  )

### Covariate Shift Options
covariate_shifts <- 
  list('no_covariate_shift' = list('shift_fx' = NULL, 
                                   'description' = 'No Covariate Shift'),
       'linear' = list('shift_fx' = 
                         list('baseline_age' = list('X' = cbind(1, matrix(1:36)), 
                                                    'beta' = c(1, -1/12, 1/12), 
                                                    'random_walk' = F, 
                                                    'probability' = F),
                              'baseline_bmi' = list('X' = cbind(1, matrix(1:36)),
                                                    'beta' = c(1, -1/12, 1/12),
                                                    'random_walk' = T,
                                                    'probability' = F, 
                                                    'walk_sd' = 1),
                              't2dm' = list('X' = cbind(1, matrix(1:36)),
                                            'beta' = c(0, (0.2 - 0.4/35), 0.4/35),
                                            'random_walk' = F,
                                            'probability' = T)),
                       'description' = 'Linear Covariate Shift'),
       'flexible' = list('shift_fx' = 
                           list('baseline_age' = list('X' = cbind(1, matrix(1:36)), 
                                                      'beta' = c(1, -1/12, 1/12), 
                                                      'random_walk' = F, 
                                                      'probability' = F),
                                'baseline_bmi' = list('X' = ns(1:36, df = 3),
                                                      'beta' = c(1, 1, -3, -1),
                                                      'random_walk' = T,
                                                      'probability' = F, 
                                                      'walk_sd' = 1),
                                't2dm' = list('X' = cbind(1, ns(1:36, df = 3)),
                                              'beta' = c(0, 0.2, 0.02, 0.6, -0.08),
                                              'random_walk' = F,
                                              'probability' = T)),
                         'description' = 'Flexible Covariate Shift'))


outcome_models <- 
  list('constant' = list('coeff' = 
                           list('(Intercept)' = 0.09787,
                                'surgery' = -0.21156,
                                'baseline_bmi' = -0.00189,
                                'gender' = 0.00471,
                                'race' = -0.00306,
                                'site[NC]' = -0.00236,
                                'site[SC]' = -0.00789,
                                'baseline_age' = -0.00061,
                                't2dm' = -0.01167,
                                'insulin' = 0.01274,
                                'hypertension' = 0.00132,
                                'hypertension_rx' = -0.00008,
                                'dyslipidemia' = -0.00088,
                                'antilipemic_rx' = 0.00036,
                                'smoking_status[former]' = -0.00258,
                                'smoking_status[never]' = 0.00046),
                         'description' = 'Constant'),
       
       'effect_mod' = list('coeff' = 
                             list('(Intercept)' = 0.09787,
                                  'surgery' = 0.66644,
                                  'baseline_age:surgery' = 0.002,
                                  't2dm:surgery' = 0.1,
                                  'baseline_bmi:surgery' = -0.025,
                                  'baseline_bmi' = -0.00189,
                                  'gender' = 0.00471,
                                  'race' = -0.00306,
                                  'site[NC]' = -0.00236,
                                  'site[SC]' = -0.00789,
                                  'baseline_age' = -0.00061,
                                  't2dm' = -0.01167,
                                  'insulin' = 0.01274,
                                  'hypertension' = 0.00132,
                                  'hypertension_rx' = -0.00008,
                                  'dyslipidemia' = -0.00088,
                                  'antilipemic_rx' = 0.00036,
                                  'smoking_status[former]' = -0.00258,
                                  'smoking_status[never]' = 0.00046),
                           'description' = 'Effect Modification'),
       
       'linear' = list('coeff' = 
                         list('(Intercept)' = 0.09787,
                              'surgery' = -0.21156,
                              'trial_id:surgery' = 0.001,
                              'baseline_bmi' = -0.00189,
                              'gender' = 0.00471,
                              'race' = -0.00306,
                              'site[NC]' = -0.00236,
                              'site[SC]' = -0.00789,
                              'baseline_age' = -0.00061,
                              't2dm' = -0.01167,
                              'insulin' = 0.01274,
                              'hypertension' = 0.00132,
                              'hypertension_rx' = -0.00008,
                              'dyslipidemia' = -0.00088,
                              'antilipemic_rx' = 0.00036,
                              'smoking_status[former]' = -0.00258,
                              'smoking_status[never]' = 0.00046),
                       'description' = 'Linear'),
       
       'spline' = list('coeff' = 
                         list('(Intercept)' = 0.09787,
                              'surgery' = -0.21156,
                              'surgery:spline_1' = 0.01,
                              'surgery:spline_2' = 0.06,
                              'surgery:spline_3' = 0.01,
                              'baseline_bmi' = -0.00189,
                              'gender' = 0.00471,
                              'race' = -0.00306,
                              'site[NC]' = -0.00236,
                              'site[SC]' = -0.00789,
                              'baseline_age' = -0.00061,
                              't2dm' = -0.01167,
                              'insulin' = 0.01274,
                              'hypertension' = 0.00132,
                              'hypertension_rx' = -0.00008,
                              'dyslipidemia' = -0.00088,
                              'antilipemic_rx' = 0.00036,
                              'smoking_status[former]' = -0.00258,
                              'smoking_status[never]' = 0.00046),
                       'description' = 'Spline'),
       
       'linear_effect_mod' = list('coeff' = 
                                    list('(Intercept)' = 0.09787,
                                         'surgery' = 0.66644,
                                         'trial_id:surgery' = 0.001,
                                         'baseline_age:surgery' = 0.002,
                                         't2dm:surgery' = 0.1,
                                         'baseline_bmi:surgery' = -0.025,
                                         'baseline_bmi' = -0.00189,
                                         'gender' = 0.00471,
                                         'race' = -0.00306,
                                         'site[NC]' = -0.00236,
                                         'site[SC]' = -0.00789,
                                         'baseline_age' = -0.00061,
                                         't2dm' = -0.01167,
                                         'insulin' = 0.01274,
                                         'hypertension' = 0.00132,
                                         'hypertension_rx' = -0.00008,
                                         'dyslipidemia' = -0.00088,
                                         'antilipemic_rx' = 0.00036,
                                         'smoking_status[former]' = -0.00258,
                                         'smoking_status[never]' = 0.00046),
                                  'description' = 'Linear + Effect Modification'),
       
       'spline_effect_mod' = list('coeff' = 
                                    list('(Intercept)' = 0.09787,
                                         'surgery' = 0.66644,
                                         'surgery:spline_1' = 0.01,
                                         'surgery:spline_2' = 0.06,
                                         'surgery:spline_3' = 0.01,
                                         'baseline_age:surgery' = 0.002,
                                         't2dm:surgery' = 0.1,
                                         'baseline_bmi:surgery' = -0.025,
                                         'baseline_bmi' = -0.00189,
                                         'gender' = 0.00471,
                                         'race' = -0.00306,
                                         'site[NC]' = -0.00236,
                                         'site[SC]' = -0.00789,
                                         'baseline_age' = -0.00061,
                                         't2dm' = -0.01167,
                                         'insulin' = 0.01274,
                                         'hypertension' = 0.00132,
                                         'hypertension_rx' = -0.00008,
                                         'dyslipidemia' = -0.00088,
                                         'antilipemic_rx' = 0.00036,
                                         'smoking_status[former]' = -0.00258,
                                         'smoking_status[never]' = 0.00046),
                                  'description' = 'Spline + Effect Modification')  
       
  )


### All Input Combinations
df_inputs <- 
  crossing('outcome_model' = names(outcome_models),
           'covariate_shift' = names(covariate_shifts))

write_csv(df_inputs, 'sim_inputs/input_summary.csv')

### Save All Inputs
for(i in 1:nrow(df_inputs)) {
  params <- base_params
  params$covariate_shift_fx <- covariate_shifts[[ df_inputs$covariate_shift[i] ]]$shift_fx
  params$covariate_shift <- !is.null(params$covariate_shift_fx)
  params$outcome_model <- outcome_models[[ df_inputs$outcome_model[i] ]]$coeff
  params$shift <- covariate_shifts[[ df_inputs$covariate_shift[i] ]]$description
  params$trt_change <- outcome_models[[ df_inputs$outcome_model[i] ]]$description
  write_rds(params, glue('sim_inputs/params_{i}.rds'))
}


### Plot Treatment Effects over Calendar Time
set.seed(562857)
df_effects <- NULL
for(i in 1:nrow(df_inputs)) {
  cat('Combo', i, 'of', nrow(df_inputs), '\n')
  params <- read_rds(glue('sim_inputs/params_{i}.rds'))
  df_truth <- 
    compute_truth(params) %>% 
    mutate('parameter_id' = i) %>% 
    select(parameter_id, everything())
  
  df_effects <- 
    df_effects %>% 
    bind_rows(df_truth)
}
write_csv(df_effects, 'sim_inputs/true_effects.csv') 


df_effects <- read_csv('sim_inputs/true_effects.csv')
ggplot(df_effects, aes(x = trial_id, y = true_ate)) + 
  facet_wrap(~fct_rev(shift)) + 
  geom_line(aes(color = trt_change)) + 
  scale_y_continuous(labels = scales::percent) + 
  labs(
    #title = 'Simulation Scenario Overview', 
    # subtitle = 'True Calendar Time-Varying Treatment Effect Values',
    x = 'Trial Index (m)',
    y = expression(paste('True ', chi[m](P))),
    color = 'Underlying Effect Structure') + 
  theme(legend.title = element_text(size = 16))
ggsave('figures/true_tv_effects.pdf', height = 4.5, width = 16/1.5)
