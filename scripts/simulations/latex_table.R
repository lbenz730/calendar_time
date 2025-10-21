library(tidyverse)
library(arrow)
library(glue)
library(xtable)
library(kableExtra)

source('scripts/helpers.R') 

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

outcome_models <- 
  list('constant' = list('coeff' = 
                           list('(Intercept)' = 0.09787,
                                'treatment' = -0.21156,
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
                                  'treatment' = 0.66644,
                                  'baseline_age:treatment' = 0.002,
                                  't2dm:treatment' = 0.1,
                                  'baseline_bmi:treatment' = -0.025,
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
                              'treatment' = -0.21156,
                              'trial_id:treatment' = 0.001,
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
                              'treatment' = -0.21156,
                              'treatment:spline_1' = 0.01,
                              'treatment:spline_2' = 0.06,
                              'treatment:spline_3' = 0.01,
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
                                         'treatment' = 0.66644,
                                         'trial_id:treatment' = 0.001,
                                         'baseline_age:treatment' = 0.002,
                                         't2dm:treatment' = 0.1,
                                         'baseline_bmi:treatment' = -0.025,
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
                                         'treatment' = 0.66644,
                                         'treatment:spline_1' = 0.01,
                                         'treatment:spline_2' = 0.06,
                                         'treatment:spline_3' = 0.01,
                                         'baseline_age:treatment' = 0.002,
                                         't2dm:treatment' = 0.1,
                                         'baseline_bmi:treatment' = -0.025,
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






### Table of Simulation Parameters
df_coeff <- 
  bind_rows(
    tibble('model' = 'Treatment',
           'model_num' = 1,
           'term_num' = 1,
           'coeff_greek' = '$\\bm \\pi_m$',
           'term' = names(base_params$treatment_model),
           'coeff_value' = unlist(base_params$treatment_model)),
    
    tibble('model' = 'Constant',
           'model_num' = 2,
           'term_num' = 1,
           'coeff_greek' = '$\\bm \\mu_m$',
           'term' = names(outcome_models$constant$coeff),
           'coeff_value' = unlist(outcome_models$constant$coeff)),
    
    tibble('model' = 'Linear',
           'model_num' = 2,
           'term_num' = 2,
           'coeff_greek' = '$\\bm \\mu_m$',
           'term' = names(outcome_models$linear$coeff),
           'coeff_value' = unlist(outcome_models$linear$coeff)),
    
    tibble('model' = 'Spline',
           'model_num' = 2,
           'term_num' = 3,
           'coeff_greek' = '$\\bm \\mu_m$',
           'term' = names(outcome_models$spline$coeff),
           'coeff_value' = unlist(outcome_models$spline$coeff)),
    
    tibble('model' = 'Effect Modification',
           'model_num' = 2,
           'term_num' = 4,
           'coeff_greek' = '$\\bm \\mu_m$',
           'term' = names(outcome_models$effect_mod$coeff),
           'coeff_value' = unlist(outcome_models$effect_mod$coeff)),
    
    tibble('model' = 'Linear + Effect Modification',
           'model_num' = 2,
           'term_num' = 5,
           'coeff_greek' = '$\\bm \\mu_m$',
           'term' = names(outcome_models$linear_effect_mod$coeff),
           'coeff_value' = unlist(outcome_models$linear_effect_mod$coeff)),
    
    tibble('model' = 'Spline + Effect Modification',
           'model_num' = 2,
           'term_num' = 6,
           'coeff_greek' = '$\\bm \\mu_m$',
           'term' = names(outcome_models$spline_effect_mod$coeff),
           'coeff_value' = unlist(outcome_models$spline_effect_mod$coeff))
  )

sci_notation <- function(x) {
  case_when(is.na(x) ~ '---',
            x == 0 ~ '0',
            abs(x) > 0.1 ~ sprintf('%0.2f', x),
            T ~ sanitize.numbers(format(x, 
                                        digits = 2,
                                        scientific = T),
                                 type = 'latex',
                                 math.style.exponents = T))
  
}

df_tbl <- 
  df_coeff %>% 
  mutate('term' = gsub('\\_', '\\\\_', term)) %>% 
  mutate('term' = gsub('\\^', '\\\\^', term)) %>% 
  mutate('term' = paste0('\\texttt{', term, '}')) %>% 
  select(-coeff_greek, -model) %>% 
  pivot_wider(names_from = c(model_num, term_num),
              values_from = coeff_value) %>% 
  mutate_at(vars(contains('_')), ~map_chr(.x, sci_notation))


latex <-
  df_tbl %>% 
  kbl(align = 'c', 
      digits = 3 , 
      escape = F, 
      format = 'latex') 


