library(tidyverse)
library(arrow)
library(glue)
library(table1)

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Custom table 1 functions to get , after the thousands
render.continuous <- function(x, ...) {
  with(stats.default(x, ...), c("",
                                "Mean (SD)"         = sprintf("%s (%s)",
                                                              signif_pad(MEAN,   3, big.mark=","),
                                                              signif_pad(SD,     3, big.mark=","))))
}

render.categorical <- function(x, ...) {
  c("", sapply(stats.apply.rounding(stats.default(x)), function(y) with(y,
                                                                        sprintf("%s (%s%%)", prettyNum(FREQ, big.mark=","), PCT))))
}

render.strat <- function (label, n, ...) {
  sprintf("<span class='stratlabel'>%s<br><span class='stratn'>(N=%s)</span></span>", 
          label, prettyNum(n, big.mark=","))
}

render.strat.subj_trials <- function (label, n, ...) {
  x <- print(unlist(str_split(names(label), '\\.')))
  n_subj <- rep(NA, length(n)) 
  for(i in 1:length(n)) {
    n_subj[i] <- n_distinct(df_trials$subject_id[df_trials$surgery == x[i] ])
  }
  
  sprintf("<span class='stratlabel'>%s<br><span class='stratn'>(N Subject Trials=%s<br>N Subjects=%s)</span></span>", 
          label, prettyNum(n, big.mark=","),  prettyNum(n_subj, big.mark=","))
}



### Load in Data
df_trials <- read_parquet(glue('{data_dir}/tv_effects/weight_trials/trials_combined.parquet')) 

df_trials <- 
  df_trials %>% 
  filter(eligible) %>% 
  mutate('bs_type' = ifelse(surgery == 0, 'CONTROL', bs_type)) %>% 
  mutate('surgery' = factor(ifelse(surgery == 1, 'Bariatric Surgery', 'No Bariatric Surgery')),
         'hypertension' = fct_rev(factor(ifelse(hypertension == 1, 'Yes', 'No'))),
         'dyslipidemia' = fct_rev(factor(ifelse(dyslipidemia == 1, 'Yes', 'No'))),
         'hypertension_rx' = fct_rev(factor(ifelse(hypertension_rx == 1, 'Yes', 'No'))),
         'antilipemic_rx' = fct_rev(factor(ifelse(antilipemic_rx == 1, 'Yes', 'No'))),
         'insulin' = fct_rev(factor(ifelse(insulin == 1, 'Yes', 'No'))),
         't2dm' =  fct_rev(factor(ifelse(t2dm == 1, 'Yes', 'No')))
         )



label(df_trials$baseline_age) <- 'Baseline Age'
label(df_trials$gender) <- 'Sex'
label(df_trials$race) <- 'Race'
label(df_trials$baseline_bmi) <- 'Baseline BMI'
label(df_trials$hypertension) <- 'Hypertension'
label(df_trials$dyslipidemia) <- 'Dyslipidemia'
label(df_trials$hypertension_rx) <- 'Antihypertensive Medication Usage'
label(df_trials$antilipemic_rx) <- 'Antilipemic Medication Usage'
label(df_trials$t2dm) <- 'Type II Diabetes Mellitus'
label(df_trials$insulin) <- 'Insulin Use'
label(df_trials$smoking_status) <- 'Self-Reported Smoking Status'


tbl_1 <- 
  table1(~baseline_age + gender + baseline_bmi + race + site + 
           hypertension + hypertension_rx + dyslipidemia + antilipemic_rx + 
           t2dm + insulin + smoking_status + 
           bs_type | surgery,
         data = df_trials, 
         render.continuous = render.continuous,
         render.strat = render.strat.subj_trials,
         render.categorical = render.categorical,
         overall = F)
tbl_1
write_lines(tbl_1, 'figures/population_summary_table.html')
