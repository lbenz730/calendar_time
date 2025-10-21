library(tidyverse)
library(lubridate)
library(arrow)
library(glue)
library(haven)
source('scripts/helpers.R')

### Directory where EHR data is stored
ehr_dir <- '/n/haneuse_ehr_l3/V1.0'
data_dir <- '/n/haneuse_ehr_l3/V1.0/clean_datasets'

### Anti-Lipemic Medications
antilipemic_cases_rx <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/antilipemics_rx.sas7bdat'))
antilipemic_controls_rx <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/antilipemics_rxfill_controls.sas7bdat'))

antilipemic_rx <- 
  bind_rows(
    antilipemic_cases_rx %>% 
      select('subject_id' = durable_studyid,
             ndc,
             'rxdate' = RXDATE, 
             'rxamt' = RXAMT,
             'rxsup' = RXSUP,
             rx_year,
             brand,
             generic),
    
    antilipemic_controls_rx %>% 
      select('subject_id' = control_studyid,
             ndc,
             rxdate,
             rxamt,
             rxsup,
             rx_year,
             'brand' = brand_name,
             'generic' = generic_name)
  )

write_parquet(antilipemic_rx, glue('{data_dir}/antilipemics_rx.parquet'))

### Hypertension Medications
hypertension_cases_rx <- read_sas(glue('{ehr_dir}/_data_for_eric_cases/hypertension_rx.sas7bdat'))
hypertension_controls_rx <- read_sas(glue('{ehr_dir}/_data_for_eric_controls_001/ht_rx_controls.sas7bdat'))


hypertension_rx <- 
  bind_rows(
    hypertension_cases_rx %>% 
      select('subject_id' = durable_studyid,
             'ndc' = NDC,
             rxdate,
             rxamt,
             rxsup,
             rx_year,
             brand,
             generic),
    
    hypertension_controls_rx %>% 
      select('subject_id' = control_studyid,
             'ndc' = NDC,
             rxdate,
             rxamt,
             rxsup,
             rx_year,
             brand,
             generic)
  )

write_parquet(hypertension_rx, glue('{data_dir}/hypertension_rx.parquet'))
