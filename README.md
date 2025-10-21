# Calendar Time-Varying Effects

Benz, L., Mukherjee, R., Wang, R., Arterburn, D., Fischer, H., Lee, C., Shortreed, S.M., Haneuse, S., and Levis, A.W. "A Statistical Framework for Understanding Causal Effects that Vary by Treatment Initiation Time in EHR-based Studies" _Under Review_ ([Pre-Print]( LINK TO COME ))


## R Scripts (`scripts/`)

* __helpers.R__: File of helper functions 

### Data (`scripts/data`)

* __build_weight_trials.R__: Build target trials for TTE dataset w/ relative % weight change as the outcome
* __process_medications.R__: Script to process some hypertension + antilipemic medications for use in the weight loss trials dataset.

### Analysis (`scripts/analysis`)

* __EIF_helpers.R__: Helper function for EIF for $\beta$ and $\chi_m$
* __EIF_weight_analysis.R__: EIF applied to weight outcomes (estimation of $\widehat\chi_m$)
* __EIF_weight_projections.R__: Projection step (estimation of $\widehat\beta$)
* __standardization_EIF_weight.R__: Standardization analysis to compute $\widehat\chi_{j,m}$
* __hypothesis_testing.R__: Run hypothesis testing for $\theta$ statistic 
* __plot_EIF_weight_results.R__: Plot results from data application (Figure 5, Table 2)
* __plot_confounder_distributions.R__: Look at patterns of covariates over time (Figures 1, S1, S2)


### Simulations (`scripts/simulations`)
* __inform_sims.R__: Get realistic range of coefficients for generating simulated data
* __generate_data.R__: Functions to generate simulated data sets and compute true effects for simulation scenarios
* __specify_inputs.R__: Specify and visualize structure of treatment effects across simulation scenarios (Figure 3)
* __compute_true_sigma_gamma_ratios.R__: Script to compute and visualize the true values of $\frac{\sigma_m^2}{\sigma_m^2 + \gamma_m^2}$ for simulation data generating process.
* __run_simulation_pipeline.R__: Wrapper to run simulation iterations
* __analyze_sim_results.R__: Analysis of simulation results (Figure 4)
* __latex_table.R__: Generate summary table of sim parameters (Table S1)

## Simulation Inputs (`sim_inputs/`)
.rds files containing the simulation parameters for each of the 18 simulation scenarios

## Figures (`figures/`)
Figures saved out from various analyses

## Shell (`shell/`)
.sh files for transferring files to and from the cluster.

## Jobs (`jobs/`)
.sh files for batch jobs on the cluster 

* __build_weight_trials.sh__: Prep data for analyses w/ weight loss trials
* __EIF_weight.sh__: Wrapper to run EIF weight analysis
* __EIF_transport.sh__: Job array for EIF weight analysis transported to various populations
* __run_sims.sh__: Wrapper to run job arrays for simulations 
    * __sims_n100.sh__: $n = 100$ subjects
    * __sims_n500.sh__: $n = 500$ subjects
    * __sims_n1000.sh__: $n = 1,000$ subjects
    * __sims_n5000.sh__: $n = 5,000$ subjects
    * __sims_n10000.sh__: $n = 10,000$ subjects
    * __sims_n25000.sh__: $n = 25,000$ subjects