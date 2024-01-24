# Paper Title: Mathematical modelling of parasite dynamics: A stochastic simulation-based approach and parameter estimation via modified sequential-type approximate Bayesian computation

The development of mathematical models for studying newly emerging and reemerging infectious diseases has gained momentum due to global events. The gyrodactylid-fish system, like many host-parasite systems, serves as a valuable resource for ecological, evolutionary, and epidemiological investigations owing to its ease of experimental manipulation and long-term monitoring. Although this system has an existing individual-based model (IBM), it falls short in capturing information about species-specific microhabitat preferences and other biological details for different *Gyrodactylus* strains across diverse fish populations. This current study introduces a new individual-based stochastic simulation model that uses a hybrid tau-leaping algorithm to incorporate this essential data, enhancing our understanding of the complexity of the gyrodactylid-fish system. We compare the infection dynamics of three gyrodactylid strains across three host populations. A modified sequential-type approximate Bayesian computation (ABC) method, based on sequential Monte Carlo (SMC) and sequential importance sampling (SIS), is developed. Additionally, we establish two penalised local-linear regression methods (based on L1 and L2 regularisations) for ABC post-processing analysis to fit our model using existing empirical data. With the support of experimental data and the fitted mathematical model, we address open biological questions for the first time and propose directions for future studies on the gyrodactylid-fish system. The adaptability of the mathematical model extends beyond the gyrodactylid-fish system to other host-parasite systems. Furthermore, the modified ABC methodologies provide efficient calibration for other multi-parameter models characterised by a large set of correlated or independent summary statistics.

**NB: All `R source codes` developed for statistical analyses and mathematical modelling, as well as the empirical data (for this study), are attached. The empirical parasite data is named `Parasite_Data.csv`, the data on the body area of fish named `Area_Fish_bodyParts.csv`, and `w0.csv` denoting the initial summary statistics weights for the ABC fitting (based on initial simulations at ABC time $t=0$).**

# Below are the labels of the main R codes/scripts of our individual-based stochastic simulation model:

1. `Simulation-single-fish-script.r`: Codes to simulate data for a single fish (given specific covariate information)
2. `Simulation-Group-fish-script.r`: Codes to simulate data for a group of fish (in a fashion similar to the observed data)
3. `Update-exactSSA-script.r`: Codes for updating exact SSA (i.e., updating simulation events across the 4 host's body regions: Tail, Lower region, Upper region, & Head)
4. `Update-tauleaping-script.r`: Codes for updating tau-leaping (i.e., updating simulation events across the 4 host's body regions: Tail, Lower region, Upper region, & Head)
5. `Descriptors-Data-script.r`: Codes for extracting experimental descriptors (given the high-dimensional observed data)
   


# Below are the labels of the main R codes/scripts of our modified sequential Monte Carlo ABC sampler (dubbed as Weighted-iterative ABC): 
1. `ABC-Importance-Sampling-Improved-script.R`: Codes for the ABC importance sampling
2. `Weighted-iterative-ABC-script.R`: Codes for the main weighted-iterative ABC or the modified ABC-SMC sampler
3. `Project-Parasite-script.r`: Codes for the function used to project parasite numbers after host mortality
4. `Post-Lasso-reg-adj-L1-script.R`: Codes for the lasso-adjusted ABC posterior correction (L1 regularisation)
5. `Post-Ridge-reg-adj-L2-script.R`: Codes for the ridge-adjusted ABC posterior correction (L2 regularisation)
6. `Sim_ABC_SMC_fit_Improved.R`: Codes for fitting the stochastic model using our modified ABC-SMC sampler
7. `Weighted-distance-script.R`: Codes for computed the weighted sum-of-squares distance
8. `combine-summary-stats-script.R`: Codes for combining the ABC summary statistics into a 2-dimensional array
9. `summary-stats-ABC-script.r`: Codes for computing the ABC summary statistics given a data
10. `MLE_catastrophe-script.r`: Codes for estimating maximum likelihood estimate of the rate of catastrophe (for the B-D-C process)
11. `GMM-1st2nd-Steps-script.r`: Codes for implementing 2-stage generalised methods of moment estimation
12. `BDC-GW-GMM-estimator-script.r`: Codes for calculating the B-D-C parameter estimates
13. `Gyro_ABC_SensitivityAnalysis.r`: Codes for ABC sensitivity analysis and model identifiability explorations (at predefined model parameter values)
14. `GyroModel_Updated ABC fitting_BOMB.r`: Codes for analysis of the ABC outputs and performing all relevant statistical analysis (given the observed empirical data)

    
