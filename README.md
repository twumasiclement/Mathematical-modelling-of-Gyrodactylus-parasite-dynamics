# Mathematical modelling of *Gyrodactylus* parasite dynamics

Understanding host-parasite systems are challenging if biologists only adopt experimental approaches, whereas mathematical models can help uncover other in-depth knowledge about host infection dynamics. The gyrodactylid-fish system, like other host-parasite systems, is widely used to investigate ecological, evolutionary, and epidemiological problems. They are particularly amenable to experimental manipulation, and there is an individual-based model (IBM) to reproduce their within-host infection dynamics. However, spatial information on species-specific microhabitat preference and other relevant biological information for different *Gyrodactylus* strains across different fish populations are yet to be fully captured in the existing IBM. Therefore, this study aims to develop a novel individual-based stochastic simulation model via a hybrid $\tau$-leaping to include these data to add to our understanding of the gyrodactylid-fish system's complexity. The infection dynamics of three gyrodactylid strains are compared across three host populations. A modified sequential-type approximate Bayesian computation (ABC) is developed for fitting this sophisticated model based on empirical data and an auxiliary stochastic model. A penalised local-linear regression for ABC post-processing analysis is proposed. A linear birth-death process with catastrophic extinction (B-D-C process) is considered the auxiliary model to refine the modified ABC's summary statistics, with other theoretical justifications and parameter estimation techniques of the B-D-C process provided. The B-D-C process simulation using $\tau$-leaping also provided additional insights on accelerating the complex simulation model. The mathematical models can be adapted for other host-parasite systems, and the modified ABC methodologies can aid in efficiently calibrating other multi-parameter models with a high-dimensional set of correlating or independent summary statistics.

**NB: All R codes developed for statistical analyses and mathematical modelling (with Jupyter notebook HTML, R source files and other saved results) as well as the empirical data (for this study) are attached. The zipped file `B-D-C_Model_ R codes.zip` contains all R codes and other results concerning the B-D-C auxiliary model. The zipped file  `SimulationModel_ABC_Rcodes.zip` contains a detailed R codes of the novel individual-based stochastic simulation model, model fitting based on the ABC, the empirical parasite data is named as `Parasite_Data.csv`, the data on the body area of fish is named `Area_Fish_bodyParts.csv`, and `w0.csv` is the initial summary statistics weights for the ABC analysis (based on initial simulations at time ABC $t=0$). `Figures_paper.zip` is a zip file of all Figures in the papers.**
 
