# Optimisation-Algorithm-1.1


Inputs: 
The thermodynamic model being optimized is formatted as a table. The results for each experiment should be written as separate files with errors, with an additional file for bulk composition data. 
Data Processing: 
Phase proportions are calculated using a mass balance equation with the inputted data. Corrections for mass balance are non-unique. The strategy employed here takes into account the size of the errors for each component and proportionally adjusts the amount of each oxide component for each phase. The mass balance equation may produce results where some phases are calculated to have no abundance. This may occur when two or more phases have very similar compositions, in which case measured phase abundances should be used instead. 
Endmember corrections remove components not present in the thermodynamic database such as CaO and Na2O in the olivine polymorphs and wustite. Charge balance is also checked to remove Fe3+ where appropriate. 
Parameter Selection: 
Only parameters with errors can be selected for optimization. For general use, parameters from the inputted thermodynamic database are used as the initial guess value. If parallel computing is possible then a random set of initial values can be used instead. This will increase the probability that the algorithm will converge on the global minimum.
Gibbs Free Energy Minimization: 
The algorithm utilizes the system function which allows external programs to be run in MATLAB. This allows users the freedom to choose alternative free energy minimization software if desired. Currently, Perple_X version 6.8.8 is the default minimization software (Connolly, 2005). If an alternative minimizer is used, some cosmetic changes to the code are required to ensure that the output can be read properly. If assistance is needed with this, please contact the corresponding author. 
Residual Calculation: 
The three components of the residual are the identity of the stable phases, phase proportions and phase compositions. Estimates for the phase proportions and phase compositions are dependent on one another as phase proportions are calculated directly from the phase compositions. Due to this a coefficient is used when phase proportions are calculated rather than measured, so that the component of the residual emanating from this is not counted twice. A large penalty is added to the residual if the difference between the expected and predicted number of stable phases exceeds 1. This discourages the algorithm from producing a more optimal fit for most of the data at the expense of one or more data points, as it is assumed that all inputted experimental data is accurate and should be reproduced by the thermodynamic model. If an experiment still has this penalty at the end of the optimization, remove it and recalculate the result. 
Bayesian Inference: 
Errors for the fitted parameters are assumed to be a truncated Gaussian distribution if a multi-start approach is not used. The mean values are the “fitted” results and the standard deviation is assumed to be the same as the prior. If a multi-start approach is used, then errors are treated as a uniform distribution between the minimum and maximum calculated values for each parameter. 
Underdetermined systems:
In underdetermined systems, systems in which there are fewer constraints than unknowns, it may prove more useful to optimize for smaller subgroups of parameters iteratively, rather than all parameters at once. This reduces the likelihood of a solution converging on a local minimum, in particular where parameters produce similar effects on the minimum free energy surface. 