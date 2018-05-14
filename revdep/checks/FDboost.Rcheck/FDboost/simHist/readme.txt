
Code for simulation study presented in 
Brockhaus, S., Melcher, M., Leisch, F. and Greven, S. (2017): 
Boosting flexible functional regression models with a high number of functional historical effects,  
Statistics and Computing, 27(4), 913-926.  


Use the run_...R files to run the different settings of the simulation study, 

- run_FDboost2.R fit the models without nuisance variables by FDboost  
- run_FDboost3.R fit the models with nuisance variables by FDboost, use stability selection
- run_FDboost4.R fit the models with nuisance variables by FDboost

- run_pffr2.R fit the models without nuisance variables by FAMM 
- run_pffr3.R try to fit the models with nuisance variables by FAMM (cannot fit the models, always gives error!)

- run_PlotModels.R get plots of the data settings and the estimated coefficients for FAMM and FDboost 

The code in analyze_results.R can be used to get plots and tables using the results. 

simUtils.R and boosting4.R contain utility functions. 


The parallelization only works on Linux, not on Windows. 