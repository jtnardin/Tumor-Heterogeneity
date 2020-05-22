# Tumor-Heterogeneity: Nonparametric Estimation

There are several directories present:

DataGeneration: This is where the data you are trying to fit goes

Figures: Generated figures will end up in this folder

Outputs: Location of .mat files after estimation is complete

PDECode: The files that solve the PDE (same as those in the NLME Section)

PrecomputedSolutions: Folder containing precomputed solutions for your parameter mesh and code to generate precomputed solutions

ProhorovFiles: Folder containing .m files for the optimization problem and using the Akaike Informaion Criteria (AIC) to choose the best parameter meshing strategy.


In addition, there are several files in the main directory:

driver.m: This is the main driver for the non-parametric estimation. It precomputes solutions, loads in the data, performs the optimization routine.

