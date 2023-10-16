## Data-adaptive-functional-dimension-reduction

Code to reproduce simulation results from “A data-adaptive dimension reduction for functional data via penalized low-rank approximation” by Yeonjoo Park, Hee-Seok Oh, and Yaeji Lim.

#Code
-simData_gen.R: 
contain source function generating simulation data sets. Based on Section 3.1, there are two options for the model (model=1 or 2) and three options for error structures (setting=1,2 or 3). We can also set true dimension K (K.true=5 or 8), missing proportions (m.prob= 0, 0.1 or 0.2), asymmetricity parameter for skewed errors (symm), sample size (n), and the variance of model contamination error (sigma.con).

-comparison_methods.R: 
contain functions for the comparison methods in simulation study, Robust FPCA, Proj-Huber, and Partial FPCA. 

-help_functions.R: 
contain remaining help functions necessary for the implementation.

-FDR.R: 
contain proposed_FDR() to implement the proposed method following Section 2.3 practical algorithm. Also, it includes proposed_FDR_L2(), the non-robust version of the proposed method using squared loss errors instead of the composite loss error function. There are several options to choose the tuning parameter $\lambda$.
-	ad-hoc approach (Section 2.4.1)  by setting adhoc.CV=TRUE
-	sparsity threshold: if we set the desired sparsity of $\boldsymbol{\Phi}$ and set sparsity.level=c, then $\lambda$ will be chosen from values among the pre-specified lambda.sequence that first achieves the sparsity equal to or higher than c.
-	The number of dimension $q$ setting at the algorithm: as described in Section 2.4.2, we set q sufficiently large. In the simulation study, when K.true=5, we set q=15, and when K.true=8, we set q=20.

-simulation.R: 
produce results displayed in Figure 2 and Table 1. The tuning parameter lambda can be chosen by adhoc CV by adhoc.CV=TRUE. Or you can set sparsity.level = 0.75 to 0.85 if the sparsity is desired to be a specific level.


#Instructions for Use
1.	Install all packages in simulation.R
2.	Load all source and help functions: simData_gen.R, comparison_methods.R, other_functions.R and FDR.R
3.	Run simulation.R with interested model and error structure options.

