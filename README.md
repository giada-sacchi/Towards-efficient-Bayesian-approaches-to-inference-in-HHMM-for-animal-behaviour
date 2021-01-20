# Parallel-Tempering-on-Hierarchical-Hidden-Markov-Models

The datasets analysed for this study can be found at https://link.springer.com/article/10.1007/s13253-017-0282-9

The R scripts are provided below in chronological order of development (i.e. as presented in the paper).

----------------------------------------

mllk : Negative log-likelihood function and needed packages (based on the code produced by Leos-Barajas)

par.vec0 : Vector of parameter estimates obtained through MLE

MH : Single-update MH algorithm

MH_fun : Function for updating the parameter vector (related to MH)

Block1 : Block-update MH algorithm (Version 1: update of the means of each variable - simultaneously for all the states- is followed by the update of the corresponding standard deviations)

Block1_fun : Function for updating the parameter vector (related to Block1)

Block2 : Block-update MH algorithm (Version 2: update of the means of each all the covariates - simultaneously for all the states- is followed by the update of the standard deviations)

Block2_fun : Function for updating the parameter vector (related to Block1)

MH1 : Single-update MH algorithm (Version 1: related to the ordering adopted in Block1)

MH1_fun : Function for updating the parameter vector (related to MH1)

MH2 : Single-update MH algorithm (Version 2: related to the ordering adopted in Block2)

MH2_fun : Function for updating the parameter vector (related to MH2)

MH_tot : Single-update MH algorithm for the complete parameter set (including transition probabilities)

MH_tot_fun : Function for updating the parameter vector (related to MH_tot)

MH_temp* : MH algorithm with parallel tempering (code not parallelised)

MH_temp_fun* : Function for updating the parameter vector, function for computing the sum of prior distributions (related to MH_temp)

Temp_parallel : MH algorithm with parallel tempering running as many chains as the number of cores minus 1 (parallelised version)

Temp_parallel_fun : Functions for updating the parameter vector, computing the sum of prior distributions, updating of the parameter vector over the total number of iterations (related to Temp_parallel)

-----------------
*Not used in order to derive inference, but as a temporary step for the implementation of the full code.
