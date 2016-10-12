# evacuation

This is the README for analysis of evacuation decision making.
------------------
Update 12 Oct 2016
------------------
New organization of directory has been implemented. In main 
directory are scripts that direct analysis:
 - evac_analysis     : initial dev and testing of functions
 - evac_bayesian_new : bayesian learning of strategies and 
 		       follow-up analysis
 - evac_grp_ind_crossval : comparison of ind and group behavior
                           and figures / analysis
 - evac_lossmat_analysis : evaluation of strategy on different
                           loss matrix combinations
 - evac_error_analysis : evaluation of error tolerances in 
                         bayesian solvers

Main directory also contains the master data sets:
 - evacuate_data.mat
 - experiment.mat

Finally, main diretory contains functions for running Bayesian
analysis on the cluster:
 - cluster_bayesian_parfor (setup, main cluster script)
 - bayesian_wrapper        (wraps parallel games and saves)
 - bayesian_game           (computes bayesian solutions)
 - bayesian_posteriors     (for use on local machine: combines
                            likelihoods from parallel computation
			    into bayesian posterior distributions)

