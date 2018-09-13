## Well-calibrated validation of biogeo package

`simulate_prep4beast.py`: call simulate_params.R, process files created by it, create .xml's for BEAST runs [output: xmls/*.xml]    
`simulate_params.R`: draw from priors, simulate trees and states with diversitree [output: prior_samples.pdf, ntips.pdf, data_param_tree.csv, discarded_param_tree.csv]    
`parse_beast_logs.py`: call parse_log.R, process files created by it, summarize into results [output on screen: the frequency at which posteriors included true parameter value]    

1) BiSSE: drawing parameter samples from prior, simulating trees and making .xmls (note that the value after -n will affect the results, even if setting the seed in R)

    ``python simulate_prep4beast.py -od csvs_plots/ -pt 'exp,exp,exp,exp,exp,exp' -pp '20;20;80;80;100;100' -pn l0,l1,m0,m1,q01,q10 -p bisse -xd bisse_xmls/ -xt bisse_beast_template.xml -n 4000 -pd /N/u/fkmendes/Carbonate/Documents/uoa/calibrated_validation/ -st 50 -b # Should ignore 10 too-large simulations, and have a n-tip median of 10; if you don't get this, something went wrong with the seeding and your beast_outputs won't match.``    

    ``python simulate_prep4beast.py -od csvs_plots/ -pt 'exp,exp,exp,exp,exp,exp' -pp '30;30;100;100;100;100' -pn l0,l1,m0,m1,q01,q10 -p bisse -xd bisse_xmls/ -xt bisse_beast_template.xml -n 2000 -pd /N/u/fkmendes/Carbonate/Documents/uoa/calibrated_validation/ -st 50 -b # 356 tree were <30spp and 14 trees were >1000spp, these were ignored``    
    
    NOTES:    
    - When any of the parameters is larger than 1, coverage will be fine, but the power analysis will show that the posterior mean has little to no correlation to the true value. Correlation occurs only when parameter values are smaller than 1-ish. 

2) CLaSSE: generate prior parameters    

    ``python generate_priors_table.py -n 3 -od csvs_plots_classe/ -s SS,V,S -pt exp,exp,exp,exp,exp -pp '80;100;20;20;20' -p classe # death, transition, lambdas``    
        
    ``python simulate_prep4beast.py -od csvs_plots_classe/ -p classe -xd classe_xmls/ -xt classe_beast_template.xml -n 2000 -pd /N/u/fkmendes/Carbonate/Documents/uoa/calibrated_validation/ -ppf csvs_plots_classe/classe_prior_params.csv -e2t spec_event_to_triplet.csv -st 50``    

    ``python generate_priors_table.py -n 3 -od csvs_plots_classe/ -s SS,V,S -pt exp,exp,exp,exp,exp -pp '80;20;20;20;20' -p classe``    

    ``python simulate_prep4beast.py -od csvs_plots_classe/ -p classe -xd classe_xmls/ -xt classe_beast_template.xml -n 2000 -pd /nesi/project/nesi00390/fkmendes/calibrated_validation/ -ppf csvs_plots_classe/classe_prior_params.csv -e2t spec_event_to_triplet.csv -st 50``    
    
    (Note above line does not work right now since xml not availbale)    
    
2) Submitting BEAST jobs (was done on cluster, I wrote a small python script to qsub all .pbs files. This step should produce a bunch of .log files that are put into beast_outputs/    

3) Parsing BEAST outputs    

    ``python parse_beast_logs.py -bd beast_outputs/ -rd ./ -cd csvs_plots/ -b 500000 -n 100 -n1 l0,l1,m0,m1,q01,q10 -n2 Lambda1,Lambda2,Mu1,Mu2,FlatQMatrix1,FlatQMatrix2``    


