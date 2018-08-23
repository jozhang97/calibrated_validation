## Well-calibrated validation of biogeo package

`simulate_prep4beast.py`: call simulate_params.R, process files created by it, create .xml's for BEAST runs [output: xmls/*.xml]    
`simulate_params.R`: draw from priors, simulate trees and states with diversitree [output: prior_samples.pdf, ntips.pdf, data_param_tree.csv, discarded_param_tree.csv]    
`parse_beast_logs.py`: call parse_log.R, process files created by it, summarize into results [output on screen: the frequency at which posteriors included true parameter value]    

1) BiSSE: drawing parameter samples from prior, simulating trees and making .xmls (note that the value after -n will affect the results, even if setting the seed in R)

    ``python simulate_prep4beast.py -od csvs_plots/ -m 0.4,0.4,0.2,0.2,0.1,0.1 -sd 0.05,0.05,0.05,0.05,0.05,0.05 -pn l0,l1,m0,m1,q01,q10 -p bisse -xd bisse_xmls/ -xt bisse_beast_template.xml -n 2000 -pd /N/u/fkmendes/Carbonate/Documents/uoa/calibrated_validation/``    
    
2) Submitting BEAST jobs (was done on cluster, I wrote a small python script to qsub all .pbs files. This step should produce a bunch of .log files that are put into beast_outputs/    

3) Parsing BEAST outputs    

    ``python parse_beast_logs.py -bd beast_outputs/ -rd ./ -cd csvs_plots/ -b 500000 -n 100 -n1 l0,l1,m0,m1,q01,q10 -n2 Lambda1,Lambda2,Mu1,Mu2,FlatQMatrix2,FlatQMatrix3``    
