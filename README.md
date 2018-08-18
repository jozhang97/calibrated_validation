## Well-calibrated validation of biogeo package

`calibrated_simulation_biogeo.py`: call simulate_params.R, process files created by it, create .xml's for BEAST runs [output: xmls/*.xml]    
`simulate_params.R`: draw from priors, simulate trees and states with diversitree [output: prior_samples.pdf, ntips.pdf, data_param_tree.csv, discarded_param_tree.csv]

1) BiSSE: drawing parameter samples from prior, simulating trees and making .xmls    

    ``python calibrated_simulation_biogeo.py -od csvs_plots/ -m 0.15 -sd 0.05 -p bisse -xd bisse_xmls/ -xt bisse_beast_template.xml -n 2000``    

1.1) (TO-DO) Submitting BEAST jobs
