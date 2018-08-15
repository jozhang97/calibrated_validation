## Well-calibrated validation of biogeo package

`calibrated_simulation_biogeo.py` - call simulate_params.R, process files created by it, create .xml's for BEAST runs [output: xmls/*.xml]
`simulate_params.R` - draw from priors, simulate trees and states with diversitree [output: prior_samples.pdf, ntips.pdf, data_param_tree.csv, discarded_param_tree.csv]

1) BiSSE: drawing parameter samples from prior, simulating trees and making .xmls    

    ``python calibrated_simulation_biogeo.py -rd ./ -od results/ -p test -n 2000 -m 0.1 -sd 0.05 -b True -st 14``    

    ( calls simulate_params.R -- which draws from prior, plots graphs, simulates with diversitree, and makes .csv --  then .xml's are made from BiSSE .xml template.)    

1.1) (TO-DO) Submitting BEAST jobs
