## Calibrated Validation

Calling `python simulate.py` will start by calling the `simulate_params.R` 
script followed by calling the `simulate_SSE.R` script


`validate_CLaSSE.R` - verifies that tree.classe outputs a tree as expected
    by comparing it to tree.bisse

`simulate_params.R` - samples N parameters NUM_SIM times and writes them
    to `all_params.csv` file. (Currently, this requires the name of the
    parameters as arguments, but it seems sufficient just to ask for the
    number of parameters. Will change soon)

`simulate_SSE.R` - this will run bisse or classe forward in time to generate
    a phylogenetic tree given a set of parameters. Note that if the 
    extinction rates are too high, no tree forms (reasonably). 
    Note that this does not work in batch. It is only for a single param set
    In future, can have this read a .csv 

`simulate.py` - calls `simulate_params` and parses the `all_params.csv` 
    Then, one parameter set at a time, we simulate bisse/classe and 
    write the data out

There are three files we output: `run.tree`, `run_tips.csv`, `run_params.csv`


Next steps
1. Get CLaSSE simulation working (currently only 2 character states)
2. Get the median number of taxa to a well defined number

