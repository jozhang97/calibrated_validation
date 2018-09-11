import csv
import argparse
import os

def generate_priors_table(output_dir, prefix, num_states, speciation_events, prior_dists, prior_params):
    """ Write _prior_params.csv so it can be read by simulate_prep4beast """

    prior_dists_list = prior_dists.split(",")
    prior_params_list = prior_params.split(";")
    lambdas = speciation_events.split(",")

    priors = list() # where we'll store everything
    
    priors = ["m"+str(i)+"|"+prior_dists_list[0]+"|"+prior_params_list[0] for i in range(1,num_states+1)] # mu priors
        

    priors.extend([l+"|"+prior_dists_list[2+i]+"|"+prior_params_list[2+i] for i,l in enumerate(lambdas)]) # lambda priors; first two are death and transition parameters

    q_priors = list()
    for i in range(num_states):
        for j in range(num_states):
            if i != j: priors.extend(["q"+str(i)+str(j)+"|"+prior_dists_list[1]+"|"+prior_params_list[1]]) # q priors

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    with open(output_dir + prefix + "_prior_params.csv", "w") as priors_file:
        priors_file.write("param_name|prior_dist|moments\n")

        for p in priors:
            priors_file.write(p + "\n")

if __name__ == "__main__":
    # Note - default params will create BiSSE model
    parser = argparse.ArgumentParser(prog="Generate priors script", description="Generates _prior_params.csv file to be used in simulate_prep4beast.")
    parser.add_argument("-od", "--output-dir", action="store", dest="outputdir", default="./", type=str, help="Full file path to directory where simulations and log files will be saved.")
    parser.add_argument("-p", "--prefix", action="store", dest="prefix", default="", type=str, help="Prefix for result files.")
    parser.add_argument("-n", "--num-states", action="store", dest="numstates", default=2, type=int, help="Number of states.")
    parser.add_argument("-s", "--speciation_events", action="store", dest="spec_events", default="S", type=str, help="Speciation events used (e.g. SS,V,S)")
    parser.add_argument("-pt", "--prior-distributions", action="store", dest="prior_dists", default=None, type=str, help="Prior distributions. The first two are for the death and transition parameters.")
    parser.add_argument("-pp", "--prior-params", action="store", dest="prior_params", default=None, type=str, help="Parameters of prior distributions (sep is ; between parameters, and comma within parameters). First two are for the death and transition parameters.")
    args = parser.parse_args()

    generate_priors_table(args.outputdir, args.prefix, args.numstates, args.spec_events, args.prior_dists, args.prior_params)
