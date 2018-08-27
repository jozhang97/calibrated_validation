# This script simulates a csv file with the following rows
# param_name, distribution_type, moments (e.g. mu, sigma or mu)
import csv
import argparse

def write_csv(priors_table, output_dir, prefix):
    file_name = output_dir + prefix + "_prior_params.csv"
    header = ["param_name", "dist_type", "moments"]
    print("Writing to prior params to " + str(file_name))
    with open(file_name, 'w') as fp:
        writer = csv.writer(fp, delimiter='|')
        writer.writerow(header)
        for row in priors_table:
            writer.writerow(row)
    return

def combine(obj_lst):
    ret_str = ""
    for obj in obj_lst:
        ret_str += str(obj)
    return ret_str

def generate_priors_table(output_dir, prefix, num_states, speciation_events):
    priors = []
    for i in range(1, num_states+1):
        priors.append([combine(["m", i]), "lnorm", "1,1"])  # mu
        for j in range(1, num_states+1):
            if i != j:
                priors.append([combine(["q", i, j]), "lnorm", "1,0.5"])  # Q

    for event in speciation_events:
        priors.append([event, "exp", "3"])

    write_csv(priors, output_dir, prefix)


if __name__ == "__main__":
    # Note - default params will create BiSSE model
    parser = argparse.ArgumentParser(prog="Generate priors script", description="Generates csv file to be used in simulate_prep4beast")
    parser.add_argument("-od", "--output-dir", action="store", dest="outputdir", default="./", type=str, help="Full file path to directory where simulations and log files will be saved.")
    parser.add_argument("-p", "--prefix", action="store", dest="prefix", default="", type=str, help="Prefix for result files.")
    parser.add_argument("-n", "--num-states", action="store", dest="numstates", default=2, type=int, help="Number of states")
    parser.add_argument("-s", "--speciation_events", action="store", dest="specevents", default="S", type=str, help="Speciation events used (e.g. SS,V,S)")
    # Tweek the prior dist and params
    args = parser.parse_args()

    spec_events = args.specevents.split(",")
    generate_priors_table(args.outputdir, args.prefix, args.numstates, spec_events)
