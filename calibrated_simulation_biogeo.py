import sys
import os
import subprocess
import csv
import argparse

def simulate(r_script_dir, output_dir, n_sims, is_bisse, prefix):
    """ Call simulation pipeline (r script 1, parsing, r script 2) """

    r_script_1 = r_script_dir + "simulate_params.R"
    r_script_2 = r_script_dir + "simulate_SSE.R"

    if is_bisse:
        subprocess.call (["Rscript", "--vanilla", r_script_1, output_dir,
                prefix, str(n_sims), "l1", "l2", "m1", "m2", "q12", "q21"])
    else:
        param_names = ["l111", "l112", "l22", "l211", "l221", "l222", "m1", "m2", "q12", "q21"]
        subprocess.call (["Rscript", "--vanilla", r_script_1, output_dir,
                prefix, str(n_sims)] + param_names)

    cmd_1 = ["Rscript", "--vanilla", r_script_2, output_dir, prefix]

    output_path = os.path.join(output_dir, prefix)
    with open(output_path + '_all_params.csv') as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')

        for i, row in enumerate(csv_reader):
            cmd_this_sim = cmd_1 + row
            cmd_this_sim[4] += str(i)
            subprocess.call(cmd_this_sim)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Simulation script", description="Script for performing well-calibrated validation of biogeo package.")
    parser.add_argument("-rd", "--rscripts-dir", action="store", dest="rscriptsdir", default="./", type=str, help="full file path to r scripts directory.")
    parser.add_argument("-od", "--output-dir", action="store", dest="outputdir", default="./", type=str, help="full file path to directory where simulations and log files will be saved.")
    parser.add_argument("-n", "--n-sims", action="store", dest="nsims", default=100, type=int, help="number of simulations.")
    parser.add_argument("-b", "--is-bisse", action="store", dest="bisse", default=True, type=bool, help="Flag for BiSSE simulations (default: True)")
    parser.add_argument("-p", "--prefix", action="store", dest="prefix", default="", type=str, help="Prefix for result files.")
    args = parser.parse_args()

    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)

    simulate(args.rscriptsdir, args.outputdir, args.nsims, args.bisse, args.prefix)

