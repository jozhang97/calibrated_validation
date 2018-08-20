import subprocess
import argparse
import collections

def parse_beast_logs(beast_output_directory, r_script_dir, csv_dir, n_sims, burnin):
    """ Call R script that parses BEAST outputs and writes .csv file """

    r_script = r_script_dir + "parse_log.R"
    cmd = ["Rscript", "--vanilla", r_script, beast_output_directory, burnin, csv_dir]
    print "\n" + " ".join(cmd) + "\n"
    subprocess.call(cmd)

def summarize_validation(csv_dir, param_names1, param_names2):
    """ Summarize original true value csv file and hpd csv file """

    # hpd_idxs = dict((names2, idx) for idx, name2 in enumerate(param_names2.split(",")))
    tup_dict = dict(zip(param_names1.split(","), param_names2.split(",")))
    hpds_dict = dict()
    with open(csv_dir + "hpds.csv", "r") as hpds_file:
        hpds_idx = dict((idx, param) for idx, param in enumerate(hpds_file.readline().rstrip().split(",")))

        for idx, line in enumerate(hpds_file):
            tokens = line.rstrip().split(",")

            if not idx in hpds_dict: hpds_dict[idx] = dict()

            for idx2, token in enumerate(tokens[:-1]):
                param = hpds_idx[idx2]
                hpds_dict[idx][param] = float(token) # nested dicts; outer keys are sim number; inner keys are parameters lower and upper limits of hpds

    with open(csv_dir + "data_param_tree.csv", "r") as true_params_file:
        true_params_list = [param for param in param_names1.split(",")] # names of params
        success_dict = dict((param,0) for param in param_names1.split(",")) # storing final result
        for true_param_name in true_params_list:
            success_dict["ContainedWithin" + true_param_name] = 0
            success_dict["ShiftLeft" + true_param_name ]  = 0
            success_dict["ShiftRight" + true_param_name ]  = 0
            success_dict["Contains" + true_param_name] = 0
            success_dict["Outside" + true_param_name] = 0
            success_dict["WidthSmaller" + true_param_name] = 0
        true_idxs = dict((param, idx) for idx, param in enumerate(true_params_file.readline().rstrip().split("|"))) # cols of params
        total_sims = 0

        for idx, line in enumerate(true_params_file):
            total_sims += 1
            tokens = line.rstrip().split("|")
            # print "Simulation " + str(idx)
            try:
                prior_lower = float(tokens[true_idxs["hdplower"]])
                prior_upper = float(tokens[true_idxs["hdpupper"]])
            except:
                print("HDP of prior not given ")
                prior_lower = 0
                prior_upper = 0

            for true_param_name in true_params_list:
                col = true_idxs[true_param_name]
                param_value = float(tokens[col])
                post_lower = hpds_dict[idx]["lower"+tup_dict[true_param_name]] # idx is the sim number
                post_upper = hpds_dict[idx]["upper"+tup_dict[true_param_name]]
                # print true_param_name, str(param_value), lower, upper

                if param_value >= post_lower and param_value <= post_upper:
                    success_dict[true_param_name] += 1

                if prior_lower < post_lower and post_upper < prior_upper:
                    success_dict["ContainedWithin" + true_param_name] += 1

                if (post_lower < prior_lower and prior_upper < post_upper):
                    success_dict["Contains" + true_param_name] += 1


                if (post_lower < prior_lower and post_upper < prior_upper and prior_lower < post_upper):
                    success_dict["ShiftLeft" + true_param_name] += 1

                if (prior_lower < post_lower and prior_upper < post_upper and post_lower < prior_upper):
                    success_dict["ShiftRight" + true_param_name] += 1

                if (prior_upper < post_lower or post_upper < prior_lower):
                    success_dict["Outside" + true_param_name] += 1

                post_width = post_upper - post_lower
                prior_width = prior_upper - prior_lower

                if post_width < prior_width:
                    success_dict["WidthSmaller" + true_param_name ]  += 1


        for k,v in sorted(success_dict.items()):
            print "Parameter "  + str(k) + ": " + str(v) + "/" + str(total_sims)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Simulation script", description="Script for parsing BEAST .log files.")
    parser.add_argument("-bd", "--beast-dir", action="store", dest="beastdir", default="./", type=str, help="Full file path to beast outputs directory.")
    parser.add_argument("-rd", "--rscripts-dir", action="store", dest="rscriptsdir", default="./", type=str, help="Full file path to r scripts directory.")
    parser.add_argument("-cd", "--csv-dir", action="store", dest="csvdir", default="./", type=str, help="Full file path to directory where .csvs from simulations were saved.")
    parser.add_argument("-n", "--n-sims", action="store", dest="nsims", default=1000, type=int, help="Number of simulations.")
    parser.add_argument("-p", "--prefix", action="store", dest="prefix", default="", type=str, help="Prefix for result files.")
    parser.add_argument("-b", "--burnin", action="store", dest="burnin", default="", type=str, help="Number of burnin samples.")
    parser.add_argument("-n1", "--names1", action="store", dest="names1", default="", type=str, help="Names of parameters in csv file with true values, separated by comma.")
    parser.add_argument("-n2", "--names2", action="store", dest="names2", default="", type=str, help="Names of matching (with -n1) parameters in hpd csv file, separated by comma.")

    args = parser.parse_args()

    parse_beast_logs(args.beastdir, args.rscriptsdir, args.csvdir, args.nsims, args.burnin)

    summarize_validation(args.csvdir, args.names1, args.names2)
