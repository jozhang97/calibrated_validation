import sys
import os
import subprocess
import csv
import argparse

class CsvInfoStash:
    def __init__(self, params, tree, n_tips, tip_states):
        self.params = params
        self.tree = tree
        self.n_tips = n_tips
        self.tip_states = tip_states

    def populate_xml(self, xml_template):
        # TODO
        self.xml = "xml"

    def write_xml(self, file_name):
        with open(file_name, "w") as f:
            f.write(self.xml)

def simulate(r_script_dir, output_dir, n_sims, is_bisse, prefix, sim_time, mu, std):
    """ Call simulation pipeline (r script 1, r script 2) """

    output_path = os.path.join(output_dir, prefix)
    param_names = list()

    # sampling 100 parameter values for all params
    r_script_1 = r_script_dir + "simulate_params.R"
    cmd_1 = ["Rscript", "--vanilla", r_script_1, output_dir, prefix, str(n_sims), str(mu), str(std), str(sim_time)]
    if is_bisse:
        param_names = ["l0", "l1", "m0", "m1", "q01", "q10"]
        cmd_bisse = cmd_1 + param_names
        print " ".join(cmd_bisse)
        subprocess.call(cmd_bisse)
        # setup_data_csv(output_path, param_names)

    else:
        param_names = ["l111", "l112", "l22", "l211", "l221", "l222", "m1", "m2", "q12", "q21"]
        cmd_classe = cmd_1 + param_names
        subprocess.call(cmd_classe)

def parse_simulations(output_dir, xml_dir, xml_template):
    """ Parse .csv file into .xmls """
    csv_info_stashs = []

    csv_tree = output_dir + "data_param_tree.csv"
    csv_inits = output_dir + "data_param_inits.csv"
    with open(csv_tree) as tree_file, open(csv_inits) as inits_file:
        tree_reader = csv.reader(tree_file, delimiter=',')
        inits_reader = csv.reader(inits_file, delimiter=',')
        for i, (tree_row, init_params) in enumerate(zip(tree_reader, inits_reader)):
            if i == 0:
                # ignore first line since its the headers
                continue
            true_params = tree_row[:-3]
            tree = tree_row[-3]
            n_tips = tree_row[-2]
            tip_states = tree_row[-1]

            # Not sure on first argument
            csv_info_stash = CsvInfoStash(init_params, tree, n_tips, tip_states)
            csv_info_stashs.append(csv_info_stash)


    for i, csv_info_stash in enumerate(csv_info_stashs):
        csv_info_stash.populate_xml(xml_template)
        xml_file_name = xml_dir + str(i) + ".xml"
        csv_info_stash.write_xml(xml_file_name)

    return


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Simulation script", description="Script for performing well-calibrated validation of biogeo package.")
    parser.add_argument("-rd", "--rscripts-dir", action="store", dest="rscriptsdir", default="./", type=str, help="full file path to r scripts directory.")
    parser.add_argument("-od", "--output-dir", action="store", dest="outputdir", default="./", type=str, help="full file path to directory where simulations and log files will be saved.")
    parser.add_argument("-xd", "--xml-dir", action="store", dest="xmldir", default="./", type=str, help="full file path to directory where .xml's will be saved.")
    parser.add_argument("-n", "--n-sims", action="store", dest="nsims", default=1000, type=int, help="number of simulations.")
    parser.add_argument("-b", "--is-bisse", action="store", dest="bisse", default=True, type=bool, help="Flag for BiSSE simulations (default: True)")
    parser.add_argument("-p", "--prefix", action="store", dest="prefix", default="", type=str, help="Prefix for result files.")
    parser.add_argument("-st", "--sim-time", action="store", dest="simtime", default=10, type=float, help="Time to run simulation.")
    parser.add_argument("-m", "--mu", action="store", dest="mu", default=0, type=float, help="Mean of lognormal dist.")
    parser.add_argument("-sd", "--std", action="store", dest="std", default=0.05, type=float, help="Stdev of lognormal dist.")
    parser.add_argument("-xt", "--xml-template", action="store", dest="xmlt", default=None, type=str, help="Full path to template of .xml file.")
    args = parser.parse_args()

    xml_str = str()
    try:
        with open(args.xmlt, "r") as xml_template:
            xml_str = xml_template.readlines()
    except IOError:
        exit("Could not find file "+xmlt+". Exiting...")

    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)

    if not os.path.exists(args.xmldir):
        os.makedirs(args.xmldir)

    simulate(args.rscriptsdir, args.outputdir, args.nsims, args.bisse, args.prefix, args.simtime, args.mu, args.std) # calls R script, produces .csv files and plots

    parse_simulations(args.outputdir, args.xmldir, args.xmlt) # parses .csv into .xml files
