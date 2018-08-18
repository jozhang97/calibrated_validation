import sys
import os
import subprocess
import csv
import argparse
import numpy as np

def convert_to_species_list(tip_states):
    def construct_taxon_xml(tip_name):
        taxon = "<taxon id=\""
        taxon += tip_name
        taxon += "\" spec=\"Taxon\"/>"
        return taxon

    taxon_set = str()
    for tip_state in tip_states.split(","):
        tip_name = tip_state.split("=")[0]
        taxon_set += construct_taxon_xml(tip_name)
        taxon_set += "\n \t \t"

    return taxon_set

def convert_to_flat_rate_matrix(q):
    # solve n * n - n = len(q)
    if len(q) == 2:  # BiSSE
        return [1e-10, q[0], q[1], 1e-10]
    raise NotImplementedError
    #roots = np.roots(1, -1, len(q))
    #n = max(roots)  # not sure
    # TODO ... this is kinda annoying, keep looking for a library

class CsvInfoStash:
    def __init__(self, init_params, prior_params, tree, n_tips, tip_states, output_dir, prefix, idx):
        self.init_params = init_params
        self.prior_params = prior_params
        self.tree = tree
        self.n_tips = n_tips
        self.tip_states = tip_states
        self.output_dir = output_dir
        self.prefix = prefix
        self.xml = ""
        self.idx = idx # which simulation it is
        
    def update_xml(self, xml_template):
        self.xml = xml_template

    def replace_in_xml(self, key, replacements, quoted=False):
        """ Fill out xml template
        
        key: Keyword to look for in xml
        replacements: list of objects (with __str__ implemented) to replace with OR a single object
        quoted: if True, adds quotes to the string
        """
        replacements_str = str()
        if type(replacements) == type([]):
            for replacement in replacements:
                replacements_str += str(replacement) + " "

        else:
            replacements_str += str(replacements)

        if quoted:
            replacements_str = "\"" + replacements_str + "\""

        self.xml = self.xml.replace(key, replacements_str)

    def populate_xml(self, xml_template):
        """
        keys = [
        "[Mean Lambda Prior]",
        "[Stdev Lambda Prior]",
        "[Mean Mu Prior]",
        "[Stdev Mu Prior]",
        "[Mean FlatQMatrix Prior]",
        "[Stdev FlatQMatrix Prior]",
        "[Tree in newick format]",
        "[List of species]",
        "[Initial transition rate values]",
        "[Initial lambda values]",
        "[Initial mu values]",
        "[Pi values]",
        "[Species=Trait State]"
        ]
        """
        self.update_xml(xml_template) # initialize xml member

        self.replace_in_xml("[Mean Lambda Prior]", self.prior_params[0], quoted=True)
        self.replace_in_xml("[Stdev Lambda Prior]", self.prior_params[1], quoted=True)
        self.replace_in_xml("[Mean Mu Prior]", self.prior_params[2], quoted=True)
        self.replace_in_xml("[Stdev Mu Prior]", self.prior_params[3], quoted=True)
        self.replace_in_xml("[Mean FlatQMatrix Prior]", self.prior_params[4], quoted=True)
        self.replace_in_xml("[Stdev FlatQMatrix Prior]", self.prior_params[5], quoted=True)
        self.replace_in_xml("[Tree in newick format]", self.tree)
        self.replace_in_xml("[List of species]", convert_to_species_list(self.tip_states))
        self.replace_in_xml("[Initial transition rate values]", convert_to_flat_rate_matrix(self.init_params[4:6]))
        self.replace_in_xml("[Initial lambda values]", self.init_params[0:2])
        self.replace_in_xml("[Initial mu values]", self.init_params[2:4])
        self.replace_in_xml("[Pi values]", "0.0 0.0 0.5 0.5")
        self.replace_in_xml("[Species=Trait State]", self.tip_states)
        self.replace_in_xml("[Simulation log file name]", "beast_outputs/" + self.prefix + str(self.idx) + "_sim_log.log")
        self.replace_in_xml("[Simulation tree file name]", "beast_outputs/" + self.prefix + str(self.idx) + "_sim_trees.trees")

    def write_xml(self, file_name):
        with open(file_name, "w") as f:
            f.write(self.xml)

def simulate(r_script_dir, output_dir, n_sims, is_bisse, prefix, sim_time, mu, std):
    """ Call simulation pipeline (create .csvs and plots) """

    output_path = os.path.join(output_dir, prefix)
    param_names = list()

    r_script = r_script_dir + "simulate_params.R"
    cmd = ["Rscript", "--vanilla", r_script, output_dir, prefix, str(n_sims), str(mu), str(std), str(sim_time)]
    if is_bisse:
        param_names = ["l0", "l1", "m0", "m1", "q01", "q10"]
        cmd_bisse = cmd + param_names
        print " ".join(cmd_bisse)
        subprocess.call(cmd_bisse)
        # setup_data_csv(output_path, param_names)

    else:
        param_names = ["l111", "l112", "l22", "l211", "l221", "l222", "m1", "m2", "q12", "q21"]
        cmd_classe = cmd + param_names
        subprocess.call(cmd_classe)

def parse_simulations(output_dir, xml_dir, xml_template_name, prefix, prior_params):
    """ Parse .csv file into .xmls """
    csv_info_stashes = list()

    csv_true = output_dir + "data_param_tree.csv"
    csv_inits = output_dir + "data_param_inits.csv"
    with open(csv_true) as tree_file, open(csv_inits) as inits_file:
        true_reader = csv.reader(tree_file, delimiter='|')
        inits_reader = csv.reader(inits_file, delimiter='|')

        # reading data and init csvs at the same time
        for i, (true_row, init_params) in enumerate(zip(true_reader, inits_reader)):
            if i == 0: continue # ignore first line since its the headers
            true_params = true_row[:6] # these will not be used as they are what we are estimating
            tree = true_row[-3]
            n_tips = true_row[-2]
            tip_states = true_row[-1]

            csv_info_stash = CsvInfoStash(init_params, prior_params, tree, n_tips, tip_states, output_dir, prefix, i)
            csv_info_stashes.append(csv_info_stash)

    with open(xml_template_name, 'r') as xml_template_file:
        xml_template = xml_template_file.read()

    # iterating over list of csv_info_stashes 
    for i, csv_info_stash in enumerate(csv_info_stashes):
        csv_info_stash.populate_xml(xml_template) # fill out xml template
        xml_file_name = xml_dir + prefix + str(i+1) + ".xml"
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
    except:
        exit("Could not find/open xml template. Exiting...")

    if not os.path.exists(args.outputdir):
        os.makedirs(args.outputdir)

    if not os.path.exists(args.xmldir):
        os.makedirs(args.xmldir)

    if not os.path.exists("beast_outputs"):
        os.makedirs("beast_outputs")
        
    # simulate(args.rscriptsdir, args.outputdir, args.nsims, args.bisse, args.prefix, args.simtime, args.mu, args.std) # calls R script, produces .csv files and plots

    prior_params = [args.mu, args.std] * 3
    parse_simulations(args.outputdir, args.xmldir, args.xmlt, args.prefix, prior_params) # parses .csv into .xml files
