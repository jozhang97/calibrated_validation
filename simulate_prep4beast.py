import sys
import os
import subprocess
import csv
import argparse
import re
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
        self.xml = str()
        self.idx = idx # which simulation it is

    def update_xml(self, xml_template):
        self.xml = xml_template

    def update_project_dir(self, project_dir):
        self.project_dir = project_dir

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

        self.xml = self.xml.replace(key, replacements_str.replace(" \"", "\""))

    def populate_xml(self, xml_template, project_dir):
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
        self.update_project_dir(project_dir) # will add cluster path to BEAST output files if defined

        self.replace_in_xml("[Mean Lambda Prior]", self.prior_params[0][0], quoted=True)
        self.replace_in_xml("[Stdev Lambda Prior]", self.prior_params[0][1], quoted=True)
        self.replace_in_xml("[Mean Mu Prior]", self.prior_params[1][0], quoted=True)
        self.replace_in_xml("[Stdev Mu Prior]", self.prior_params[1][1], quoted=True)
        self.replace_in_xml("[Mean FlatQMatrix Prior]", self.prior_params[2][0], quoted=True)
        self.replace_in_xml("[Stdev FlatQMatrix Prior]", self.prior_params[2][1], quoted=True)
        self.replace_in_xml("[Tree in newick format]", self.tree)
        self.replace_in_xml("[List of species]", convert_to_species_list(self.tip_states))
        self.replace_in_xml("[Initial transition rate values]", convert_to_flat_rate_matrix(self.init_params[4:6]))
        self.replace_in_xml("[Initial lambda values]", self.init_params[0:2])
        self.replace_in_xml("[Initial mu values]", self.init_params[2:4])
        self.replace_in_xml("[Pi values]", "0.0 0.0 0.5 0.5")
        self.replace_in_xml("[Species=Trait State]", self.tip_states)

        beast_output_path = "beast_outputs/"
        if self.project_dir:
            beast_output_path = project_dir + "beast_outputs/"
        self.replace_in_xml("[Simulation log file name]", beast_output_path + self.prefix + str(self.idx) + "_sim.log")
        self.replace_in_xml("[Simulation tree file name]", beast_output_path + self.prefix + str(self.idx) + "_sim.trees")

    def write_xml(self, file_name):
        with open(file_name, "w") as f:
            f.write(self.xml)

def simulate(r_script_dir, output_dir, n_sims, is_bisse, prefix, sim_time, pnames, mu, std):
    """ Call simulation pipeline (create .csvs and plots) """
    
    output_path = os.path.join(output_dir, prefix)
    param_names = list()

    r_script = r_script_dir + "simulate_params.R"
    cmd = ["Rscript", "--vanilla", r_script, output_dir, prefix, str(n_sims), mu, std, pnames, str(sim_time)]
    subprocess.call(cmd)
    print " ".join(cmd)
    
def parse_simulations(output_dir, xml_dir, xml_template_name, prefix, prior_params, project_dir):
    """ Parse .csv file into .xmls """
    num_params = len(prior_params)*2

    csv_info_stashes = list()

    csv_true = output_dir + "data_param_tree.csv"
    csv_inits = output_dir + "data_param_inits.csv"
    with open(csv_true) as tree_file, open(csv_inits) as inits_file:
        true_reader = csv.reader(tree_file, delimiter='|')
        inits_reader = csv.reader(inits_file, delimiter='|')

        # reading data and init csvs at the same time
        for i, (true_row, init_params) in enumerate(zip(true_reader, inits_reader)):
            if i == 0: continue # ignore first line since its the headers
            true_params = true_row[:num_params] # these will not be used as they are what we are estimating
            # TODO this indexing may change in the future
            tree = true_row[num_params]
            n_tips = true_row[num_params + 1]
            tip_states = true_row[num_params + 2]

            csv_info_stash = CsvInfoStash(init_params, prior_params, tree, n_tips, tip_states, output_dir, prefix, i)
            csv_info_stashes.append(csv_info_stash)

    with open(xml_template_name, 'r') as xml_template_file:
        xml_template = xml_template_file.read()

    # iterating over list of csv_info_stashes
    for i, csv_info_stash in enumerate(csv_info_stashes):
        csv_info_stash.populate_xml(xml_template, project_dir) # fill out xml template
        xml_file_name = xml_dir + prefix + str(i+1) + ".xml"
        csv_info_stash.write_xml(xml_file_name)

    return

def write_pbs(xml_dir, prefix, project_dir):
    """ Write one .pbs script per .xml """

    sim_n_regex = re.compile("[0-9]+")
    xml_file_names = [f for f in os.listdir(xml_dir) if f.endswith(".xml")]
    for xml_file_name in xml_file_names:
        sim_n = re.findall(sim_n_regex, xml_file_name.split("_")[0])[0]

        with open("pbs_scripts/" + prefix + sim_n + ".PBS", "w") as pbs_file:
            pbs_file.write("#!/bin/bash\n#PBS -N beast_" + sim_n + \
                           "\n#PBS -l nodes=1:ppn=1,walltime=96:00:00\n#PBS -M fkmendes@iu.edu\n#PBS -m abe\n\njava -jar " + project_dir + "biogeo.jar " + \
                           project_dir + xml_dir + xml_file_name
            )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Simulation script", description="Script for performing well-calibrated validation of biogeo package.")
    parser.add_argument("-rd", "--rscripts-dir", action="store", dest="rscriptsdir", default="./", type=str, help="Full file path to r scripts directory.")
    parser.add_argument("-od", "--output-dir", action="store", dest="outputdir", default="./", type=str, help="Full file path to directory where simulations and log files will be saved.")
    parser.add_argument("-xd", "--xml-dir", action="store", dest="xmldir", default="./", type=str, help="Full file path to directory where .xml's will be saved.")
    parser.add_argument("-n", "--n-sims", action="store", dest="nsims", default=1000, type=int, help="Number of simulations.")
    parser.add_argument("-b", "--is-bisse", action="store", dest="bisse", default=True, type=bool, help="Flag for BiSSE simulations (default: True)")
    parser.add_argument("-p", "--prefix", action="store", dest="prefix", default="", type=str, help="Prefix for result files.")
    parser.add_argument("-st", "--sim-time", action="store", dest="simtime", default=10, type=float, help="Time to run simulation.")
    parser.add_argument("-pn", "--param-names", action="store", dest="pnames", default=None, type=str, help="Parameter names (order matters for means and stdevs of priors).")
    parser.add_argument("-m", "--mu", action="store", dest="mu", default=None, type=str, help="Means of prior distns, in the order parameter names are passed.")
    parser.add_argument("-sd", "--std", action="store", dest="std", default=None, type=str, help="Standard deviations of prior distns, in the order parameter names are passed.")
    parser.add_argument("-xt", "--xml-template", action="store", dest="xmlt", default=None, type=str, help="Full path to template of .xml file.")
    parser.add_argument("-pd", "--project-dir", action="store", dest="projdir", default=None, type=str, help="Full path to calibration folder if -pbs.")
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

    simulate(args.rscriptsdir, args.outputdir, args.nsims, args.bisse, args.prefix, args.simtime, args.pnames, args.mu, args.std) # calls R script, produces .csv files and plots

    mus = args.mu.split(",")
    stds = args.std.split(",")
    idxs = [i for i in xrange(0, len(args.pnames.split(",")),2)]
    prior_params = zip([mus[i:(i+2)] for i in idxs], [stds[i:(i+2)] for i in idxs])
    parse_simulations(args.outputdir, args.xmldir, args.xmlt, args.prefix, prior_params, args.projdir) # parses .csv into .xml files

    if args.projdir:
        if not os.path.exists("pbs_scripts"):
            os.makedirs("pbs_scripts")
        write_pbs(args.xmldir, args.prefix, args.projdir)
