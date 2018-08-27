import sys
import os
import subprocess
import csv
import argparse
import re
import numpy as np

def  parse_param_event_files(prior_params_file_name, event_to_triplet):
    """
    Parses the prior params and event triplet files to supply the data in req data format.
    Will convert the speciation event in prior params to the specific state changes in event_to_triplet
    """
    # Preparing mapping
    e2t_mapping = {}
    with open(event_to_triplet) as e2t_file:
        e2t = csv.reader(e2t_file, delimiter='|')
        for i, event_triplet in enumerate(e2t):
            if i == 0: continue
            event_name, triplet = event_triplet
            if event_name in e2t_mapping:
                e2t_mapping[event_name].append(triplet)
            else:
                e2t_mapping[event_name] = [triplet]


    # Read the table to write the param list
    param_names = []
    prior_dists = []
    prior_params = []
    with open(prior_params_file_name) as prior_params_file:
        prior_params_r = csv.reader(prior_params_file, delimiter='|')
        for i, prior in enumerate(prior_params_r):
            if i == 0: continue
            param_name, dist_type, moments = prior
            if param_name.startswith("m") or param_name.startswith("q"):
                param_names.append(param_name)
                prior_dists.append(dist_type)
                prior_params.append(moments)
            else:
                if param_name in e2t_mapping:
                    for triplet in e2t_mapping[param_name]:
                        spec_name = "l" + triplet.replace(",", "")
                        param_names.append(spec_name)
                        prior_dists.append(dist_type)
                        prior_params.append(moments)

    # Add unused speciation params
    # ?

    # Order the params
    def sort_by_list(list_to_sort, index_lst):
        ret = [list_to_sort[i] for i in index_lst]
        return ret

    idxs = np.argsort(param_names)
    param_names = sort_by_list(param_names, idxs)
    prior_dists = sort_by_list(prior_dists, idxs)
    prior_params = sort_by_list(prior_params, idxs)

    return param_names, prior_dists, prior_params

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

def stringfy_prior_params(prior_dist, prior_params, param_name, quoted=False):
    """ Generate prior distribution parameter string to put in BEAST .xml """

    # prior_dist = e.g. ["lnorm", "lnorm"]
    # prior_params = e.g. ["0.1,0.05" "0.1,0.05"] they will always be repeated
    prior_param_string = str()
    prior_params_list = prior_params[0].split(",")
    prior_dist = prior_dist[0].split(",")[0]

    # lognormal prior
    if prior_dist == "lnorm":
        mu = prior_params_list[0]
        std = prior_params_list[1]

        if quoted: mu = "\"" + mu + "\""; std = "\"" + std + "\""

        prior_param_string = "<distr id=\"LogNormal." + param_name + "\" spec=\"beast.math.distributions.LogNormalDistributionModel\" offset=\"0.0\" meanInRealSpace=\"false\">\n"
        prior_param_string += "\t\t<parameter name=\"M\" value=" + mu + " estimate=\"false\"/>\n"
        prior_param_string += "\t\t<parameter name=\"S\" value=" + std + " estimate=\"false\"/>\n\t  </distr>\n"

    # exponential prior
    elif prior_dist == "exp":
        rate = 1/float(prior_params_list[0])

        if quoted: rate = "\"" + str(rate) + "\""
        
        prior_param_string = "\t<distr id=\"Exponential." + param_name + "\" spec=\"beast.math.distributions.Exponential\" offset=\"0.0\" mean=" + rate + "/>\n"

    return prior_param_string

class CsvInfoStash:
    def __init__(self, init_params, param_names, prior_dists, prior_params, tree, n_tips, tip_states, output_dir, prefix, idx):
        self.init_params = init_params
        self.param_names = param_names
        self.prior_dists = prior_dists
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
        "[Lambda Prior Parameters]",
        "[Mu Prior Parameters]",
        "[FlatQMatrix Prior Parameters]",
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

        self.replace_in_xml("[Lambda Prior Parameters]", stringfy_prior_params(self.prior_dists[:2], self.prior_params[:2], "Lambda", quoted=True)) # note that quoted applies to stringfy function
        self.replace_in_xml("[Mu Prior Parameters]", stringfy_prior_params(self.prior_dists[3:5], self.prior_params[3:5], "Mu", quoted=True))
        self.replace_in_xml("[FlatQMatrix Prior Parameters]", stringfy_prior_params(self.prior_dists[4:], self.prior_params[4:], "FlatQMatrix", quoted=True))
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

def simulate(r_script_dir, output_dir, n_sims, is_bisse, prefix, sim_time, pnames, prior_dists, prior_params):
    """ Call simulation pipeline (create .csvs and plots) """

    output_path = os.path.join(output_dir, prefix)
    param_names = list()

    r_script = r_script_dir + "simulate_params.R"
    cmd = ["Rscript", "--vanilla", str(r_script), output_dir, prefix, str(n_sims), pnames, prior_dists, prior_params, str(sim_time)]
    print "\nR command call: " + " ".join(cmd) + "\n"
    subprocess.call(cmd)

def parse_simulations(output_dir, xml_dir, xml_template_name, prefix, param_names, prior_dists, prior_params, project_dir):
    """ Parse .csv file into .xmls """

    num_params = len(param_names)
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

            csv_info_stash = CsvInfoStash(init_params, param_names, prior_dists, prior_params, tree, n_tips, tip_states, output_dir, prefix, i)
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
                           "\n#PBS -l nodes=1:ppn=1,walltime=72:00:00\n#PBS -M fkmendes@iu.edu\n#PBS -m abe\n\njava -jar " + project_dir + "biogeo.jar " + \
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
    parser.add_argument("-pt", "--prior-distributions", action="store", dest="prior_dists", default=None, type=str, help="Prior distributions (one per parameter).")
    parser.add_argument("-pp", "--prior-params", action="store", dest="prior_params", default=None, type=str, help="Parameters of prior distributions (sep is | between parameters, and comma within parameters).")
    parser.add_argument("-m", "--mu", action="store", dest="mu", default=None, type=str, help="Means of prior distns, in the order parameter names are passed.")
    parser.add_argument("-sd", "--std", action="store", dest="std", default=None, type=str, help="Standard deviations of prior distns, in the order parameter names are passed.")
    parser.add_argument("-xt", "--xml-template", action="store", dest="xmlt", default=None, type=str, help="Full path to template of .xml file.")
    parser.add_argument("-pd", "--project-dir", action="store", dest="projdir", default=None, type=str, help="Full path to calibration folder if -pbs.")
    parser.add_argument("-ppf", "--prior-params-file", action="store", dest="prior_params_file", default=None, type=str, help="Full path to prior parameters file.")
    parser.add_argument("-e2t", "--event-to-triplet", action="store", dest="event_to_triplet", default=None, type=str, help="Full path to event to triplet mapping file.")

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

    if args.prior_params_file and args.event_to_triplet:
        param_names, prior_dists, prior_params = parse_param_event_files(args.prior_params_file, args.event_to_triplet)
        pnames = ",".join(param_names)
        pdists = ",".join(prior_dists)
        pparams = ";".join(prior_params)
    else:
        pnames = args.pnames
        pdists = args.prior_dists
        pparams = args.prior_params
        param_names = args.pnames.split(",")
        prior_dists = args.prior_dists.split(",")
        prior_params = args.prior_params.split(";")

    simulate(args.rscriptsdir, args.outputdir, args.nsims, args.bisse, args.prefix, args.simtime, pnames, pdists, pparams) # calls R script, produces .csv files and plots

    parse_simulations(args.outputdir, args.xmldir, args.xmlt, args.prefix, param_names, prior_dists, prior_params, args.projdir) # parses .csv into .xml files

    if args.projdir:
        if not os.path.exists("pbs_scripts"):
            os.makedirs("pbs_scripts")
        write_pbs(args.xmldir, args.prefix, args.projdir)
