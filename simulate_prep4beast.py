import sys
import os
import subprocess
import csv
import argparse
import numpy as np
from utils_lib import *

def parse_param_event_files(prior_params_file_name, event_to_triplet):
    """
    Parse the prior params and event triplet files to supply simulate() with the data in req data format.
    Will associate the correct prior parameters (depending on speciation event type) to each triplet, as specified in event triplet file.
    """

    param_names, prior_dists, prior_params = list(), list(), list()

    # preparing mapping
    e2t_dict = dict()
    with open(event_to_triplet) as e2t_file:
        e2t = csv.reader(e2t_file, delimiter='|')

        for i, event_triplet in enumerate(e2t):
            if i == 0: continue # ignore header
            event_name, triplet = event_triplet
            triplet = "l_" + "".join(triplet.split(","))

            try:
                e2t_dict[event_name].append(triplet)
            except:
                e2t_dict[event_name] = [triplet]
    
    # reading file of param names and priors
    with open(prior_params_file_name) as prior_params_file:
        prior_params_line = csv.reader(prior_params_file, delimiter='|')

        for i, param_info in enumerate(prior_params_line):
            if i == 0: continue # ignore header
            param_name, dist_type, moments = param_info

            # death and transition parameters
            if param_name.startswith("m") or param_name.startswith("q"):
                param_names.append(param_name)
                prior_dists.append(dist_type)
                prior_params.append(moments)

            else:
                if param_name in e2t_dict:
                    for l in e2t_dict[param_name]:
                        param_names.append(l)
                        prior_dists.append(dist_type)
                        prior_params.append(moments)

    sorted_param_names = natural_sort(param_names)
    list_of_idx_tuples = find_matching_index(sorted_param_names, param_names) # (sorted_idx, unsorted_idx)
    sorted_prior_dists = [prior_dists[unsorted_idx] for sorted_idx, unsorted_idx in list_of_idx_tuples]
    sorted_prior_params = [prior_params[unsorted_idx] for sorted_idx, unsorted_idx in list_of_idx_tuples]
    # print sorted_param_names
    # print sorted_prior_params
    return sorted_param_names, sorted_prior_dists, sorted_prior_params, e2t_dict

def convert_to_flat_rate_matrix(q):
    # solve n * n - n = len(q)
    if len(q) == 2:  # BiSSE
        return " ".join([q[0], q[1]])
    raise NotImplementedError
    #roots = np.roots(1, -1, len(q))
    #n = max(roots)  # not sure
    # TODO ... this is kinda annoying, keep looking for a library

class CsvInfoStash:
    def __init__(self, init_params, param_names, prior_dists, prior_params, tree, n_tips, tip_states, output_dir, prefix, idx, is_bisse, prior_params_file_name = None, e2t_dict = None):
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
        self.is_bisse = is_bisse
        self.prior_params_file_name = prior_params_file_name
        self.event_triplet_dict = e2t_dict

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
        self.replace_in_xml("[Tree in newick format]", self.tree)
        self.replace_in_xml("[List of species]", convert_to_species_list(self.tip_states))
        self.replace_in_xml("[Species=Trait State]", self.tip_states)

        if (self.is_bisse):
            self.replace_in_xml("[Lambda Prior Parameters]", stringfy_prior_params(self.prior_dists[:2], self.prior_params[:2], "Lambda", quoted=True)) # note that quoted applies to stringfy function
            self.replace_in_xml("[Mu Prior Parameters]", stringfy_prior_params(self.prior_dists[3:5], self.prior_params[3:5], "Mu", quoted=True))
            self.replace_in_xml("[FlatQMatrix Prior Parameters]", stringfy_prior_params(self.prior_dists[4:], self.prior_params[4:], "FlatQMatrix", quoted=True))
            self.replace_in_xml("[Initial lambda values]", self.init_params[0:2])
            self.replace_in_xml("[Initial mu values]", self.init_params[2:4])
            self.replace_in_xml("[Initial transition rate values]", convert_to_flat_rate_matrix(self.init_params[4:6]))
            self.replace_in_xml("[Pi values]", "0.0 0.0 0.5 0.5")

        else:
            # replacing prior on death rates
            m_params_idxs = [idx for idx, p in enumerate(self.param_names) if p.startswith("m")]
            m_dists = [self.prior_dists[i] for i in m_params_idxs]
            m_params = [self.prior_params[i] for i in m_params_idxs]
            self.replace_in_xml("[Mu Prior Parameters]", stringfy_prior_params(m_dists, m_params, "Mu", quoted=True))

            # start: replacing prior on transition rates # 
            q_params_idxs = [idx for idx, p in enumerate(self.param_names) if p.startswith("q")]
            q_dists = [self.prior_dists[i] for i in q_params_idxs]
            q_params = [self.prior_params[i] for i in q_params_idxs]
            self.replace_in_xml("[FlatQMatrix Prior Parameters]", stringfy_prior_params(m_dists, m_params, "FlatQMatrix", quoted=True))
            # end: replacing prior on transition rates # 

            # start: dealing with info on starting values of transition rates, and pi
            n_states = len(m_params_idxs)
            self.replace_in_xml("[Number of transition rates minus diagonals]", n_states*n_states-n_states)
            self.replace_in_xml("[Number of states]", n_states)
            self.replace_in_xml("[Twice the number of states]", n_states*2)
            self.replace_in_xml("[Pi values]", " ".join([str(0) for i in range(n_states)] + [str(float(1)/n_states) for i in range(n_states)]))
            # end: dealing with info on starting values of transition rates

            # start: replacing prior on lambdas in template #
            event_dist_dict = dict()
            with open(self.prior_params_file_name, "r") as prior_params_file:
                for idx, line in enumerate(prior_params_file):
                    if idx == 0: continue
                    if line.startswith("m") or line.startswith("q") or line.startswith("IG"): continue
                    event, dist, params = line.rstrip().split("|")
                    event_dist_dict[event] = [dist, params]

            self.replace_in_xml("[Sympatric Prior Parameters]", stringfy_prior_params(event_dist_dict["S"][0], event_dist_dict["S"][1], "SympatricRate", quoted=True)) # note that quoted applies to stringfy function
            self.replace_in_xml("[Subsympatric Prior Parameters]", stringfy_prior_params(event_dist_dict["SS"][0], event_dist_dict["SS"][1], "SubsympatricRate", quoted=True))
            self.replace_in_xml("[Vicariant Prior Parameters]", stringfy_prior_params(event_dist_dict["S"][0], event_dist_dict["V"][1], "VicariantRate", quoted=True))
            # end: replacing prior on lambdas in template #

            # start: filling out triplets #
            event_type_list, parent_state_list, left_state_list, right_state_list = list(), list(), list(), list()
            event_matching_triplet_for_init = dict()
            for event, triplet_list in self.event_triplet_dict.items():
                if event == "IGNORE": continue

                for triplet_str in triplet_list:
                    event_type_list.append(event)

                    # filling out dict that will be used for initialization values in .xml
                    if not event in event_matching_triplet_for_init:
                        event_matching_triplet_for_init[event] = triplet_list[0]
                        
                    this_triplet_list = list(triplet_str.split("_")[1]) # does not work if there are more than 9 states!
                    parent_state_list.append(this_triplet_list[0])
                    left_state_list.append(this_triplet_list[1])
                    right_state_list.append(this_triplet_list[2])
            
            self.replace_in_xml("[Cladogenetic triplets]", stringfy_triplets(parent_state_list, left_state_list, right_state_list, event_type_list))
            # end: filling out triplets # 

            # start: filling out initialization lambda values #
            init_csv_tokens = list()
            with open(self.output_dir + "data_param_inits.csv", "r") as init_file:
                for idx, line in enumerate(init_file):
                    if idx == 0: # header is what we want
                        init_csv_tokens = line.rstrip().split("|")

            for event, info in event_matching_triplet_for_init.items():
                idx = init_csv_tokens.index(info) # which col of init params to grab
                init_value = self.init_params[idx]

                if event == "S":
                    self.replace_in_xml("[Initial sympatric rate value]", self.init_params[idx])
                                        
                elif event == "SS":
                    self.replace_in_xml("[Initial subsympatric rate value]", self.init_params[idx])
                    
                elif event == "V":
                    self.replace_in_xml("[Initial vicariant rate value]", self.init_params[idx])
            # end: filling out initialization values #

            # start: filling out initialization mu and q-matrix values #
            m_init_string = str()
            q_init_string = str()
            for idx, t in enumerate(init_csv_tokens):
                if t.startswith("m"):
                    m_init_string += self.init_params[idx]+" "

                elif t.startswith("q"):
                    q_init_string += self.init_params[idx]+" "

            m_init_string.rstrip(" "); q_init_string.rstrip(" ")

            self.replace_in_xml("[Initial mu values]", m_init_string)
            self.replace_in_xml("[Initial transition rate values]", q_init_string)
            # end: filling out initialization mu and pi values #
            
        beast_output_path = self.prefix + "_beast_outputs/"
        if self.project_dir:
            beast_output_path = project_dir + self.prefix + "_beast_outputs/"
        self.replace_in_xml("[Simulation log file name]", beast_output_path + self.prefix + str(self.idx) + "_sim.log")
        self.replace_in_xml("[Simulation tree file name]", beast_output_path + self.prefix + str(self.idx) + "_sim.trees")
        
    def write_xml(self, file_name):
        with open(file_name, "w") as f:
            f.write(self.xml)

def simulate(r_script_dir, output_dir, n_sims, prefix, sim_time, pnames, prior_dists, prior_params):
    """ Call simulation pipeline (create .csvs and plots) """

    output_path = os.path.join(output_dir, prefix)
    param_names = list()
    r_script_name = "simulate_params.R"

    r_script = r_script_dir + r_script_name
    cmd = ["Rscript", "--vanilla", str(r_script), output_dir, prefix, str(n_sims), pnames, prior_dists, prior_params, str(sim_time)]
    print "\nR command call: " + " ".join(cmd) + "\n"
    subprocess.call(cmd)

def parse_simulations(output_dir, xml_dir, xml_template_name, prefix, param_names, prior_dists, prior_params, project_dir, is_bisse, e2t_dict):
    """ Parse .csv file into .xmls """

    num_params = len(param_names)
    csv_info_stashes = list()
    csv_true = output_dir + "data_param_tree.csv"
    csv_inits = output_dir + "data_param_inits.csv"
    csv_classe_prior_params = str()
    if not is_bisse:
        csv_classe_prior_params = output_dir + "classe_prior_params.csv"
        
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

            csv_info_stash = CsvInfoStash(init_params, param_names, prior_dists, prior_params, tree, n_tips, tip_states, output_dir, prefix, i, is_bisse, csv_classe_prior_params, e2t_dict)
            csv_info_stashes.append(csv_info_stash)

    with open(xml_template_name, 'r') as xml_template_file:
        xml_template = xml_template_file.read()

    # iterating over list of csv_info_stashes
    for i, csv_info_stash in enumerate(csv_info_stashes):
        csv_info_stash.populate_xml(xml_template, project_dir) # fill out xml template
        xml_file_name = xml_dir + prefix + str(i+1) + ".xml"
        csv_info_stash.write_xml(xml_file_name)

    return

if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Simulation script", description="Script for performing well-calibrated validation of biogeo package.")
    parser.add_argument("-rd", "--rscripts-dir", action="store", dest="rscriptsdir", default="./", type=str, help="Full file path to r scripts directory.")
    parser.add_argument("-od", "--output-dir", action="store", dest="outputdir", default="./", type=str, help="Full file path to directory where simulations and log files will be saved.")
    parser.add_argument("-xd", "--xml-dir", action="store", dest="xmldir", default="./", type=str, help="Full file path to directory where .xml's will be saved.")
    parser.add_argument("-n", "--n-sims", action="store", dest="nsims", default=1000, type=int, help="Number of simulations.")
    parser.add_argument("-b", "--is-bisse", action="store_true", dest="bisse", help="Flag for BiSSE simulations.")
    parser.add_argument("-p", "--prefix", action="store", dest="prefix", default="", type=str, help="Prefix for result files.")
    parser.add_argument("-st", "--sim-time", action="store", dest="simtime", default=10, type=float, help="Time to run simulation.")
    parser.add_argument("-pn", "--param-names", action="store", dest="pnames", default=None, type=str, help="Parameter names (order matters for means and stdevs of priors).")
    parser.add_argument("-pt", "--prior-distributions", action="store", dest="prior_dists", default=None, type=str, help="Prior distributions (one per parameter).")
    parser.add_argument("-pp", "--prior-params", action="store", dest="prior_params", default=None, type=str, help="Parameters of prior distributions (sep is ; between parameters, and comma within parameters).")
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

    # CLaSSE
    e2t_dict = dict()
    if args.prior_params_file and args.event_to_triplet:
        param_names, prior_dists, prior_params, e2t_dict = parse_param_event_files(args.prior_params_file, args.event_to_triplet)
        pnames = ",".join(param_names)
        pdists = ",".join(prior_dists)
        pparams = ";".join(prior_params)
        print pnames, pdists, pparams

    # BiSSE
    else:
        pnames = args.pnames
        pdists = args.prior_dists
        pparams = args.prior_params
        param_names = args.pnames.split(",")
        prior_dists = args.prior_dists.split(",")
        prior_params = args.prior_params.split(";")

    # simulate(args.rscriptsdir, args.outputdir, args.nsims, args.prefix, args.simtime, pnames, pdists, pparams) # calls R script, produces .csv files and plots

    parse_simulations(args.outputdir, args.xmldir, args.xmlt, args.prefix, param_names, prior_dists, prior_params, args.projdir, args.bisse, e2t_dict) # parses .csv into .xml files

    if args.projdir:
        if not os.path.exists(args.prefix + "_pbs_scripts"):
            os.makedirs(args.prefix + "_pbs_scripts")
        write_pbs(args.xmldir, args.prefix, args.projdir)

        if not os.path.exists("shell_scripts"):
            os.makedirs("shell_scripts")
        write_sbatch(args.xmldir, args.prefix, args.projdir)
