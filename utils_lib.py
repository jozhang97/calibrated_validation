import re
import os

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

def stringfy_prior_params(prior_dist, prior_params, param_name, quoted=False):
    """ Generate prior distribution parameter string to put in BEAST .xml """

    # prior_dist = e.g. ["lnorm", "lnorm"]
    # prior_params = e.g. ["0.1,0.05" "0.1,0.05"] they will always be repeated
    prior_param_string = str()
    prior_dist_str = str()
    prior_params_list = list()

    if type(prior_params) is list:
        prior_params_list = prior_params[0].split(",")
        prior_dist_str = prior_dist[0].split(",")[0]
    else:
        prior_params_list = [prior_params]
        prior_dist_str = prior_dist
        
    # lognormal prior
    if prior_dist_str == "lnorm":
        mu = prior_params_list[0]
        std = prior_params_list[1]

        if quoted: mu = "\"" + mu + "\""; std = "\"" + std + "\""

        prior_param_string = "<distr id=\"LogNormal." + param_name + "\" spec=\"beast.math.distributions.LogNormalDistributionModel\" offset=\"0.0\" meanInRealSpace=\"false\">\n"
        prior_param_string += "\t\t<parameter name=\"M\" value=" + mu + " estimate=\"false\"/>\n"
        prior_param_string += "\t\t<parameter name=\"S\" value=" + std + " estimate=\"false\"/>\n\t  </distr>\n"

    # exponential prior
    elif prior_dist_str == "exp":
        rate = 1/float(prior_params_list[0])

        if quoted: rate = "\"" + str(rate) + "\""
        
        prior_param_string = "\t<distr id=\"Exponential." + param_name + "\" spec=\"beast.math.distributions.Exponential\" offset=\"0.0\" mean=" + rate + "/>\n"

    return prior_param_string

def stringfy_triplets(parent_state_list, left_state_list, right_state_list, event_type_list):
    """ Generate triplet strings to put in BEAST .xml when using CLaSSE """

    triplet_string = str()
    for i, event_type in enumerate(event_type_list):
        if event_type == "S": event_type = "SYMPATRY"
        if event_type == "SS": event_type = "SUBSYMPATRY"
        if event_type == "V": event_type = "VICARIANCE"
        ith_triplet_string = "<CladoTriplets id=\"CladoTriplet" + str(i+1) + "\" spec=\"biogeo.CladoTriplet\" LeftChildState=\"" + left_state_list[i] + "\" ParentState=\"" + parent_state_list[i] + "\" RightChildState=\"" + right_state_list[i] + "\" SpeciationType=\"" + event_type + "\"/>\n\t    "
        triplet_string += ith_triplet_string

    return triplet_string.rstrip('\t    \n')

def write_pbs(xml_dir, prefix, project_dir):
    """ Write one .pbs script per .xml """

    sim_n_regex = re.compile("[0-9]+")
    xml_file_names = [f for f in os.listdir(xml_dir) if f.endswith(".xml")]
    for xml_file_name in xml_file_names:
        sim_n = re.findall(sim_n_regex, xml_file_name.split("_")[0])[0]

        with open(prefix + "_pbs_scripts/" + prefix + sim_n + ".PBS", "w") as pbs_file:
            pbs_file.write("#!/bin/bash\n#PBS -N beast_" + sim_n + \
                           "\n#PBS -l nodes=1:ppn=1,walltime=70:00:00\n#PBS -M fkmendes@iu.edu\n#PBS -m abe\n\njava -jar " + project_dir + "biogeo.jar " + \
                           project_dir + xml_dir + xml_file_name
            )

def natural_sort(l): 
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 

    return sorted(l, key = alphanum_key)

def find_matching_index(sorted_list, unsorted_list):
    unsorted_idxs = { element: index for index, element in enumerate(unsorted_list) }

    return [(sorted_idx, unsorted_idxs[element])
        for sorted_idx, element in enumerate(sorted_list)]

def write_sbatch(xml_dir, prefix, project_dir):
    """ Write one .sh script per .xml for running on NeSI """

    sim_n_regex = re.compile("[0-9]+")
    xml_file_names = [f for f in os.listdir(xml_dir) if f.endswith(".xml")]
    for xml_file_name in xml_file_names:
        sim_n = re.findall(sim_n_regex, xml_file_name.split("_")[0])[0]

        with open(prefix + "_shell_scripts/" + prefix + sim_n + ".sh", "w") as shell_file:
            shell_file.write("#!/bin/bash -e\n#SBATCH -J beast_" + sim_n + "\n" +\
                             "#SBATCH -A nesi00390\n" + \
                             "#SBATCH --time=70:00:00\n" + \
                             "#SBATCH --mem-per-cpu=12288\n" + \
                             "#SBATCH --cpus-per-task=1\n" + \
                             "#SBATCH --ntasks=1\n" + \
                             "#SBATCH --hint=nomultithread\n" + \
                             "#SBATCH -D ./\n" + \
                             "#SBATCH -o beast_" + sim_n + "_out.txt\n" + \
                             "#SBATCH -e beast_" + sim_n + "_err.txt\n\n" + \

                             "srun /nesi/project/nesi00390/fkmendes/programs/jdk-10.0.2/bin/java -jar " + project_dir + "biogeo.jar " + \
                             project_dir + xml_dir + xml_file_name
            )
