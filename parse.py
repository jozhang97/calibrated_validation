import argparse
import csv

class Parser:
    def __init__(self, params, tree, n_tips, tip_states):
        self.params = params
        self.tree = tree
        self.n_tips = n_tips
        self.tip_states = tip_states

    def load_to_xml(self):
        pass

    def __str__(self):
        return str(self.params) + str(self.tree)

def parse_csv(file_name):
    with open(file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for i, row in enumerate(csv_reader):
            params = row[:-3]
            tree = row[-3]
            n_tips = row[-2]
            tip_states = row[-1]

            p = Parser(params, tree, n_tips, tip_states)
            print(p)


if __name__ == "__main__":
    parse_csv("/Users/jeff/Documents/Research/Phylogenetics/calibrated_validation/sim_data_3/bisse_data.csv")
    pass
    parser = argparse.ArgumentParser(prog="Simulation script", description="Script for performing well-calibrated validation of biogeo package.")
    parser.add_argument("-rd", "--input-dir", action="store", dest="rscriptsdir", default="./", type=str, help="full file path to r scripts directory.")
    parser.add_argument("-od", "--output-dir", action="store", dest="outputdir", default="./", type=str, help="full file path to directory where simulations and log files will be saved.")
    parser.add_argument("-n", "--n-sims", action="store", dest="nsims", default=1000, type=int, help="number of simulations.")
    parser.add_argument("-b", "--is-bisse", action="store", dest="bisse", default=True, type=bool, help="Flag for BiSSE simulations (default: True)")
    parser.add_argument("-p", "--prefix", action="store", dest="prefix", default="", type=str, help="Prefix for result files.")
    parser.add_argument("-st", "--sim-time", action="store", dest="simtime", default=10, type=float, help="Time to run simulation.")
    parser.add_argument("-m", "--mu", action="store", dest="mu", default=0, type=float, help="Mean of lognormal dist.")
    parser.add_argument("-sd", "--std", action="store", dest="std", default=0.05, type=float, help="Stddev of lognormal dist.")
    args = parser.parse_args()
