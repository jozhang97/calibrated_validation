import argparse
import csv

class Parser:
    def __init__(self, params, tree, n_tips, tip_states):
        self.params = params
        self.tree = tree
        self.n_tips = n_tips
        self.tip_states = tip_states

    def populate_xml(self):
        # TODO
        pass

    def write_xml(self, output_dir):
        # TODO
        pass

def parse_csv(file_name, output_dir):
    with open(file_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for i, row in enumerate(csv_reader):
            params = row[:-3]
            tree = row[-3]
            n_tips = row[-2]
            tip_states = row[-1]

            p = Parser(params, tree, n_tips, tip_states)


if __name__ == "__main__":
    parse_csv("/Users/jeff/Documents/Research/Phylogenetics/calibrated_validation/sim_data_3/bisse_data.csv")
    pass
    parser = argparse.ArgumentParser(prog="CSV parser", description="Script for parsing simulation data to BEAST XML.")
    parser.add_argument("-cp", "--csv-path", action="store", dest="csvpath", default="./", type=str, help="full file path to csv.")
    parser.add_argument("-od", "--output-dir", action="store", dest="outputdir", default="./", type=str, help="full file path to directory where XML files will be saved.")
    args = parser.parse_args()

    parse_csv(args.csvpath, args.outputdir)
