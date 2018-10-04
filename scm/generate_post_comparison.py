import csv
import argparse
import matplotlib.pyplot as plt
import pylab

def read_csv_1(csv_name):
    with open(csv_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 1:
                post = [float(r) for r in row]
            line_count += 1
    return post

def read_csv_2(csv_name):
    with open(csv_name) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
            post = [float(r) for r in row]
    return post


def construct_post_comparison(csv_1_name, csv_2_name):
    csv_1 = read_csv_1(csv_1_name)
    csv_2 = read_csv_2(csv_2_name)

    plt.scatter(csv_1, csv_2)
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    pylab.savefig("asr.png")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Generate asr validation", description="Generates ancestral state reconstruction figures that validate method.")
    parser.add_argument("-c1", "--csv-1-name", action="store", dest="csv1", default=None, type=str, help="CSV file containing post data.")
    parser.add_argument("-c2", "--csv-2-name", action="store", dest="csv2", default=None, type=str, help="CSV file containing post data.")
    args = parser.parse_args()

    construct_post_comparison(args.csv1, args.csv2)
