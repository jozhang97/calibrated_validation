from absl import flags
import subprocess
import csv
import sys
import os

FLAGS = flags.FLAGS
flags.DEFINE_string("save_dir", "sim_data", "Where to store data")
flags.DEFINE_string("run_name", "run1", "Name for the run")
flags.DEFINE_integer("num_sims", 100, "Number of simulations")
FLAGS(sys.argv)

path="/Users/jeff/Documents/Research/Phylogenetics/calibrated_validation"
path_param = path + "/simulate_params.R"
path_sse = path + "/simulate_SSE.R"


subprocess.call (["Rscript", "--vanilla", path_param, FLAGS.save_dir,
                  FLAGS.run_name, str(FLAGS.num_sims), "l1", "l2", "m1", "m2", "q12", "q21"])

path = os.path.join(FLAGS.save_dir, FLAGS.run_name)
with open(path + '_params.csv') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=',')
  import ipdb; ipdb.set_trace()
  line_count = 0
  for row in csv_reader:
    if line_count == 0:
      line_count += 1
    else:
      print(f'\t{row[0]} works in the {row[1]} department, and was born in {row[2]}.')
      line_count += 1
      print(f'Processed {line_count} lines.')



subprocess.call (["Rscript", "--vanilla", path_sse, "sim_data",
                 "run1", ".2", ".4", ".01", ".1", ".1", ".4"])
