from absl import flags
import subprocess
import csv
import sys
import os

FLAGS = flags.FLAGS
flags.DEFINE_string("save_dir", "sim_data", "Where to store data")
flags.DEFINE_string("run_name", "run", "Name for the run")
flags.DEFINE_integer("num_sims", 100, "Number of simulations")
flags.DEFINE_boolean("is_bisse", True, "True - Bisse, False - classe")
FLAGS(sys.argv)

if not os.path.exists(FLAGS.save_dir):
  os.makedirs(FLAGS.save_dir)

path="/Users/jeff/Documents/Research/Phylogenetics/calibrated_validation"
path_param = path + "/simulate_params.R"
path_sse = path + "/simulate_SSE.R"
call_sse = ["Rscript", "--vanilla", path_sse, FLAGS.save_dir, FLAGS.run_name]


if FLAGS.is_bisse:
    subprocess.call (["Rscript", "--vanilla", path_param, FLAGS.save_dir,
                  FLAGS.run_name, str(FLAGS.num_sims), "l1", "l2", "m1", "m2", "q12", "q21"])
else:
    param_names = ["l111", "l112", "l22", "l211", "l221", "l222", "m1",
                    "m2", "q12", "q21"]
    subprocess.call (["Rscript", "--vanilla", path_param, FLAGS.save_dir,
                  FLAGS.run_name, str(FLAGS.num_sims)] + param_names)


path = os.path.join(FLAGS.save_dir, FLAGS.run_name)
with open(path + '_all_params.csv') as csv_file:
  csv_reader = csv.reader(csv_file, delimiter=',')
  for i, row in enumerate(csv_reader):
    print(f'\t{row[0]} works in the {row[1]} department, and was born in {row[2]}.')
    call = call_sse + row
    call[4] += str(i)
    subprocess.call(call)

