import sys
import os
import subprocess

from logger import Log, Status

PATH_TO_TESTS = "/home/frol/PROG/HydraRepos2/HydraAPI-tests"
ENABLE_GPU    = True

############################################################################################################
############################################################################################################

def run_sample(test_name, on_gpu=False, gpu_id=0):
    Log().info("Running sample: {0}, gpu={1}".format(test_name, on_gpu))
    full = PATH_TO_TESTS + "/" + test_name + "/statex_00001.xml"
    outp = "z_out.bmp"
    args = ["./cmake-build-release/hydra", "-in", full, "-out", outp, "-method", "mispt"]
    args = args + ["--gpu_id", str(gpu_id)]  # for single launch samples
    if on_gpu:
        args = args + ["--gpu"]

    try:
        res = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception as e:
        Log().status_info("Failed to launch sample {0} : {1}".format(test_name, e), status=Status.FAILED)
        return -1
    if res.returncode != 0:
        Log().status_info("{}: launch".format(test_name), status=Status.FAILED)
        Log().save_std_output(test_name, res.stdout.decode(), res.stderr.decode())
        return -1

    return 0

############################################################################################################
############################################################################################################

class REQ:
  def __init__(self, name, tests):
    self.name  = name
    self.tests = tests

reqs = []
reqs.append( REQ("mat: lambert", ["tests_f/test_101"]) )
reqs.append( REQ("mat: lambert_texture", ["tests_f/test_103"]) )

Log().set_workdir(".")
Log().info("PATH_TO_TESTS = {}".format(PATH_TO_TESTS))

os.chdir('..') # use HydraCore3 root dir as current
for req in reqs:
  Log().info("REQ: {}".format(req.name))
  for test in req.tests:
    run_sample(test, ENABLE_GPU)
