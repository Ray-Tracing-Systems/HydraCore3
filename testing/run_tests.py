import sys
import os
import subprocess

os.environ["OPENCV_IO_ENABLE_OPENEXR"]="1"
import cv2

from logger import Log, Status
from colorama import Fore

PATH_TO_TESTS = "/home/frol/PROG/HydraRepos3/HydraAPI-tests"
ENABLE_GPU    = True
GLOBAL_SPP    = 1024

SKIP_NAIVE  = 1
SKIP_SHADOW = 2


############################################################################################################
############################################################################################################
############################################################################################################

def run_sample(test_name, imsize, skip, on_gpu=False, gpu_id=0):
    Log().info("  rendering scene: {0}, gpu={1}".format(test_name, on_gpu))
    full = PATH_TO_TESTS + "/tests_f/" + test_name + "/statex_00001.xml"
    outp = PATH_TO_TESTS + "/tests_images/" + test_name + "/z_out.bmp"
    args = ["./cmake-build-release/hydra", "-in", full, "-out", outp, "-integrator", "all", "-spp-naive-mul", str(16)]
    args = args + ["-gpu_id", str(gpu_id)]  # for single launch samples
    args = args + ["-width", str(imsize[0]), "-height", str(imsize[1])]
    if on_gpu:
        args = args + ["--gpu"]
    #print(args)
    try:
        res = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        image_ref    = cv2.imread(PATH_TO_TESTS + "/tests_images/" + test_name + "/w_ref.png",         cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
        image_naive  = cv2.imread(PATH_TO_TESTS + "/tests_images/" + test_name + "/z_out_naivept.bmp", cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
        image_shadow = cv2.imread(PATH_TO_TESTS + "/tests_images/" + test_name + "/z_out_shadowpt.bmp",cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
        image_mis    = cv2.imread(PATH_TO_TESTS + "/tests_images/" + test_name + "/z_out_mispt.bmp",   cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
        metrics      = [cv2.PSNR(image_ref,image_naive), cv2.PSNR(image_ref,image_shadow), cv2.PSNR(image_ref,image_mis)]
        if skip == SKIP_NAIVE:
          metrics = metrics[1:]
        elif skip == SKIP_SHADOW:
          metrics = [metrics[0], metrics[2]]
        minPSNR      = min(metrics) 
        color        = Fore.GREEN
        message      = "[PASSED]"
        if minPSNR < 35.0: (color,message) = (Fore.YELLOW,"[PASSED]") 
        if minPSNR < 30.0: (color,message) = (Fore.RED, "[FAILED]") 
        if skip != 0:
          methods = "shadow,mis" if skip == SKIP_NAIVE else "naive,mis"
          Log().print_colored_text("  {}: PSNR({}) = ({:.2f},{:.2f})".format(message,methods,metrics[0], metrics[1]), color = color)
        else:
          Log().print_colored_text("  {}: PSNR(naive,shadow,mis) = ({:.2f},{:.2f},{:.2f})".format(message,metrics[0], metrics[1], metrics[2]), color = color)
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
  def __init__(self, name, tests, imsize = [(512,512)], skip = 0):
    self.name   = name
    self.tests  = tests
    self.imsize = imsize
    self.skip   = skip

reqs = []
reqs.append( REQ("mat: lambert", ["test_101"]) )
reqs.append( REQ("mat: mirror",  ["test_102"], [(1024, 1024)], SKIP_SHADOW) )
reqs.append( REQ("mat: lambert_texture", ["test_103"]) )

Log().set_workdir(".")
Log().info("PATH_TO_TESTS = {}".format(PATH_TO_TESTS))

os.chdir('..') # use HydraCore3 root dir as current
Log().set_workdir("testing")

for req in reqs:
  Log().info("REQ: {}".format(req.name))
  for (test,imsize) in zip(req.tests, req.imsize):
    run_sample(test, imsize, req.skip, ENABLE_GPU)
