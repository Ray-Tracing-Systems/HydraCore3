import sys
import os
import subprocess

os.environ["OPENCV_IO_ENABLE_OPENEXR"]="1"
import cv2

from logger import Log, Status
from colorama import Fore

TEST_CPU             = False
PATH_TO_HYDRA2_TESTS = "/home/frol/PROG/HydraRepos/HydraAPI-tests"
PATH_TO_HYDRA3_SCENS = "/home/frol/PROG/gitflic/comparisonrender"

############################################################################################################
############################################################################################################
############################################################################################################

def alignIntegratorName(name): # "naivept","shadowpt","mispt" ==> "naivept ","shadowpt","mispt   " 
  if name == "mispt":
    return name + "   "
  elif  name == "naivept":
    return name + " "
  else:
    return name

############################################################################################################
############################################################################################################

class REQ:
  def __init__(self, name, tests, imsize = (512,512), inregrators = ["naivept","shadowpt","mispt"], naivemul = 4):
    self.name   = name
    self.tests  = tests
    self.imsize = imsize
    self.integs = inregrators
    self.naivem = naivemul
  
  def test(req, gpu_id=0):
    pass

class REQ_H2(REQ):
  def __init__(self, name, tests, imsize = (512,512), inregrators = ["naivept","shadowpt","mispt"], naivemul = 4):
    self.name   = name
    self.tests  = tests
    self.imsize = imsize
    self.integs = inregrators
    self.naivem = naivemul

  def test(req, gpu_id=0):
    for test_name in req.tests: 
      image_ref = cv2.imread(PATH_TO_HYDRA2_TESTS + "/tests_images/" + test_name + "/w_ref.png", cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
      full = PATH_TO_HYDRA2_TESTS + "/tests_f/" + test_name + "/statex_00001.xml"
      devices = ["gpu"] if not TEST_CPU else ["gpu", "cpu"]
      for dev_type in devices:
        Log().info("  rendering scene: '{0}', dev_type='{1}', dev_id = '{2}'".format(test_name, dev_type, gpu_id))
        for inregrator in req.integs:
          outp = PATH_TO_HYDRA2_TESTS + "/tests_images/" + test_name + "/z_" + dev_type + inregrator + ".bmp"
          args = ["./cmake-build-release/hydra", "-in", full, "-out", outp, "-integrator", inregrator, "-spp-naive-mul", str(req.naivem), "-gamma", "2.2"]
          args = args + ["-gpu_id", str(gpu_id)]  # for single launch samples
          args = args + ["-width", str(req.imsize[0]), "-height", str(req.imsize[1])]
          args = args + ["--" + dev_type]
          #print(args)
          try:
            res = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)      
            image_mis    = cv2.imread(outp, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
            PSNR         = cv2.PSNR(image_ref,image_mis)
            color        = Fore.GREEN
            message      = "[PASSED]"
            if PSNR < 35.0: (color,message) = (Fore.YELLOW,"[PASSED]") 
            if PSNR < 30.0: (color,message) = (Fore.RED, "[FAILED]") 
            Log().print_colored_text("  {}: PSNR({}) = {:.2f}".format(message,alignIntegratorName(inregrator),PSNR), color = color)
          except Exception as e:
            Log().status_info("Failed to launch sample {0} : {1}".format(test_name, e), status=Status.FAILED)
            return
          if res.returncode != 0:
            Log().status_info("{}: launch".format(test_name), status=Status.FAILED)
            Log().save_std_output(test_name, res.stdout.decode(), res.stderr.decode())
            return
      if TEST_CPU:
        Log().info("  compare CPU/GPU:")
        for inregrator in req.integs:
          out_gpu = PATH_TO_HYDRA2_TESTS + "/tests_images/" + test_name + "/z_gpu" + inregrator + ".bmp"
          out_cpu = PATH_TO_HYDRA2_TESTS + "/tests_images/" + test_name + "/z_cpu" + inregrator + ".bmp"
          image1  = cv2.imread(out_gpu, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
          image2  = cv2.imread(out_cpu, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
          PSNR    = cv2.PSNR(image1, image2)
          color   = Fore.GREEN
          message = "[PASSED]"
          if PSNR < 40.0: (color,message) = (Fore.YELLOW,"[PASSED]") 
          if PSNR < 35.0: (color,message) = (Fore.RED, "[FAILED]") 
          Log().print_colored_text("  {}: PSNR({}, CPU/GPU) = {:.2f}".format(message,alignIntegratorName(inregrator),PSNR), color = color) 

class REQ_H3(REQ):
  def __init__(self, name, scn_path, ref_path, imsize = (1024,1024), inregrators = ["naivept","shadowpt","mispt"], naivemul = 4):
    self.name   = name
    self.scn_path = scn_path
    self.ref_path = ref_path
    self.imsize = imsize
    self.integs = inregrators
    self.naivem = naivemul

  def test(req, gpu_id=0):
    image_ref  = cv2.imread(PATH_TO_HYDRA3_SCENS + "/Report/Images/" + req.ref_path, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
    scene_path = PATH_TO_HYDRA3_SCENS + "/" + req.scn_path
    devices = ["gpu"] if not TEST_CPU else ["gpu", "cpu"]
    for dev_type in devices:
      Log().info("  rendering scene: '{0}', dev_type='{1}', dev_id = '{2}'".format(req.name, dev_type, gpu_id))
      for inregrator in req.integs:
        outp = PATH_TO_HYDRA3_SCENS + "/Report/Images/" + req.name  + "/z_" + dev_type + inregrator + ".bmp"
        args = ["./cmake-build-release/hydra", "-in", scene_path, "-out", outp, "-integrator", inregrator, "-spp-naive-mul", str(req.naivem), "-gamma", "2.2"]
        args = args + ["-gpu_id", str(gpu_id)]  # for single launch samples
        args = args + ["-width", str(req.imsize[0]), "-height", str(req.imsize[1])]
        args = args + ["--" + dev_type]
        #print(args)
        try:
          res = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)      
          image_mis    = cv2.imread(outp, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
          PSNR         = cv2.PSNR(image_ref,image_mis)
          color        = Fore.GREEN
          message      = "[PASSED]"
          if PSNR < 35.0: (color,message) = (Fore.YELLOW,"[PASSED]") 
          if PSNR < 30.0: (color,message) = (Fore.RED, "[FAILED]") 
          Log().print_colored_text("  {}: PSNR({}) = {:.2f}".format(message,alignIntegratorName(inregrator),PSNR), color = color)
        except Exception as e:
          Log().status_info("Failed to launch sample {0} : {1}".format(test_name, e), status=Status.FAILED)
          return
        if res.returncode != 0:
          Log().status_info("{}: launch".format(test_name), status=Status.FAILED)
          Log().save_std_output(test_name, res.stdout.decode(), res.stderr.decode())
          return
  
reqs = []

reqs.append( REQ_H2("mat: lambert", ["test_101"], naivemul = 16) )
reqs.append( REQ_H2("mat: mirror",  ["test_102"], inregrators = ["naivept","mispt"]) )
reqs.append( REQ_H2("mat: lambert_texture", ["test_103"]) )

reqs.append( REQ_H3("SmoothPlastic", "smooth_plastic/SmoothPlastic_sphere_hydra.xml",  "SmoothPlastic/SmoothPlastic_sphere_mitsuba.png",  naivemul = 16) )
reqs.append( REQ_H3("SmoothPlastic", "smooth_plastic/SmoothPlastic_cornell_hydra.xml", "SmoothPlastic/SmoothPlastic_cornell_mitsuba.png", naivemul = 16) )

Log().set_workdir(".")
Log().info("PATH_TO_TESTS = {}".format(PATH_TO_HYDRA2_TESTS))

os.chdir('..') # use HydraCore3 root dir as current
Log().set_workdir("testing")

for req in reqs:
  Log().info("REQ: {}".format(req.name))
  req.test()


