import sys
import os
import subprocess

os.environ["OPENCV_IO_ENABLE_OPENEXR"]="1"
import cv2

from logger import Log, Status
from colorama import Fore

TEST_CPU             = False
PATH_TO_HYDRA2_TESTS = "/home/frol/PROG/HydraRepos/HydraAPI-tests"
PATH_TO_HYDRA3_SCENS = "/home/frol/PROG/HydraRepos/comparisonrender"

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
  
  def run(req, test_name, args, image_ref, outp, inregrator):
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
        Log().info("  rendering scene: '{0}', dev_type='{1}', scene = '{2}'".format(test_name, dev_type, full))
        for inregrator in req.integs:
          outp = PATH_TO_HYDRA2_TESTS + "/tests_images/" + test_name + "/z_" + dev_type + inregrator + ".bmp"
          args = ["./cmake-build-release/hydra", "-in", full, "-out", outp, "-integrator", inregrator, "-spp-naive-mul", str(req.naivem), "-gamma", "2.2"]
          args = args + ["-gpu_id", str(gpu_id)]  # for single launch samples
          args = args + ["-width", str(req.imsize[0]), "-height", str(req.imsize[1])]
          args = args + ["--" + dev_type]
          req.run(test_name, args, image_ref, outp, inregrator)
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

class REQ_HX(REQ):
  def __init__(self, name, scn_path, ref_path, imsize = [(1024,1024)], inregrators = ["naivept","shadowpt","mispt"], naivemul = 4):
    self.name   = name
    self.scn_path = scn_path
    self.ref_path = ref_path
    self.imsize = imsize
    self.integs = inregrators
    self.naivem = naivemul

  def test(req, gpu_id=0):
    for (scnp, imgp, id, imsize2) in zip(req.scn_path, req.ref_path, range(0,len(req.scn_path)), req.imsize):
      image_ref   = cv2.imread(imgp, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
      scene_path  = scnp
      folder_path = os.path.dirname(imgp)
      devices   = ["gpu"] if not TEST_CPU else ["gpu", "cpu"]
      test_name = req.name
      for dev_type in devices:
        Log().info("  rendering scene: '{0}', dev_type='{1}', scene = '{2}'".format(test_name, dev_type, scene_path))
        for inregrator in req.integs:
          outp = folder_path + "/y" + str(id) + "_" + dev_type + inregrator + ".bmp"
          args = ["./cmake-build-release/hydra", "-in", scene_path, "-out", outp, "-integrator", inregrator, "-spp-naive-mul", str(req.naivem), "-gamma", "2.2"]
          args = args + ["-gpu_id", str(gpu_id)]  # for single launch samples
          args = args + ["-width", str(imsize2[0]), "-height", str(imsize2[1])]
          args = args + ["--" + dev_type]
          #print(args)
          req.run(test_name, args, image_ref, outp, inregrator)
  
reqs = []


reqs.append( REQ_HX("geo_inst_remap_list", [PATH_TO_HYDRA2_TESTS + "/tests/test_078/statex_00001.xml", 
                                            PATH_TO_HYDRA2_TESTS + "/tests/test_078/statex_00002.xml", 
                                            PATH_TO_HYDRA2_TESTS + "/tests/test_079/statex_00001.xml", 
                                            PATH_TO_HYDRA2_TESTS + "/tests/test_079/statex_00002.xml"],

                                           [PATH_TO_HYDRA2_TESTS + "/tests_images/test_078/w_ref.png", 
                                            PATH_TO_HYDRA2_TESTS + "/tests_images/test_078/w_ref2.png",
                                            PATH_TO_HYDRA2_TESTS + "/tests_images/test_079/w_ref.png", 
                                            PATH_TO_HYDRA2_TESTS + "/tests_images/test_079/w_ref2.png"], 
                                            imsize = [(512,512), (512,512), (512,512), (512,512)], naivemul = 1))


reqs.append( REQ_HX("mat_lambert", [PATH_TO_HYDRA2_TESTS + "/tests_f/test_101/statex_00001.xml",  
                                    PATH_TO_HYDRA3_SCENS + "/Lambert/Lambert_cornell_hydra.xml"],
                                   [PATH_TO_HYDRA2_TESTS + "/tests_images/test_101/w_ref.png", 
                                    PATH_TO_HYDRA3_SCENS + "/Report/Images/Lambert/Lambert_cornell_mitsuba.png"], 
                                    imsize = [(512,512), (1024,1024)], naivemul = 16))

reqs.append( REQ_HX("mat_emission", [PATH_TO_HYDRA2_TESTS + "/tests_f/test_123/statex_00001.xml",  
                                     PATH_TO_HYDRA2_TESTS + "/tests_f/test_123/statex_00002.xml",
                                     PATH_TO_HYDRA2_TESTS + "/tests_f/test_123/statex_00003.xml"],
                                    [PATH_TO_HYDRA2_TESTS + "/tests_images/test_123/w_ref.png", 
                                     PATH_TO_HYDRA2_TESTS + "/tests_images/test_123/w_ref2.png",
                                     PATH_TO_HYDRA2_TESTS + "/tests_images/test_123/w_ref3.png"], 
                                     imsize = [(512,512), (512,512), (512,512)],
                                     naivemul = 16, inregrators = ["naivept","mispt"]))

reqs.append( REQ_H2("mat_mirror",  ["test_102"], inregrators = ["naivept","mispt"]) )

reqs.append( REQ_HX("mat_smooth_plastic", [PATH_TO_HYDRA3_SCENS + "/smooth_plastic/SmoothPlastic_sphere_hydra.xml",  
                                           PATH_TO_HYDRA3_SCENS + "/smooth_plastic/SmoothPlastic_cornell_hydra.xml"],
                                          [PATH_TO_HYDRA3_SCENS + "/Report/Images/SmoothPlastic/SmoothPlastic_sphere_mitsuba.png", 
                                           PATH_TO_HYDRA3_SCENS + "/Report/Images/SmoothPlastic/SmoothPlastic_cornell_mitsuba.png"], naivemul = 16))
                                           
reqs.append( REQ_H2("mat_lambert_texture", ["test_103"]) )


Log().set_workdir(".")
Log().info("PATH_TO_TESTS = {}".format(PATH_TO_HYDRA2_TESTS))

os.chdir('..') # use HydraCore3 root dir as current
Log().set_workdir("testing")

for req in reqs:
  Log().info("REQ: {}".format(req.name))
  req.test()


