import sys
import os
import subprocess
import re

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
  def __init__(self, name, tests, imsize = (512,512), integrators = ["naivept","shadowpt","mispt"], naivemul = 4):
    self.name   = name
    self.tests  = tests
    self.imsize = imsize
    self.integs = integrators
    self.naivem = naivemul
    self.times  = []
  
  def test(req, gpu_id=0):
    pass

  def need_render_time(req):
    return False
  
  def run(req, test_name, args, image_ref, outp, inregrator):
    try:
      res = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
      output  = res.stdout.decode('utf-8')
      #print(output)
      pattern = r'PathTraceBlock\(MIS-PT\) = (\d+\.\d+) ms'
      match   = re.search(pattern, output) 
      if match:
        execution_time_ms = round(float(match.group(1)))
        req.times.append(execution_time_ms)
        if req.need_render_time():
          print(f"  render-time: {execution_time_ms} ms")
      image_mis = cv2.imread(outp, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
      PSNR      = cv2.PSNR(image_ref,image_mis)
      color     = Fore.GREEN
      message   = "[PASSED]"
      if PSNR < 35.0: (color,message) = (Fore.YELLOW,"[PASSED]") 
      if PSNR < 30.0: (color,message) = (Fore.RED, "[FAILED]") 
      Log().print_colored_text("  {}: PSNR({}) = {:.2f}".format(message,alignIntegratorName(inregrator),PSNR), color = color)
    except Exception as e:
      Log().status_info("Failed to launch sample {0} : {1}".format(test_name, e), status=Status.FAILED)
      return
    if res.returncode != 0:
      Log().status_info("{}: launch, returncode = {}".format(test_name,res.returncode), status=Status.FAILED)
      Log().save_std_output(test_name, res.stdout.decode(), res.stderr.decode())
      return

class REQ_H2(REQ):
  def __init__(self, name, tests, imsize = (512,512), integrators = ["naivept","shadowpt","mispt"], naivemul = 4):
    self.name   = name
    self.tests  = tests
    self.imsize = imsize
    self.integs = integrators
    self.naivem = naivemul
    self.times  = []

  def test(req, gpu_id=0):
    for test_name in req.tests: 
      image_ref = cv2.imread(PATH_TO_HYDRA2_TESTS + "/tests_images/" + test_name + "/w_ref.png", cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
      full = PATH_TO_HYDRA2_TESTS + "/tests_f/" + test_name + "/statex_00001.xml"
      devices = ["gpu"] if not TEST_CPU else ["gpu", "cpu"]
      for dev_type in devices:
        Log().info("  rendering scene: '{0}', dev_type='{1}', scene = '{2}'".format(test_name, dev_type, full))
        for inregrator in req.integs:
          outp = PATH_TO_HYDRA2_TESTS + "/tests_images/" + test_name + "/z_" + dev_type + inregrator + ".bmp"
          args = ["./bin-release/hydra", "-in", full, "-out", outp, "-integrator", inregrator, "-spp-naive-mul", str(req.naivem)]
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
  def __init__(self, name, scn_path, ref_path, imsize = [(1024,1024)], integrators = ["naivept","shadowpt","mispt"], naivemul = 4, is_spectral = False):
    self.name   = name
    self.scn_path = scn_path
    self.ref_path = ref_path
    self.imsize = imsize
    self.integs = integrators
    self.naivem = naivemul
    self.spectral = is_spectral
    self.times  = []

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
          args = ["./bin-release/hydra", "-in", scene_path, "-out", outp, "-integrator", inregrator, "-spp-naive-mul", str(req.naivem)]
          args = args + ["-gpu_id", str(gpu_id)]  # for single launch samples
          args = args + ["-width", str(imsize2[0]), "-height", str(imsize2[1])]
          args = args + ["--" + dev_type]
          if scene_path.find(PATH_TO_HYDRA3_SCENS) != -1:
            args = args + ["-scn_dir", PATH_TO_HYDRA3_SCENS]
          if req.spectral:
            args = args + ["--spectral"]
          #print(args)
          req.run(test_name, args, image_ref, outp, inregrator)
          #print("finished")


class REQ_HP(REQ):
  def __init__(self, name, scn_path, ref_path, imsize):
    self.name     = name
    self.scn_path = scn_path
    self.ref_path = ref_path
    self.imsize   = imsize
    self.integs   = ["mispt"]
    self.naivem   = 1
    self.times  = []
    self.names    = ["test_102", "sponza", "cry_sponza"]

  def need_render_time(req):
    return True

  def test(req, gpu_id=0):
    for (scnp, imgp, id, imsize2) in zip(req.scn_path, req.ref_path, range(0,len(req.scn_path)), req.imsize):
      image_ref   = cv2.imread(imgp, cv2.IMREAD_ANYCOLOR | cv2.IMREAD_ANYDEPTH)
      scene_path  = scnp
      folder_path = os.path.dirname(imgp)
      devices   = ["gpu"] if not TEST_CPU else ["gpu", "cpu"]
      test_name = req.name
      for dev_type in devices:
        Log().info("  rendering scene: '{0}', dev_type='{1}', scene = '{2}'".format(test_name, dev_type, scene_path))
        for integrator in req.integs:
          outp = folder_path + "/y" + str(id) + "_" + dev_type + integrator + ".bmp"
          print(outp)
          args = ["./bin-release/hydra", "-in", scene_path, "-out", outp, "-integrator", integrator, "-spp", "1024"]
          args = args + ["-gpu_id", str(gpu_id)]  # for single launch samples
          args = args + ["-width", str(imsize2[0]), "-height", str(imsize2[1])]
          args = args + ["--" + dev_type]
          if scene_path.find(PATH_TO_HYDRA3_SCENS) != -1:
            args = args + ["-scn_dir", PATH_TO_HYDRA3_SCENS]
          # print(args)
          req.run(test_name, args, image_ref, outp, integrator)
    print(req.names)
    print(req.times)


reqs = []

'''
reqs.append( REQ_HP("perf_test", [PATH_TO_HYDRA2_TESTS + "/tests_f/test_102/statex_00001.xml",  
                                  "/home/frol/PROG/msu-graphics-group/scenes/03_classic_scenes/01_sponza/statex_00001.xml",
                                  "/home/frol/PROG/msu-graphics-group/scenes/03_classic_scenes/02_cry_sponza/statex_00001.xml"],
                                 [PATH_TO_HYDRA2_TESTS + "/tests_images/test_102/z0_gpumispt.bmp", 
                                  PATH_TO_HYDRA2_TESTS + "/tests_images/test_102/z1_gpumispt.bmp",
                                  PATH_TO_HYDRA2_TESTS + "/tests_images/test_102/z2_gpumispt.bmp"],
                                  [(1024,1024), (1024,1024), (1024,1024)]))

'''

reqs.append( REQ_H2("mat_mirror",           ["test_102"], integrators = ["naivept","mispt"]) )
reqs.append( REQ_H2("mat_lambert_texture",  ["test_103"]) )
reqs.append( REQ_H2("mat_texture_matrices", ["test_110"]) )
reqs.append( REQ_H2("mat_emission_texture", ["test_124"], integrators = ["naivept","mispt"]) )
reqs.append( REQ_H2("mat_normal_bump",      ["test_127"], naivemul = 4, imsize = (1024,768)) )

reqs.append( REQ_H2("lgt_sphere",          ["test_201"]) )
reqs.append( REQ_H2("lgt_point_omni",      ["test_213"], integrators = ["mispt"]) )
reqs.append( REQ_H2("lgt_area4_transform", ["test_215"]) )
reqs.append( REQ_H2("lgt_area_rotate",     ["test_223"]) )
reqs.append( REQ_H2("lgt_area_rotate",     ["test_224"]) )

reqs.append( REQ_H2("lgt_point_ies",       ["test_228"], integrators = ["mispt"]) )
reqs.append( REQ_H2("lgt_area_ies",        ["test_206", "test_207", "test_208", "test_216", "test_232"], integrators = ["mispt"]) )
reqs.append( REQ_H2("lgt_area_disk",       ["test_246"], naivemul = 4) )

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
                                    PATH_TO_HYDRA3_SCENS + "/Tests/Lambert/0001/Lambert_cornell_hydra2.xml"],
                                   [PATH_TO_HYDRA2_TESTS + "/tests_images/test_101/w_ref.png",
                                    PATH_TO_HYDRA3_SCENS + "/Tests/Lambert/0001/Images/Lambert_cornell_mitsuba.png"],
                                    imsize = [(512,512), (1024,1024)], naivemul = 4))

reqs.append( REQ_HX("mat_emission", [PATH_TO_HYDRA2_TESTS + "/tests_f/test_123/statex_00001.xml",
                                     PATH_TO_HYDRA2_TESTS + "/tests_f/test_123/statex_00002.xml",
                                     PATH_TO_HYDRA2_TESTS + "/tests_f/test_123/statex_00003.xml"],
                                    [PATH_TO_HYDRA2_TESTS + "/tests_images/test_123/w_ref.png",
                                     PATH_TO_HYDRA2_TESTS + "/tests_images/test_123/w_ref2.png",
                                     PATH_TO_HYDRA2_TESTS + "/tests_images/test_123/w_ref3.png"],
                                     imsize = [(512,512), (512,512), (512,512)],
                                     naivemul = 16, integrators = ["naivept","mispt"]))


#reqs.append( REQ_HX("mat_smooth_plastic", [PATH_TO_HYDRA3_SCENS + "/Tests/Plastic_smooth/0001/PlasticSmooth_sphere_hydra2.xml",
#                                           PATH_TO_HYDRA3_SCENS + "/Tests/Plastic_smooth/0002/PlasticSmooth_cornell_hydra2.xml"],
#                                          [PATH_TO_HYDRA3_SCENS + "/Tests/Plastic_smooth/0001/Images/PlasticSmooth_sphere_mitsuba.png",
#                                           PATH_TO_HYDRA3_SCENS + "/Tests/Plastic_smooth/0002/Images/PlasticSmooth_cornell_mitsuba.png"],
#                                           imsize = [(1024, 1024), (1024, 1024)], naivemul = 4))

reqs.append( REQ_HX("mat_conductor",
                    [
                      PATH_TO_HYDRA3_SCENS + "/Tests/Conductor/0001/Smooth-eta1.5-sphere-hydra3.xml",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Conductor/0002/Rough-uv01-sphere-hydra3.xml",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Conductor/0003/Rough-u025-v001-plane-hydra3.xml",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Conductor/0004/Rough-u001-v025-plane-hydra3.xml",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Conductor/0006/Rough-texture-eta1.5-sphere-hydra3.xml"
                    ],
                    [
                      PATH_TO_HYDRA3_SCENS + "/Tests/Conductor/0001/Images/Smooth-eta1.5-sphere-mitsuba.png",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Conductor/0002/Images/Rough-uv01-sphere-mitsuba.png",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Conductor/0003/Images/Rough-u025-v001-plane-mitsuba.png",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Conductor/0004/Images/Rough-u001-v025-plane-mitsuba.png",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Conductor/0006/Images/Rough-texture-eta1.5-sphere-mitsuba.png"
                    ],
                    imsize = [(1024, 1024), (1024, 1024), (1024, 1024), (1024, 1024), (1024, 1024)],
                    naivemul = 16))

reqs.append( REQ_HX("spectral",
                    [
                      PATH_TO_HYDRA3_SCENS + "/Tests/Spectral/0001/Spectral-ior-sphere-hydra3.xml",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Spectral/0002/Spectral-ior-model-hydra3.xml",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Spectral/0003/Spectral-diffuse-sphere-hydra3.xml",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Spectral/0004/spectral_cornell_hydra3.xml",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Spectral/0005/spectral_cornell_hydra3.xml",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Spectral/0007/PlasticRough-025_sphere_hydra3.xml"
                    ],
                    [
                      PATH_TO_HYDRA3_SCENS + "/Tests/Spectral/0001/Images/Spectral-ior-sphere-mitsuba.png",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Spectral/0002/Images/Spectral-ior-model-mitsuba.png",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Spectral/0003/Images/Spectral-diffuse-sphere-mitsuba.png",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Spectral/0004/Images/spectral_cornell_mitsuba.png",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Spectral/0005/Images/spectral_cornell_mitsuba.png",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Spectral/0007/Images/PlasticRough-025_sphere_mitsuba.png"
                    ],
                    imsize = [(1024, 1024), (1024, 1024), (1024, 1024), (1024, 1024), (1024, 1024), (1024, 1024)],
                    naivemul = 16, integrators = ["mispt"], is_spectral = True))

reqs.append( REQ_HX("blend",
                    [
                      PATH_TO_HYDRA3_SCENS + "/Tests/Blend/0001/spectral-blend-sphere-hydra3.xml",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Blend/0002/spectral-blend-sphere-hydra3.xml"
                    ],
                    [
                      PATH_TO_HYDRA3_SCENS + "/Tests/Blend/0001/Images/spectral-blend-sphere-mitsuba.png",
                      PATH_TO_HYDRA3_SCENS + "/Tests/Blend/0002/Images/spectral-blend-sphere-mitsuba.png"
                    ],
                    imsize = [(1024, 1024), (1024, 1024)],
                    naivemul = 16, integrators = ["mispt"], is_spectral = True))

reqs.append( REQ_HX("mat_smooth_glass", [PATH_TO_HYDRA3_SCENS + "/Tests/Glass/0001/Glass-sphere_rough-0_cornell_hydra3.xml",
                                         PATH_TO_HYDRA3_SCENS + "/Tests/Glass/0002/Glass_rough-0_cornell_hydra3.xml",
                                         PATH_TO_HYDRA3_SCENS + "/Tests/Glass/0003/quartz-prism-rough-0-cornell-hydra3.xml"],
                                        [PATH_TO_HYDRA3_SCENS + "/Tests/Glass/0001/Images/Glass-sphere_rough-0_cornell_mitsuba.png",
                                         PATH_TO_HYDRA3_SCENS + "/Tests/Glass/0002/Images/Glass_rough-0_cornell_mitsuba.png",
                                         PATH_TO_HYDRA3_SCENS + "/Tests/Glass/0003/Images/quartz-prism-rough-0-cornell-mitsuba.png"],
                                         imsize = [(1024, 1024), (1024, 1024), (1024, 1024)], naivemul = 8))

Log().set_workdir(".")
Log().info("PATH_TO_TESTS = {}".format(PATH_TO_HYDRA2_TESTS))

os.chdir('..') # use HydraCore3 root dir as current
Log().set_workdir("testing")

for req in reqs:
  Log().info("REQ: {}".format(req.name))
  req.test()

###
