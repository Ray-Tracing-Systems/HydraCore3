import os, sys, subprocess, re, csv

#def list_scenes(directory):
#  return [(entry, directory + "/" + entry + "/statex_00001.xml") for entry in os.listdir(directory) if os.path.isdir(os.path.join(directory, entry))]

def list_scenes(directory):
    entries = [(entry, directory + "/" + entry + "/statex_00001.xml") 
              for entry in os.listdir(directory) 
              if os.path.isdir(os.path.join(directory, entry))]
    return sorted(entries, key=lambda x: x[0])

def run(test_name, args, time_list):
  try:
    res = subprocess.run(args, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output  = res.stdout.decode('utf-8')
    #print(output)
    pattern = r'PathTraceBlock\(exec\) = (\d+\.\d+) ms'
    match   = re.search(pattern, output)
    if match:
      execution_time_ms = round(float(match.group(1)))
      time_list.append(execution_time_ms)
  except Exception as e:
    print("Failed to run scene {0} : {1}".format(test_name, e),)
    return

if __name__ == '__main__':

  HYDRA3_PATH = ""
  SPP = 1024

  os.chdir('..') # use HydraCore3 root dir as current
  if sys.platform == 'win32':
    HYDRA3_PATH = "./bin-release/Release/hydra.exe"
  else:  # if sys.platform == 'linux':
    HYDRA3_PATH = "./bin-release/hydra"
  
  for SQRRES in [False, True]:
  
    time_list = []
    scenes = list_scenes("../HydraScenes")
    for scn in scenes:
        args = [HYDRA3_PATH, "-in", scn[1], "-out",  "y_" + scn[0] + ".png", "-spp", str(SPP), "--gpu"]
        if(SQRRES):
          args = args + ["-width", "1024", "-height", "1024"]
        print(args)  
        run(scn[0], args, time_list)
    
    with open('testing/data.csv', 'a', newline='', encoding='utf-8') as file:
      writer = csv.writer(file)
      
      # Если файл пустой, добавляем заголовки
      if file.tell() == 0:
          writer.writerow(["Scene", "Time"])  # Заголовки столбцов
      
      # Записываем данные построчно
      for name, value in zip(scenes, time_list):
          suffix = "(1024)" if SQRRES else ""
          writer.writerow([name[0] + suffix, value])
    
    print(time_list)
