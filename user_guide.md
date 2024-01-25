# Working with HydraCore3

Render a simple test scene from command line after building hydra executable:

`./bin-release/hydra -in scenes/test_035/statex_00001.xml -out test_035.exr -spp 256`

## Command line parameters
* `-in path` - path to input scene, example:

`-in scenes/test_035/statex_00001.xml`

* `-out path` - path to rendered image, example:

`-out myrender.exr`

`-out myrender.bmp`
* `-scn_dir *path*` - path to scene root directory (needed when scene description uses relative paths to assets), example:

`-scn_dir my_scenes/scenes_hydra`

* `-integrator type` - rendering algorithm, can be one of: *naivept* (naive path tracing), *shadowpt* (shadow path tracing), *mispt* (path tracing with multiple importance sampling), *prt* (primary rays only), *raytracing* (whitted), example:

`-integrator mispt`

* `-spp N_SAMPLES` - number of samples per pixel, example:

`-spp 1024`

* `--spectral` - spectral rendering mode
* `--gpu` - render on GPU
* `-gpu_id ID` - GPU device ID, example:

`-gpu_id 1`

* `--cpu` - render on CPU
* `-width WIDTH -height HEIGHT` - override rendering resolution, example:

`-width 1920 -height 1080`

* `-channels CHANNEL_NUM` - number of framebuffer channels (for example, spectral bands), 1 (monochrome), 4 (RGB) or >4 (spectral) example:

`-channels 16`

* `-gamma VAL` - gamma for output LDR image, example:

`-gamma 2.2`
`-gamma srgb`

* `-look_at MATRIX_STR` - override view transform with a look at matrix in a row major order, example:

`-look_at 0.707107 -0.408248 0.57735 0 0 0.816497 0.57735 0 -0.707107 -0.408248 0.57735 0 0 0 -1.73205 1`

