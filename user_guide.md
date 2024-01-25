# Working with HydraCore3

## Command line parameters
* `-in *path*` - path to input scene
`-in scenes/test_035/statex_00001.xml`
* `-out *path*` - path to rendered image
`-out myrender.exr`
* `-scn_dir *path*` - path to scene root directory (needed when scene description uses relative paths to assets)
* `-integrator *type*` - rendering algorithm, can be one of: *naivept* (naive path tracing), *shadowpt* (shadow path tracing), *mispt* (path tracing with multiple importance sampling), *prt* (primary rays only), *raytracing* (whitted)


