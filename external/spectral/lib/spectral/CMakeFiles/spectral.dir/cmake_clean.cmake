file(REMOVE_RECURSE
  "../../../../bin-debug/libspectral.pdb"
  "../../../../bin-debug/libspectral.so"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/spectral.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
