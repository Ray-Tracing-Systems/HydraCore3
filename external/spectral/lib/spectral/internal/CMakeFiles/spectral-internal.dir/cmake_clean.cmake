file(REMOVE_RECURSE
  "../../../../../bin-debug/libspectral-internal.a"
  "../../../../../bin-debug/libspectral-internal.pdb"
)

# Per-language clean rules from dependency scanning.
foreach(lang CXX)
  include(CMakeFiles/spectral-internal.dir/cmake_clean_${lang}.cmake OPTIONAL)
endforeach()
