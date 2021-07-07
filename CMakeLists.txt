cmake_minimum_required(VERSION 3.4)

project(PolyChordLite Fortran)

# languages
enable_language(Fortran)
enable_language(CXX)

# flags for gfortran todo: compiler dependence
set (CMAKE_Fortran_FLAGS "-ffree-line-length-none -cpp -fPIC -fno-stack-arrays")

# this will be module building directory, allows linking
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/generated/mods)

# on installation we get this prefix as well
set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})

# target directories -> stuff installed in these directories
set(target_lib_dir ${PROJECT_BINARY_DIR}/lib)
set(target_bin_dir ${PROJECT_BINARY_DIR}/bin)

# directories where we want to make things
add_subdirectory(src)
add_subdirectory(likelihoods)
