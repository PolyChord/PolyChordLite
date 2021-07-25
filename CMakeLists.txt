cmake_minimum_required(VERSION 3.4)

project(PolyChordLite Fortran)

# languages
enable_language(CXX)
enable_language(Fortran)

# enforce matching of the CXX and Fortran compilers
if (NOT "${CMAKE_Fortran_COMPILER_ID}" MATCHES "${CMAKE_CXX_COMPILER_ID}")
    message(FATAL_ERROR "You need to use the same vendor for your C++ and Fortran compiler")
endif ()

# flags for all three compiler types
# todo: inspect the flags and implement in a nicer way
if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Intel")
    message("Using INTEL compilers")

    set(CMAKE_Fortran_FLAGS "-fpp -fpic -assume noold_maxminloc -heap-arrays -ipo -O3 -no-prec-div -xHost -w -vec-report0 -qopt-report0")
    set(CMAKE_CXX_FLAGS "-std=c++11 -fpic -ipo -O3 -no-prec-div -xHost -w -vec-report0 -qopt-report0")

    add_compile_options("-nofor-main")
    link_libraries("-nofor-main")

elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    message("Using GNU compilers")

    set(CMAKE_Fortran_FLAGS "-ffree-line-length-none -cpp -fPIC -fno-stack-arrays -Ofast")
    set(CMAKE_CXX_FLAGS "-std=c++11 -fPIC -Ofast")

    # for GitHub actions MacOS -- debug only
    add_compile_options("-lstdc++")
    link_libraries("-lstdc++")

elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Cray")
    message("Using Cray compilers")

    set(CMAKE_Fortran_FLAGS "-fpp -fpic -qopenmp -dynamic")
    set(CMAKE_CXX_FLAGS "-fpic -qopenmp -dynamic")

endif ()

# find MPI & link to everything from now on
find_package(MPI REQUIRED)
message(STATUS "MPI_Fortran_COMPILER: ${MPI_Fortran_COMPILER}")
link_libraries(MPI::MPI_Fortran)

# this will be module building directory, allows linking
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/generated/mods)

# on installation we get this prefix as well
set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})

# target directories -> stuff installed in these directories
set(target_lib_dir ${PROJECT_BINARY_DIR}/lib)
set(target_bin_dir ${PROJECT_BINARY_DIR}/bin)
set(target_python_dir ${PROJECT_BINARY_DIR}/python)

# directories where we want to make things
add_subdirectory(src)
add_subdirectory(likelihoods)

# PyPolyChord stuff
message(STATUS "Enabling Python package build")
find_package (Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

# find out the suffix for the current python version
EXECUTE_PROCESS(
    COMMAND ${Python3_EXECUTABLE} "-c" "from importlib.machinery import EXTENSION_SUFFIXES;print(EXTENSION_SUFFIXES[0])"
    OUTPUT_VARIABLE python_so_suffix
)
STRING(REPLACE "\n" "" python_so_suffix ${python_so_suffix})
message(STATUS "Python shared object suffix is: ${python_so_suffix}")

# the extension module
add_library(_pypolychord SHARED pypolychord/_pypolychord.cpp)
set_target_properties(_pypolychord
    PROPERTIES
        OUTPUT_NAME _pypolychord
        PREFIX ""
        SUFFIX "${python_so_suffix}"
    )
# dependencies and includes
target_link_libraries(_pypolychord chord_shared Python3::Python Python3::NumPy)
target_include_directories(_pypolychord PUBLIC src/polychord) # PolyChord headers for python extensions

# place all outputs in one place
install(TARGETS _pypolychord
        RUNTIME DESTINATION ${target_python_dir}
        LIBRARY DESTINATION ${target_python_dir}
        ARCHIVE DESTINATION ${target_python_dir})