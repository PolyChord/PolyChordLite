cmake_minimum_required(VERSION 3.4)

project(PolyChordLite Fortran)

# languages
enable_language(CXX)
enable_language(Fortran)

# options
option(MPI "Build with MPI" ON)
option(python "Build PyPolyChord package" ON)
option(python_user_install "Install Python package for user (--user)" ON)

# enforce matching of the CXX and Fortran compilers
if (NOT "${CMAKE_Fortran_COMPILER_ID}" MATCHES "${CMAKE_CXX_COMPILER_ID}"
        AND NOT ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU" AND
        "${CMAKE_CXX_COMPILER_ID}" MATCHES "AppleClang"))
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

    if (CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER 10.0)
        # GCC 10 is complaining about MPI interfaces, turn this to warnings
        SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
    endif ()

elseif ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Cray")
    message("Using Cray compilers")

    set(CMAKE_Fortran_FLAGS "-fpp -fpic -qopenmp -dynamic")
    set(CMAKE_CXX_FLAGS "-fpic -qopenmp -dynamic")

endif ()

# find MPI & link to everything from now on
if (MPI)
    find_package(MPI REQUIRED)
    message(STATUS "MPI_Fortran_COMPILER: ${MPI_Fortran_COMPILER}")
    link_libraries(MPI::MPI_Fortran)
    link_libraries(MPI::MPI_CXX)

    # definitions
    add_compile_definitions(USE_MPI)
endif (MPI)

# this will be module building directory, allows linking
set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/generated/mods)

# on installation we get this prefix as well
set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})

# target directories -> stuff installed in these directories
set(target_lib_dir ${PROJECT_SOURCE_DIR}/lib)
set(target_bin_dir ${PROJECT_SOURCE_DIR}/bin)

# directories where we want to make things
add_subdirectory(src)
add_subdirectory(likelihoods)

if (python)
    #################################################################
    # The Compiled .so file
    #################################################################

    # PyPolyChord stuff
    message(STATUS "Enabling Python package build")
    find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

    # find out the suffix for the current python version
    EXECUTE_PROCESS(
            COMMAND ${Python3_EXECUTABLE} "-c" "import sysconfig;print(sysconfig.get_config_var('EXT_SUFFIX'))"
            OUTPUT_VARIABLE python_so_suffix
    )
    STRING(REPLACE "\n" "" python_so_suffix ${python_so_suffix})
    message(STATUS "Python shared object suffix is: ${python_so_suffix}")

    # the extension module
    add_library(_pypolychord SHARED pypolychord/_pypolychord.cpp $<TARGET_OBJECTS:objlib_chord>)
    set_target_properties(_pypolychord
            PROPERTIES
            OUTPUT_NAME _pypolychord
            PREFIX ""
            SUFFIX "${python_so_suffix}"
            )
    # dependencies and includes
    target_link_libraries(_pypolychord Python3::Python Python3::NumPy)
    target_include_directories(_pypolychord PUBLIC src/polychord) # PolyChord headers for python extensions

    # remember the place of this
    set(target_python_so ${CMAKE_CURRENT_BINARY_DIR}/_pypolychord${python_so_suffix})

    # place all outputs in one place
    install(TARGETS _pypolychord
            LIBRARY DESTINATION ${target_lib_dir})

    #################################################################
    # The simple python module
    #################################################################

    # hack ofr version
    set(PACKAGE_VERSION "1.18.2")

    set(SETUP_PY_IN "${CMAKE_CURRENT_SOURCE_DIR}/setup_cmake_template.py")
    set(SETUP_PY "${CMAKE_CURRENT_BINARY_DIR}/setup.py")
    set(DEPS "${CMAKE_CURRENT_SOURCE_DIR}/pypolychord/__init__.py")
    set(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/build/timestamp") # dummy target for dependence

    # fill in the variables of the setup.py file
    configure_file(${SETUP_PY_IN} ${SETUP_PY})

    # custom target
    add_custom_command(OUTPUT ${OUTPUT}
            COMMAND ${Python3_EXECUTABLE} ${SETUP_PY} build
            COMMAND ${CMAKE_COMMAND} -E touch ${OUTPUT}
            DEPENDS ${DEPS} _pypolychord)
    add_custom_target(target ALL DEPENDS ${OUTPUT})

    # installation with setup.py
    if (python_user_install)
        install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} ${SETUP_PY} install)")
    else (python_user_install)
        install(CODE "execute_process(COMMAND ${Python3_EXECUTABLE} ${SETUP_PY} install --user)")
    endif (python_user_install)

endif (python)
