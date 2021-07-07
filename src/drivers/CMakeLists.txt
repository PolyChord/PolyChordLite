project(polychord_examples)

# build the object file -> use elsewhere, eg. in examples
add_library(${PROJECT_NAME} OBJECT polychord_examples.f90)
target_include_directories(${PROJECT_NAME} INTERFACE ${CMAKE_Fortran_MODULE_DIRECTORY})

# standard: where to put the targets
install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION ${target_bin_dir}
        LIBRARY DESTINATION ${target_lib_dir}
        ARCHIVE DESTINATION ${target_lib_dir})
