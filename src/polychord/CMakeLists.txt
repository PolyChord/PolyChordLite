project(chord Fortran)

# file glob -> sources -> object -> static & shared lib as well
file(GLOB ${PROJECT_NAME}_sources "*.f90" "*.F90" "*.cpp")
message(STATUS "Sources for project ${PROJECT_NAME}: ${${PROJECT_NAME}_sources}")

# this is the "object library" target: compiles the sources only once
add_library(objlib OBJECT ${${PROJECT_NAME}_sources})

# shared libraries need PIC
set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

# libcoord.a
add_library(${PROJECT_NAME} STATIC $<TARGET_OBJECTS:objlib>)

# libcoord.so as well -> needs a different target name, but the same lib name
add_library(${PROJECT_NAME}_shared SHARED $<TARGET_OBJECTS:objlib>)
set_target_properties(${PROJECT_NAME}_shared PROPERTIES OUTPUT_NAME ${PROJECT_NAME} CLEAN_DIRECT_OUTPUT 1)

# standard: where to put the targets
install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION ${target_bin_dir}
        LIBRARY DESTINATION ${target_lib_dir}
        ARCHIVE DESTINATION ${target_lib_dir})
