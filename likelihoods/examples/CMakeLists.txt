project(examples Fortran)

# parse each example separately
file(GLOB ${PROJECT_NAME}_sources "*.f90" "*.F90")
set(${PROJECT_NAME}_library_list "") # hold all the targets we make
foreach (filename ${${PROJECT_NAME}_sources})
    # strip the file name
    get_filename_component(archive_name ${filename} NAME_WE)

    if (${archive_name} STREQUAL "polychord_examples")
        # make sure this is not turned into an object & executable alone
        continue()
    endif ()

    # compile once into object -> use in lib and exec as well
    add_library(objlib_${archive_name} OBJECT ${filename})
    set_target_properties(objlib_${archive_name} PROPERTIES POSITION_INDEPENDENT_CODE 1 LINKER_LANGUAGE Fortran)
    target_link_libraries(objlib_${archive_name} chord)

    # static library
    add_library(${archive_name} STATIC $<TARGET_OBJECTS:objlib_${archive_name}>)

    # executable + linking, uses the program and one of the likelihood implementations
    add_executable(${archive_name}_exec $<TARGET_OBJECTS:objlib_${archive_name}> polychord_examples.f90)
    target_link_libraries(${archive_name}_exec chord) # link libchoord.a to it
    set_target_properties(${archive_name}_exec PROPERTIES OUTPUT_NAME ${archive_name} CLEAN_DIRECT_OUTPUT 1 LINKER_LANGUAGE Fortran)

    # add to target list
    list(APPEND ${PROJECT_NAME}_library_list ${archive_name})
    list(APPEND ${PROJECT_NAME}_library_list ${archive_name}_exec)

    # tell the user
    message("Adding example: ${archive_name} (from ${filename})")
endforeach ()


# standard: where to put the targets
install(TARGETS ${${PROJECT_NAME}_library_list}
        RUNTIME DESTINATION ${target_bin_dir}
        LIBRARY DESTINATION ${target_lib_dir}
        ARCHIVE DESTINATION ${target_lib_dir})
