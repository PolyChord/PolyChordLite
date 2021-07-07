project(examples)

# parse each example separately
file(GLOB ${PROJECT_NAME}_sources "*.f90" "*.F90")
set(${PROJECT_NAME}_library_list "") # hold all the targets we make
foreach (filename ${${PROJECT_NAME}_sources})
    # strip the file name
    get_filename_component(archive_name ${filename} NAME_WE)

    # compile once into object -> use in lib and exec as well
    add_library(objlib_${archive_name} OBJECT ${filename})
    set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

    # static library
    add_library(${archive_name} STATIC $<TARGET_OBJECTS:objlib_${archive_name}>)

    # executable + linking
    add_executable(${archive_name}_exec $<TARGET_OBJECTS:objlib_${archive_name}> $<TARGET_OBJECTS:polychord_examples>)
    target_link_libraries(${archive_name}_exec chord) # link libchoord.a to it
    set_target_properties(${archive_name}_exec PROPERTIES OUTPUT_NAME ${archive_name} CLEAN_DIRECT_OUTPUT 1)

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
