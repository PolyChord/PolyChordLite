project(examples)

# parse each example separately
file(GLOB ${PROJECT_NAME}_sources "*.f90" "*.F90")
set(${PROJECT_NAME}_library_list "") # hold all the targets we make
foreach(filename ${${PROJECT_NAME}_sources})
    # strip the file name
    get_filename_component(archive_name ${filename} NAME_WE)
    add_library(${archive_name} ${filename})

    # add to target list
    list(APPEND ${PROJECT_NAME}_library_list ${archive_name})

    # tell the user
    message("Adding example: ${archive_name} (from ${filename})")
endforeach()


# standard: where to put the targets
install(TARGETS ${${PROJECT_NAME}_library_list}
        RUNTIME DESTINATION ${target_bin_dir}
        LIBRARY DESTINATION ${target_lib_dir}
        ARCHIVE DESTINATION ${target_lib_dir})
