project(polychord Fortran)

# manual add one
add_library(${PROJECT_NAME} STATIC utils.F90 abort.F90)
message("manually added STATIC utils.F90 into `polychord` module in ${CMAKE_CURRENT_SOURCE_DIR}")

# this seems to be standard stuff
install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static)
