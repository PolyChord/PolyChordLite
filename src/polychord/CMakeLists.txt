project(chord Fortran)

## original: glob -> sources -> static library
file(GLOB ${PROJECT_NAME}_sources "*.f90" "*.F90" "*.cpp")
message(STATUS "Sources for project ${PROJECT_NAME}: ${${PROJECT_NAME}_sources}")
add_library(${PROJECT_NAME} STATIC ${${PROJECT_NAME}_sources})

# this seems to be standard stuff
install(TARGETS ${PROJECT_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib/static)
