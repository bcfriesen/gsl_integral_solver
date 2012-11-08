CMAKE_MINIMUM_REQUIRED (VERSION 2.8)
PROJECT (gsl_integral_solver)
ENABLE_LANGUAGE (C)

#----------------------------BEGIN ADD EXTRA MODULES----------------------------
# Appends the cmake/modules path inside the MAKE_MODULE_PATH variable, which
# stores the directories of additional CMake modules
# (e.g. MacroOutOfSourceBuild.cmake):
SET (CMAKE_MODULE_PATH "${${PROJECT_NAME}_SOURCE_DIR}/cmake/modules" ${CMAKE_MODULE_PATH})

# The macro below forces the build directory to be different from source
# directory
INCLUDE (MacroOutOfSourceBuild)
macro_ensure_out_of_source_build("${PROJECT_NAME} requires an out of source build.")
#-----------------------------END ADD EXTRA MODULES-----------------------------

SET (CMAKE_VERBOSE_MAKEFILE "True")
SET (CMAKE_BUILD_TYPE "Debug")

ADD_SUBDIRECTORY (src)
