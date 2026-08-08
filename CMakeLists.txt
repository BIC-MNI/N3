# CMakeFiles.txt for N3
#
# Vladimir S. FONOV - vladimir.fonov@gmail.com

CMAKE_MINIMUM_REQUIRED(VERSION 3.10)

SET(N3_PACKAGE_VERSION_MAJOR 1)
SET(N3_PACKAGE_VERSION_MINOR 12)
SET(N3_PACKAGE_VERSION_PATCH 00)

IF(NOT MINC_TOOLKIT_BUILD)
  PROJECT(N3)
  #SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake-modules")
  ENABLE_TESTING()
  INCLUDE(CTest)

  # Packaging defines
  SET(CPACK_GENERATOR TGZ)
  SET(CPACK_PACKAGE_VERSION_MAJOR ${N3_PACKAGE_VERSION_MAJOR})
  SET(CPACK_PACKAGE_VERSION_MINOR ${N3_PACKAGE_VERSION_MINOR})
  SET(CPACK_PACKAGE_VERSION_PATCH ${N3_PACKAGE_VERSION_PATCH})
  INCLUDE(CPack)
ENDIF(NOT MINC_TOOLKIT_BUILD)

SET(PACKAGE "N3")

SET(PACKAGE_BUGREPORT "a.janke@gmail.com")

SET(PACKAGE_NAME "N3")
SET(PACKAGE_VERSION "${N3_PACKAGE_VERSION_MAJOR}.${N3_PACKAGE_VERSION_MINOR}.${N3_PACKAGE_VERSION_PATCH}")
SET(PACKAGE_STRING  "${PACKAGE_NAME} ${PACKAGE_VERSION}")

SET(CMAKE_CXX_STANDARD 17)
SET(CMAKE_C_STANDARD 99)

#setup variables only needed inside scripts
SET(version      ${PACKAGE_VERSION})
SET(long_version ${PACKAGE_STRING})
#todo setup: date user machine host 

IF(NOT MINC_TOOLKIT_BUILD)
  FIND_PACKAGE( LIBMINC REQUIRED )
  FIND_PACKAGE( EBTKS REQUIRED )
ELSE(NOT MINC_TOOLKIT_BUILD)
  # FIND_PACKAGE(EBTKS) is skipped under a superbuild, so EBTKSConfig.cmake is never read
  # and these never arrive that way. Empty default matches EBTKS's own bundled default: a
  # superbuild that built its own EBTKS with a non-bundled EBTKS_BLAS_BACKEND passes the
  # same two values here directly, alongside the EBTKS_LIBRARIES/EBTKS_USE_FILE it already
  # passes (/app/PLAN.md §10, "The superbuild hooks already exist").
  SET(EBTKS_LAPACK_LIBRARIES "" CACHE STRING
      "External LAPACK/BLAS legacy/N3 links, matching the EBTKS_BLAS_BACKEND the outer \
superbuild built EBTKS with. Empty if EBTKS was built bundled.")
  SET(EBTKS_LAPACK_LIBRARY_DIRS "" CACHE STRING
      "Extra -L directories for EBTKS_LAPACK_LIBRARIES, when its entries are bare names.")
ENDIF(NOT MINC_TOOLKIT_BUILD)

# EBTKS's "lapacke" backend provides dsysv_ under the name EBTKS_dsysv, not dsysv_ itself
# (see EBTKS/shim/dsysv_lapacke.c -- a real dsysv_ here would collide with the one the
# linked LAPACKE library uses internally). TBSpline.cc picks which symbol to call via this
# define; every other backend is unaffected and keeps calling dsysv_ directly as before.
IF(EBTKS_BLAS_BACKEND STREQUAL "lapacke")
  SET_SOURCE_FILES_PROPERTIES(src/Splines/TBSpline.cc PROPERTIES
    COMPILE_DEFINITIONS EBTKS_DSYSV_LAPACKE_SHIM)
ENDIF()

# LIBMINCConfig.cmake records the absolute paths of the machine libminc was built
# on.  Where such a path no longer exists, look the same file up by name in
# LIBMINC_LIBRARY_DIRS instead.  A no-op on an install whose paths are intact.
FOREACH(_var LIBMINC_LIBRARIES VOLUME_IO_LIBRARIES)
  SET(_repaired)
  FOREACH(_lib ${${_var}})
    IF(IS_ABSOLUTE "${_lib}" AND NOT EXISTS "${_lib}")
      GET_FILENAME_COMPONENT(_name "${_lib}" NAME)
      UNSET(_found CACHE)
      FIND_FILE(_found ${_name} PATHS ${LIBMINC_LIBRARY_DIRS} NO_DEFAULT_PATH)
      IF(_found)
        MESSAGE(STATUS "LIBMINC: ${_lib} -> ${_found}")
        LIST(APPEND _repaired ${_found})
      ELSE(_found)
        MESSAGE(FATAL_ERROR "LIBMINC names ${_lib}, which does not exist, and "
                            "${_name} is not in ${LIBMINC_LIBRARY_DIRS}")
      ENDIF(_found)
    ELSE(IS_ABSOLUTE "${_lib}" AND NOT EXISTS "${_lib}")
      LIST(APPEND _repaired ${_lib})
    ENDIF(IS_ABSOLUTE "${_lib}" AND NOT EXISTS "${_lib}")
  ENDFOREACH(_lib)
  SET(${_var} ${_repaired})
ENDFOREACH(_var)
UNSET(_found CACHE)

FIND_PACKAGE( Perl REQUIRED )
INCLUDE( ${LIBMINC_USE_FILE} )
INCLUDE( ${EBTKS_USE_FILE} )

ADD_DEFINITIONS(-DHAVE_CONFIG_H)
IF(NOT APPLE)
  ADD_DEFINITIONS(-D_XOPEN_SOURCE=500)
ENDIF(NOT APPLE)  
ADD_DEFINITIONS(-DUSE_COMPMAT)
ADD_DEFINITIONS(-DUSE_DBLMAT)
ADD_DEFINITIONS(-DUSE_FCOMPMAT)

# aliases
SET(VERSION "${PACKAGE_VERSION}")


# check for prereqs
INCLUDE(CheckFunctionExists)
INCLUDE(CheckSymbolExists)
INCLUDE(CheckIncludeFiles)


if(NOT MSVC)
  set(CMAKE_REQUIRED_LIBRARIES m)
endif()

set(CMAKE_REQUIRED_LIBRARIES "m")
set(CMAKE_REQUIRED_DEFINITIONS -D_GNU_SOURCE=1 )
CHECK_SYMBOL_EXISTS(finite    "math.h" HAVE_FINITE)
CHECK_SYMBOL_EXISTS(isfinite  "math.h" HAVE_ISFINITE)

if(NOT MSVC)
    unset(CMAKE_REQUIRED_LIBRARIES)
endif()

unset(CMAKE_REQUIRED_DEFINITIONS)


CHECK_INCLUDE_FILES(float.h     HAVE_FLOAT_H)
CHECK_INCLUDE_FILES(sys/dir.h   HAVE_SYS_DIR_H)
CHECK_INCLUDE_FILES(sys/ndir.h  HAVE_SYS_NDIR_H)
CHECK_INCLUDE_FILES(sys/stat.h  HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILES(sys/wait.h  HAVE_SYS_WAIT_H)
CHECK_INCLUDE_FILES(values.h    HAVE_VALUES_H)
CHECK_INCLUDE_FILES(unistd.h    HAVE_UNISTD_H)
CHECK_INCLUDE_FILES(dirent.h    HAVE_DIRENT_H)
CHECK_INCLUDE_FILES(memory.h    HAVE_MEMORY_H)
CHECK_INCLUDE_FILES(malloc.h    HAVE_MALLOC_H)
CHECK_INCLUDE_FILES(stdlib.h    HAVE_STDLIB_H)
CHECK_INCLUDE_FILES(fcntl.h     HAVE_FCNTL_H)
CHECK_INCLUDE_FILES(vfork.h     HAVE_VFORK_H)
CHECK_INCLUDE_FILES(inttypes.h  HAVE_INTTYPES_H)
CHECK_INCLUDE_FILES(string.h    HAVE_STRING_H)
CHECK_INCLUDE_FILES(strings.h   HAVE_STRINGS_H)

# HACKY HACK
IF(APPLE)
  SET(HAVE_ISFINITE ON)
ENDIF(APPLE)


# others
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/epm-header.in ${CMAKE_CURRENT_BINARY_DIR}/epm-header)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/include/version.h.in ${CMAKE_CURRENT_BINARY_DIR}/version.h)


# to fix a bug with diverging N3
#STRING(REGEX REPLACE "-O3" "-O2" CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
#STRING(REGEX REPLACE "-O3" "-O2" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
# 
# STRING(REGEX REPLACE "-O2" "-O1" CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
# STRING(REGEX REPLACE "-O2" "-O1" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})
# 
# STRING(REGEX REPLACE "-O1" "-O0" CMAKE_C_FLAGS_RELEASE ${CMAKE_C_FLAGS_RELEASE})
# STRING(REGEX REPLACE "-O1" "-O0" CMAKE_CXX_FLAGS_RELEASE ${CMAKE_CXX_FLAGS_RELEASE})



INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_BINARY_DIR} ${CMAKE_SOURCE_DIR}/include )

ADD_LIBRARY(mincprog STATIC src/MincProg/print_version.c)

LINK_LIBRARIES(mincprog ${EBTKS_LIBRARIES} ${EBTKS_LAPACK_LIBRARIES} ${VOLUME_IO_LIBRARIES} ${LIBMINC_LIBRARIES})

ADD_EXECUTABLE(correct_field src/CorrectField/correctField.cc)

ADD_EXECUTABLE(evaluate_field  src/EvaluateField/evaluateField.cc 
        src/EvaluateField/evaluateFieldArgs.cc 
        src/SplineSmooth/fieldIO.cc 
        src/Splines/Spline.cc 
        src/Splines/TBSpline.cc)
        
ADD_EXECUTABLE(extracttag  src/ExtractTag/extracttag.c)

ADD_EXECUTABLE(sharpen_hist  
        src/SharpenHist/sharpen_hist.cc 
        src/SharpenHist/args.cc )

ADD_EXECUTABLE(spline_smooth  
        src/SplineSmooth/splineSmooth.cc 
        src/SplineSmooth/splineSmoothArgs.cc 
        src/SplineSmooth/fieldIO.cc 
        src/Splines/Spline.cc 
        src/Splines/TBSpline.cc)
        
ADD_EXECUTABLE(volume_hist  
        src/VolumeHist/minchist.cc 
        src/VolumeHist/DHistogram.cc
        src/VolumeHist/WHistogram.cc
        src/VolumeHist/GHistogram.cc
        src/VolumeHist/args.cc)

ADD_EXECUTABLE(volume_stats  
        src/VolumeStats/volumeStats.cc 
        src/VolumeStats/VolumeStatsArgs.cc)



# The in-memory C++ reimplementation of the Perl drivers.  Built as an OBJECT
# library so that both nu_correct_cxx and the per-routine tests link the same
# compilation.
# sharpen_hist.cc keeps its algorithm inside main(), so it is compiled with
# both main() and its args class renamed and linked for the three free
# functions it does export.  args.cc comes along only to satisfy the
# references that unused main() makes; renaming the class as well is what
# lets a second such translation unit coexist -- SharpenHist, SplineSmooth and
# VolumeHist each define a different global `class args`.
ADD_LIBRARY(sharpenhist_core OBJECT
        src/SharpenHist/sharpen_hist.cc
        src/SharpenHist/args.cc)
TARGET_COMPILE_DEFINITIONS(sharpenhist_core PRIVATE
        main=n3_sharpen_hist_main args=n3_sharpen_hist_args)

# splineSmooth.cc carries volume_domain, reduced_domain and the two fit loops,
# and a main() and an args class alongside them.  Same treatment.
ADD_LIBRARY(splinesmooth_core OBJECT
        src/SplineSmooth/splineSmooth.cc
        src/SplineSmooth/splineSmoothArgs.cc)
TARGET_COMPILE_DEFINITIONS(splinesmooth_core PRIVATE
        main=n3_spline_smooth_main args=n3_spline_smooth_args)

ADD_LIBRARY(n3pipeline_core OBJECT
        src/N3Pipeline/Buffers.cc
        src/N3Pipeline/Histogram.cc
        src/N3Pipeline/Sharpen.cc
        src/N3Pipeline/MincTools.cc
        src/N3Pipeline/FitField.cc
        src/N3Pipeline/SmoothField.cc
        src/N3Pipeline/Stopping.cc
        src/N3Pipeline/NuEstimate.cc
        src/N3Pipeline/NuEvaluate.cc
        src/SplineSmooth/fieldIO.cc
        src/Splines/Spline.cc
        src/Splines/TBSpline.cc
        src/VolumeHist/DHistogram.cc
        src/VolumeHist/WHistogram.cc
        src/VolumeHist/GHistogram.cc)

# The in-memory driver: nu_correct and nu_estimate collapsed into one binary,
# chosen by argv[0] (an installed name of nu_estimate*) or -estimate_only.  The
# two installed names are the same source, so the estimate-only program is a
# second executable target rather than a rename.
ADD_EXECUTABLE(nu_correct_cxx src/N3Pipeline/nu_correct_cxx.cc
  $<TARGET_OBJECTS:n3pipeline_core> $<TARGET_OBJECTS:sharpenhist_core> $<TARGET_OBJECTS:splinesmooth_core>)
ADD_EXECUTABLE(nu_estimate_cxx src/N3Pipeline/nu_correct_cxx.cc
  $<TARGET_OBJECTS:n3pipeline_core> $<TARGET_OBJECTS:sharpenhist_core> $<TARGET_OBJECTS:splinesmooth_core>)


SET(PERL ${PERL_EXECUTABLE})


CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/NUcorrect/field2imp.in              ${CMAKE_CURRENT_BINARY_DIR}/field2imp @ONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/NUcorrect/imp2field.in              ${CMAKE_CURRENT_BINARY_DIR}/imp2field @ONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/NUcorrect/make_template.in          ${CMAKE_CURRENT_BINARY_DIR}/make_template @ONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/NUcorrect/nu_estimate.in            ${CMAKE_CURRENT_BINARY_DIR}/nu_estimate @ONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/NUcorrect/nu_estimate.in            ${CMAKE_CURRENT_BINARY_DIR}/nu_correct @ONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/NUcorrect/nu_estimate_np_and_em.in  ${CMAKE_CURRENT_BINARY_DIR}/nu_estimate_np_and_em @ONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/NUcorrect/nu_evaluate.in      ${CMAKE_CURRENT_BINARY_DIR}/nu_evaluate @ONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/NUcorrect/resample_labels.in  ${CMAKE_CURRENT_BINARY_DIR}/resample_labels @ONLY)
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/NUcorrect/sharpen_volume.in   ${CMAKE_CURRENT_BINARY_DIR}/sharpen_volume @ONLY)

INSTALL(TARGETS 
 correct_field evaluate_field extracttag sharpen_hist spline_smooth volume_hist volume_stats
 nu_correct_cxx nu_estimate_cxx
  DESTINATION bin)

INSTALL(PROGRAMS
  ${CMAKE_CURRENT_BINARY_DIR}/field2imp
  ${CMAKE_CURRENT_BINARY_DIR}/imp2field
  ${CMAKE_CURRENT_BINARY_DIR}/make_template
  ${CMAKE_CURRENT_BINARY_DIR}/nu_estimate
  ${CMAKE_CURRENT_BINARY_DIR}/nu_correct
  ${CMAKE_CURRENT_BINARY_DIR}/nu_estimate_np_and_em
  ${CMAKE_CURRENT_BINARY_DIR}/nu_evaluate
  ${CMAKE_CURRENT_BINARY_DIR}/resample_labels
  ${CMAKE_CURRENT_BINARY_DIR}/sharpen_volume
  PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ  WORLD_EXECUTE WORLD_READ
   DESTINATION bin )

#install masks
INSTALL(FILES
    model_data/N3/average_305_mask_1mm.mnc.gz
    model_data/N3/icbm_avg_152_t1_tal_nlin_symmetric_VI_mask.mnc.gz
  DESTINATION 
    share/N3 )

# testing
IF(BUILD_TESTING)
  ADD_SUBDIRECTORY( testing )
ENDIF(BUILD_TESTING)
