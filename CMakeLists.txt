# CMakeFiles.txt for N3
#
# Vladimir S. FONOV - vladimir.fonov@gmail.com

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

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
ENDIF(NOT MINC_TOOLKIT_BUILD)

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

LINK_LIBRARIES(mincprog ${EBTKS_LIBRARIES} ${VOLUME_IO_LIBRARIES} ${LIBMINC_LIBRARIES})

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
        src/VolumeHist/args.cc)

ADD_EXECUTABLE(volume_stats  
        src/VolumeStats/volumeStats.cc 
        src/VolumeStats/VolumeStatsArgs.cc)



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
