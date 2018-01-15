# - Find the PacBio BAM includes and library
#
# This module searches libpng, the library for working with PBBAM files.
#
# It defines the following variables
#  PBBAM_INCLUDE_DIRS, where to find png.h, etc.
#  PBBAM_LIBRARIES, the libraries to link against to use PBBAM.
#  PBBAM_DEFINITIONS - You should add_definitons(${PBBAM_DEFINITIONS}) before compiling code that includes png library files.
#  PBBAM_FOUND, If false, do not try to use PBBAM.
#  PBBAM_VERSION_STRING - the version of the PBBAM library found (since CMake 2.8.8)
# Also defined, but not for general use are
#  PBBAM_LIBRARY, where to find the PBBAM library.
# For backward compatiblity the variable PBBAM_INCLUDE_DIR is also set. It has the same value as PBBAM_INCLUDE_DIRS.
#
# Since PBBAM depends on the Zlib libraries, none of the above will be
# defined unless ZLib can be found.

if(PBBAM_FIND_QUIETLY)
  set(_FIND_ZLIB_ARG QUIET)
endif()
find_package(ZLIB ${_FIND_ZLIB_ARG})

if(ZLIB_FOUND)
  find_path(PBBAM_INCLUDE_DIR pbbam/PbiFile.h
		/usr/local/include/             # OpenBSD
  )

  find_library(PBBAM_LIBRARY NAMES pbbam )

  if (PBBAM_LIBRARY AND PBBAM_INCLUDE_DIR)
      # png.h includes zlib.h. Sigh.
      set(PBBAM_INCLUDE_DIRS ${PBBAM_INCLUDE_DIR} ${ZLIB_INCLUDE_DIR} )
      set(PBBAM_INCLUDE_DIR ${PBBAM_INCLUDE_DIRS} ) # for backward compatiblity
      set(PBBAM_LIBRARIES ${PBBAM_LIBRARY} ${ZLIB_LIBRARY})

#       if (CYGWIN)
#         if(BUILD_SHARED_LIBS)
#            # No need to define PBBAM_USE_DLL here, because it's default for Cygwin.
#         else()
#           set (PBBAM_DEFINITIONS -DPBBAM_STATIC)
#         endif()
#       endif ()

  endif ()

  if (PBBAM_INCLUDE_DIR AND EXISTS "${PBBAM_INCLUDE_DIR}/png.h")
      file(STRINGS "${PBBAM_INCLUDE_DIR}/png.h" png_version_str REGEX "^#define[ \t]+PBBAM_LIBPBBAM_VER_STRING[ \t]+\".+\"")

      string(REGEX REPLACE "^#define[ \t]+PBBAM_LIBPBBAM_VER_STRING[ \t]+\"([^\"]+)\".*" "\\1" PBBAM_VERSION_STRING "${png_version_str}")
      unset(png_version_str)
  endif ()
endif()

# handle the QUIETLY and REQUIRED arguments and set PBBAM_FOUND to TRUE if
# all listed variables are TRUE
include(${CMAKE_CURRENT_LIST_DIR}/FindPackageHandleStandardArgs.cmake)
find_package_handle_standard_args(PBBAM
                                  REQUIRED_VARS PBBAM_LIBRARY PBBAM_INCLUDE_DIR
                                  VERSION_VAR PBBAM_VERSION_STRING)

mark_as_advanced(PBBAM_INCLUDE_DIR PBBAM_LIBRARY )
