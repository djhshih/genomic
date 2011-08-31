# Locate FTGL library
# This module defines
#  FTGL_FOUND
#  FTGL_LIBRARIES
#  FTGL_INCLUDE_DIR

FIND_PATH(FTGL_INCLUDE_DIR ftgl.h
	PATH_SUFFIXES include/FTGL include
	PATH
	/usr
	/usr/local
)

FIND_LIBRARY(FTGL_LIBRARY ftgl
	PATH_SUFFIXES lib64 lib
	PATHS
	/usr
	/usr/local
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FTGL DEFAULT_MSG FTGL_LIBRARY FTGL_INCLUDE_DIR)

MARK_AS_ADVANCED(FTGL_INCLUDE_DIR FTGL_LIBRARY)

