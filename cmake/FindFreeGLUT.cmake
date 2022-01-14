# Copyright (C) 2007-2009 LuaDist.
# Created by Peter Kapec <kapecp@gmail.com>
# Redistribution and use of this file is allowed according to the terms of the MIT license.
# For details see the COPYRIGHT file distributed with LuaDist.
#	Note:
#		Searching headers and libraries is very simple and is NOT as powerful as scripts
#		distributed with CMake, because LuaDist defines directories to search for.
#		Everyone is encouraged to contact the author with improvements. Maybe this file
#		becomes part of CMake distribution sometimes.

# - Find FreeGLUT
# Find the native FreeGLUT libraries.
#
#  FREEGLUT_LIBRARIES    - List of libraries when using FreeGLUT.
#  FREEGLUT_FOUND        - True if FreeGLUT found.


FIND_LIBRARY(FREEGLUT_LIBRARY NAMES glut
	PATHS
	/opt/homebrew/lib
	/usr/lib)

# Look for the library.
if(APPLE)
	FIND_PATH(FREEGLUT_INCLUDE_DIR NAMES glut.h
	    PATHS
	    ${FREEGLUT_LIBRARY}/Headers)
else()
	FIND_PATH(FREEGLUT_INCLUDE_DIR NAMES GL/glut.h
	    PATHS
	    /usr/include)
endif()


# Handle the QUIETLY and REQUIRED arguments and set FREEGLUT_FOUND to TRUE if all listed variables are TRUE.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(FREEGLUT DEFAULT_MSG FREEGLUT_LIBRARY FREEGLUT_INCLUDE_DIR)

# Copy the results to the output variables.
IF(FREEGLUT_FOUND)
	SET(FREEGLUT_LIBRARIES ${FREEGLUT_LIBRARY})
	SET(FREEGLUT_INCLUDE_DIRS ${FREEGLUT_INCLUDE_DIR})
ELSE(FREEGLUT_FOUND)
	SET(FREEGLUT_LIBRARIES)
	SET(FREEGLUT_INCLUDE_DIRS)
ENDIF(FREEGLUT_FOUND)

MARK_AS_ADVANCED(FREEGLUT_INCLUDE_DIRS FREEGLUT_LIBRARIES)
