cmake_minimum_required (VERSION 3.0.0)

if(CMAKE_CXX_COMPILER MATCHES "/em\\+\\+")
  SET (CMAKE_C_COMPILER_WORKS 1)
  SET (CMAKE_CXX_COMPILER_WORKS 1)
endif()

project (APIC2D)

set(CMAKE_INSTALL_PREFIX "${CMAKE_SOURCE_DIR}/install")

# Initialize the build type (Release, Debug, etc)
if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: Debug Release."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

add_definitions (-DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}")

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffast-math")

# Add warnings to the compiler flags
if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions (-DDEBUG)
endif (CMAKE_BUILD_TYPE MATCHES Debug)

if (CMAKE_BUILD_TYPE MATCHES Release)
  add_definitions (-DNDEBUG)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif (CMAKE_BUILD_TYPE MATCHES Release)

# Add directory with macros
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# don't build in the source directory
if ("${CMAKE_BINARY_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")
  message (SEND_ERROR "Do not build in the source directory.")
  message (FATAL_ERROR "Remove the created \"CMakeCache.txt\" file and the \"CMakeFiles\" directory, then create a build directory and call \"${CMAKE_COMMAND} <path to the sources>\".")
endif ("${CMAKE_BINARY_DIR}" STREQUAL "${CMAKE_SOURCE_DIR}")

# finds all files with a given extension
macro (append_files files ext)
  foreach (dir ${ARGN})
    file (GLOB _files "${dir}/*.${ext}")
    list (APPEND ${files} ${_files})
  endforeach (dir)
endmacro (append_files)

# Eigen3 library is required
include_directories (${PROJECT_SOURCE_DIR}/include/eigen)

if (WIN32)
  include(ExternalProject)
  find_package(Git REQUIRED)
  ExternalProject_Add(thirdparty_win64_common
    GIT_REPOSITORY https://github.com/nepluno/thirdparty_win64_common.git
    GIT_TAG apic2d
    PREFIX ${PROJECT_SOURCE_DIR}/thirdparty
    CONFIGURE_COMMAND ""
    BUILD_COMMAND ""
    INSTALL_COMMAND ""
  )
  add_custom_target(DllCopy ALL
    COMMAND ${CMAKE_COMMAND} -E copy_directory
      "${PROJECT_SOURCE_DIR}/thirdparty/src/thirdparty_win64_common/bin/"
      "${CMAKE_BINARY_DIR}/apic2d/$<CONFIG>/"
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    VERBATIM)
  add_dependencies (DllCopy thirdparty_win64_common)
endif ()

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory (apic2d)
