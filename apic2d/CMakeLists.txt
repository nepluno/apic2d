# APIC2D Executable

append_files (Headers "h" .)
append_files (Sources "cpp" .)
append_files (Sources "c" .)

if(CMAKE_CXX_COMPILER MATCHES "/em\\+\\+(-[a-zA-Z0-9.])?$")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -s LEGACY_GL_EMULATION=1 -s GL_FFP_ONLY=1")
  set(CMAKE_EXECUTABLE_SUFFIX ".html")
else()
  # Locate OpenGL
  find_package (OpenGL REQUIRED)
  if (OPENGL_FOUND)
    include_directories (${OPENGL_INCLUDE_DIR})
    set (APIC2D_LIBRARIES ${APIC2D_LIBRARIES} ${OPENGL_LIBRARIES})
  else (OPENGL_FOUND)
    message (SEND_ERROR "Unable to locate OpenGL")
  endif (OPENGL_FOUND)

  # Locate GLUT
  find_package (FreeGLUT REQUIRED)
  if (FREEGLUT_FOUND)
    include_directories (${FREEGLUT_INCLUDE_DIR})
    set (APIC2D_LIBRARIES ${APIC2D_LIBRARIES} ${FREEGLUT_LIBRARIES})
  else ()
    message (SEND_ERROR "Unable to locate FreeGLUT")
  endif ()
endif()

add_executable (APIC2D ${Headers} ${Templates} ${Sources})
set_property(TARGET APIC2D PROPERTY CXX_STANDARD 14)
if (WIN32)
  add_dependencies (APIC2D DllCopy)
endif ()
target_link_libraries (APIC2D ${APIC2D_LIBRARIES} )
