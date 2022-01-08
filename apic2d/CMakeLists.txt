# APIC2D Executable

append_files (Headers "h" . pcgsolver)
append_files (Sources "cpp" . pcgsolver)
append_files (Sources "c" . pcgsolver)

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

add_executable (APIC2D ${Headers} ${Templates} ${Sources})
if (WIN32)
  add_dependencies (APIC2D DllCopy)
endif ()
target_link_libraries (APIC2D ${APIC2D_LIBRARIES} )
