# APIC2D Executable

append_files (Headers "h" .)
append_files (Sources "cpp" .)
append_files (Sources "c" .)

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

# Locate TBB
find_package (TBB REQUIRED)
if (TBB_FOUND)
  include_directories (${TBB_INCLUDE_DIRS})
  set (APIC2D_LIBRARIES ${APIC2D_LIBRARIES} ${TBB_LIBRARIES})
else (TBB_FOUND)
  message (SEND_ERROR "Unable to locate TBB")
endif (TBB_FOUND)

add_executable (APIC2D ${Headers} ${Templates} ${Sources})
set_property(TARGET APIC2D PROPERTY CXX_STANDARD 14)
if (WIN32)
  add_dependencies (APIC2D DllCopy)
endif ()
target_link_libraries (APIC2D ${APIC2D_LIBRARIES} )
