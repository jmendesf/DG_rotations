cmake_minimum_required(VERSION 3.3.0)
project(eval3D)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
set(OpenGL_GL_PREFERENCE "GLVND")

find_package(OpenGL REQUIRED)
find_package(Qt4 REQUIRED)
FIND_PACKAGE(DGtal REQUIRED 0.9.1)

SET(QT_USE_QTOPENGL TRUE)
INCLUDE(${QT_USE_FILE})
INCLUDE_DIRECTORIES(${QT_INCLUDES})

find_path(QGLVIEWER_INCLUDE_DIR qglviewer.h
    /usr/include/QGLViewer
    /opt/local/include/QGLViewer
    /usr/local/include/QGLViewer
    /sw/include/QGLViewer
    ENV QGLVIEWERROOT
)

find_library(QGLVIEWER_LIBRARY_RELEASE
  NAMES qglviewer QGLViewer qglviewer-qt4 QGLViewer-qt4
  PATHS /usr/lib
        /usr/local/lib
        /opt/local/lib
        /sw/lib
        ENV QGLVIEWERROOT
        ENV LD_LIBRARY_PATH
        ENV LIBRARY_PATH
  PATH_SUFFIXES QGLViewer QGLViewer/release
)

find_library(QGLVIEWER_LIBRARY_DEBUG
  NAMES dqglviewer dQGLViewer dqglviewer-qt4 dQGLViewer-qt4 QGLViewerd2
  PATHS /usr/lib
        /usr/local/lib
        /opt/local/lib
        /sw/lib
        ENV QGLVIEWERROOT
        ENV LD_LIBRARY_PATH
        ENV LIBRARY_PATH
  PATH_SUFFIXES QGLViewer QGLViewer/debug
)

if(QGLVIEWER_LIBRARY_RELEASE)
  if(QGLVIEWER_LIBRARY_DEBUG)
    set(QGLVIEWER_LIBRARY optimized ${QGLVIEWER_LIBRARY_RELEASE} debug ${QGLVIEWER_LIBRARY_DEBUG})
  else()
    set(QGLVIEWER_LIBRARY ${QGLVIEWER_LIBRARY_RELEASE})
  endif()
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QGLVIEWER DEFAULT_MSG QGLVIEWER_INCLUDE_DIR QGLVIEWER_LIBRARY)

ADD_DEFINITIONS(${QT_DEFINITIONS})
include_directories(${DGTAL_INCLUDE_DIRS})
include_directories(${QT_QTOPENGL_INCLUDE_DIR} ${OPENGL_INCLUDE})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

set(SOURCE_FILES src/main.cpp)
add_executable(eval3D ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(eval3D ${DGTAL_LIBRARIES} 
 		      ${QGLVIEWER_LIBRARY} 
 		      ${QT_QTXML_LIBRARY} 
   		      ${QT_QTOPENGL_LIBRARY} 
                      ${QT_QTGUI_LIBRARY} 
                      ${QT_QTCORE_LIBRARY} 
                      ${OPENGL_gl_LIBRARY} 
                      ${OPENGL_glu_LIBRARY})
