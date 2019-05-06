cmake_minimum_required(VERSION 3.3.0)
project(rotation3D)

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

find_package(OpenGL REQUIRED)
find_package(Qt5 COMPONENTS Core Gui OpenGL Xml Widgets REQUIRED)
FIND_PACKAGE(DGtal REQUIRED 0.9.1)
set(QT_USE_QTOPENGL TRUE)

find_path(QGLVIEWER_INCLUDE_DIR qglviewer.h
    /usr/include/QGLViewer
    /opt/local/include/QGLViewer
    /usr/local/include/QGLViewer
    /sw/include/QGLViewer
    ENV QGLVIEWERROOT
)

find_library(QGLVIEWER_LIBRARY_RELEASE
  NAMES qglviewer QGLViewer qglviewer-qt5 QGLViewer-qt5
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
  NAMES dqglviewer dQGLViewer dqglviewer-qt5 dQGLViewer-qt5 QGLViewerd2
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
include_directories(${QT_QTOPENGL_INCLUDE_DIR} ${OPENGL_INCLUDE_DIR})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

set(SOURCE_FILES src/main.cpp)
add_executable(rotation3D ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(rotation3D ${DGTAL_LIBRARIES} Qt5::OpenGL ${OpenGL_LIBRARIES}  Qt5::Gui Qt5::Xml)
