cmake_minimum_required(VERSION 2.8)
project(rotation2D)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -std=c++17")
FIND_PACKAGE(DGtal REQUIRED 0.9.1)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

set(SOURCE_FILES src/main.cpp src/tools.cpp src/images.cpp src/rotations.cpp src/digital_objects.cpp)
add_executable(rotation2D ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(rotation2D ${DGTAL_LIBRARIES})
