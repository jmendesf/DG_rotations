cmake_minimum_required(VERSION 2.8)
project(rotation2D)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -O3 -std=c++17")
FIND_PACKAGE(DGtal REQUIRED 0.9.1)
#FIND_PACKAGE(DGtal REQUIRED 1.0)
INCLUDE_DIRECTORIES(${DGTAL_INCLUDE_DIRS})
LINK_DIRECTORIES(${DGTAL_LIBRARY_DIRS})

set(SOURCE_FILES src/main.cpp src/tools.cpp src/images.cpp src/rotations.cpp src/digital_objects.cpp)
add_executable(rotation2D ${SOURCE_FILES})
TARGET_LINK_LIBRARIES(rotation2D ${DGTAL_LIBRARIES})

# Fix GCC C++17 and Boost.ICL incompatibility (needed to build dynarmic)
# See https://bugzilla.redhat.com/show_bug.cgi?id=1485641#c1
if (CMAKE_COMPILER_IS_GNUCC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-new-ttp-matching")
endif()

