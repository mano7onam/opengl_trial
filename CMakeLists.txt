cmake_minimum_required(VERSION 3.14)
project(OpenGL)

set(CMAKE_CXX_STANDARD 14)

add_executable(OpenGL main1.cpp ply_loader.h)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

target_link_libraries(OpenGL ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )