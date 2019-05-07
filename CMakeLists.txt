cmake_minimum_required(VERSION 3.13)
project(OpenGL)

set(CMAKE_CXX_STANDARD 14)

#add_executable(OpenGL main.cpp kdtree.h SetPairs.h DataStructure.h Triangulation.h DelaunayTriangulation.cpp Triangle.cpp Vector3D.cpp DSU.h)
add_executable(OpenGL main_nurbs.cpp)

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
include_directories( ${OPENGL_INCLUDE_DIRS}  ${GLUT_INCLUDE_DIRS} )

target_link_libraries(OpenGL ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} )