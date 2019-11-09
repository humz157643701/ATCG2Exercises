cmake_minimum_required(VERSION 3.0)

project(external_dependencies)

## build and link all external dependencies to a single library target
add_library(atcg2p1_external_dependencies INTERFACE)
set(link_targets "")
set(dep_includes "")

## glm
# should make the interface library target "glm" available
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/glm-master")
list(APPEND link_targets glm)

## glfw
# build settings for glfw (keep the extra stuff to a minimum)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/glfw-3.2.1")
list(APPEND link_targets glfw ${GLFW_LIBRARIES})
list(APPEND dep_includes "${CMAKE_CURRENT_SOURCE_DIR}/glfw-3.2.1/include")

## glew
# build settings for glew
set(BUILD_UTILS OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

# should make a target named "glew_s" available
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/glew-2.1.0/build/cmake")
list(APPEND link_targets glew_s)
list(APPEND dep_includes "${CMAKE_CURRENT_SOURCE_DIR}/glew-2.1.0/include" "${CMAKE_CURRENT_SOURCE_DIR}/glew-2.1.0/include/GL")

## opengl system libs
if(NOT OPENGL_FOUND)
    find_package(OpenGL REQUIRED)
    if(NOT OPENGL_FOUND)
        message(FATAL_ERROR "OpenGL not found")
    else()
        list(APPEND link_targets ${OPENGL_LIBRARIES})
    endif()
else()
    message(STATUS "OpenGL is already included by GLEW or GLFW")
endif()

## stb
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/stb")
list(APPEND link_targets stb)

## opencv
#find_package(OpenCV REQUIRED)
#list(APPEND link_targets  ${OpenCV_LIBS})

## assimp
# note: set assimp_DIR to the install/lib/cmake directory to allow cmake to find assimps package config
# because assimp is retarded, you must append the .lib directories to your system path o.O
# find_package(assimp REQUIRED)
# list(APPEND link_targets ${ASSIMP_LIBRARIES})
# list(APPEND dep_includes ${ASSIMP_INCLUDE_DIRS})

## freetype
find_package(freetype CONFIG)
if(NOT freetype_FOUND)
    message(STATUS "freetype not found in CONFIG mode. Trying MODULE mode...")
    find_package(freetype REQUIRED)
endif()
list(APPEND link_targets freetype)

## AxfReader
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/AxfReader")
list(APPEND link_targets AxfReader)

target_link_libraries(atcg2p1_external_dependencies INTERFACE ${link_targets})
target_include_directories(atcg2p1_external_dependencies INTERFACE ${dep_includes})

