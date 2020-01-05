cmake_minimum_required(VERSION 3.0)

project(external_dependencies)

## build and link all external dependencies to a single library target
add_library(atcg2p2_external_dependencies INTERFACE)
set(link_targets "")
set(dep_includes "")

# glew
## glew
# build settings for glew
set(BUILD_UTILS OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

# glad
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/glad")
list(APPEND link_targets glad)

# nanoflann
add_library(nanoflann INTERFACE)
target_include_directories(nanoflann INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/nanoflann-master/include")
list(APPEND link_targets nanoflann)

# should make a target named "glew_s" available
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/glew-2.1.0/build/cmake")
list(APPEND link_targets glew_s)
list(APPEND dep_includes "${CMAKE_CURRENT_SOURCE_DIR}/glew-2.1.0/include" "${CMAKE_CURRENT_SOURCE_DIR}/glew-2.1.0/include/GL")

# glfw
# build settings for glfw (keep the extra stuff to a minimum)
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)

add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/glfw-3.3")
list(APPEND link_targets glfw ${GLFW_LIBRARIES})
list(APPEND dep_includes "${CMAKE_CURRENT_SOURCE_DIR}/glfw-3.3/include")

# eigen
add_library(eigen INTERFACE)
target_include_directories(eigen INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/eigen")
list(APPEND link_targets eigen)

# libigl
add_library(libigl INTERFACE)
target_include_directories(libigl INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/libigl-master/include")
list(APPEND link_targets libigl)

target_link_libraries(atcg2p2_external_dependencies INTERFACE ${link_targets})
target_include_directories(atcg2p2_external_dependencies INTERFACE ${dep_includes})

