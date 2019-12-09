cmake_minimum_required(VERSION 3.0)

project(external_dependencies)

## build and link all external dependencies to a single library target
add_library(atcg2p2_external_dependencies INTERFACE)
set(link_targets "")
set(dep_includes "")

# libigl
#add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/libigl-master")
add_library(libigl INTERFACE)
target_include_directories(libigl INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/libigl-master/include")
list(APPEND link_targets libigl)

target_link_libraries(atcg2p2_external_dependencies INTERFACE ${link_targets})
target_include_directories(atcg2p2_external_dependencies INTERFACE ${dep_includes})

