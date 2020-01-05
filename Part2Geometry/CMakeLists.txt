cmake_minimum_required(VERSION 3.0)
project(ATCG2_Part2_Geometry)

set(SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
set(INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/src")

##-------------------------------add your source files here-------------------------------------------------------------
list(APPEND INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/src")
list(APPEND INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/include")

list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/include/icp.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/include/symmetry.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/include/meshsamplers.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/include/Octree.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/include/integral_invariant_signatures.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/include/mesh_saliency.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/include/mesh.h")

list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/icp.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Octree.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/mesh.cpp")

##--------------------------------build source groups for visual studio-------------------------------------------------

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCES})

##--------------------------------external dependencies-----------------------------------------------------------------
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/libs")

##--------------------------------executable target---------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 14)

add_executable(ATCG2P2Geometry ${SOURCES})
target_include_directories(
        ATCG2P2Geometry
        PRIVATE ${INCLUDES}
)

target_link_libraries(ATCG2P2Geometry PUBLIC atcg2p2_external_dependencies)

##-------------------------------copy assets to output------------------------------------------------------------------

#file(COPY "assets" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")