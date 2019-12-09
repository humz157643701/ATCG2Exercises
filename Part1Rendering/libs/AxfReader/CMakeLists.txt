project(AxfReader)
cmake_minimum_required(VERSION 3.5)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# hope that shit works
find_package(HDF5)

add_executable(AxfTest AxfTest.cpp)

target_include_directories(AxfTest PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/HighFive/include/
    ${HDF5_INCLUDE_DIRS}
)

target_link_libraries(AxfTest
    ${HDF5_LIBRARIES}
)

# add interface library for easy integration
add_library(AxfReader INTERFACE)

target_include_directories(AxfReader INTERFACE
    ${HDF5_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/HighFive/include/
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(AxfReader INTERFACE
    ${HDF5_LIBRARIES}
)
