cmake_minimum_required(VERSION 3.0)
project(ATCG2_Part1_Rendering)

set(SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
set(INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/src")

##-------------------------------add your source files here-------------------------------------------------------------
list(APPEND INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/src")

list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Camera.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Camera.cpp")
# list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/DepthPeelRenderer.h")
# list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/DepthPeelRenderer.cpp")
# list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/WeightedBlendedOITRenderer.h")
# list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/WeightedBlendedOITRenderer.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/IRenderer.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Material.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Material.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Mesh.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Mesh.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Model.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Model.cpp")
# list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/OpaqueRenderer.h")
# list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/OpaqueRenderer.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Primitives.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Primitives.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Scene.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Scene.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Window.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Window.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Texture.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Texture.cpp")
# list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Framebuffer.h")
# list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Framebuffer.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Shader.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Shader.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Transform.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Lights.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/Lights.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/SceneLoader.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/SceneLoader.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/WardRenderer.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/WardRenderer.cpp")


#linked list version of the algorithms
#list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/INDLLWBT.h")
#list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/INDLLWBT.cpp")

#list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/LinkedListWBTv2.h")
#list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/LinkedListWBTv2.cpp")

#layered rendering version of the algortihms
#list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/INDLRWBT.h")
#list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/INDLRWBT.cpp")

#particle stuff
#list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/ParticleEmitter.h")
#list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/ParticleEmitter.cpp")

#bspline class
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/BSpline.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/BSpline.cpp")

#perfdisplay
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/PerfDisplay.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/PerfDisplay.cpp")

##---------------------------------framework files----------------------------------------------------------------------
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/framework/CommonTypes.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/framework/fw_config.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/framework/GameWindow.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/framework/GameWindow.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/framework/glerror.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/framework/glerror.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/framework/Input.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/framework/Input.cpp")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/framework/libheaders.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/framework/OBJLoader.h")
list(APPEND SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/framework/OBJLoader.cpp")
# add that directory to include list:
list(APPEND INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/framework")

##--------------------------------build source groups for visual studio-------------------------------------------------

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SOURCES})

##--------------------------------external dependencies-----------------------------------------------------------------
add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/libs")

##--------------------------------executable target---------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 14)

add_executable(ATCG2P1Rendering ${SOURCES})
target_include_directories(
        ATCG2P1Rendering
        PRIVATE ${INCLUDES}
)

target_link_libraries(ATCG2P1Rendering PUBLIC atcg2p1_external_dependencies)

##-------------------------------copy assets to output------------------------------------------------------------------

file(COPY "assets" DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")