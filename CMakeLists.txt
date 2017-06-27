cmake_minimum_required(VERSION 3.2.3)
project (FlexiLib)

# Setup include directories
if (COMMAND cmake_policy)
    cmake_policy (SET CMP0003 NEW)
    cmake_policy (SET CMP0022 NEW)
    cmake_policy (SET CMP0026 OLD)
    cmake_policy (SET CMP0042 NEW)
endif ()

set (CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/CMake/Modules)
include (UrhoCommon)
include_directories (${CMAKE_SOURCE_DIR}/Source)

# Libraries
add_subdirectory (Source/CharacterAnimator)
add_subdirectory (Source/CharacterController)

# Tools
add_subdirectory (Source/CharacterAnimationTool)
add_subdirectory (Source/Player)
