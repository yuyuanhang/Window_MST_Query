cmake_minimum_required(VERSION 3.10)

project(MST)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-O3)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

# Add the executable target
add_executable(MST src/main.cpp src/graph_io.cpp src/graph.cpp src/mst.cpp src/window_generator.cpp src/quad_tree.h
        src/hybrid.h
        src/dynamic_mst.h
        src/dynamic_mst.cpp)

# set include directories
target_include_directories(MST PRIVATE src/)
target_link_libraries(MST stdc++fs)