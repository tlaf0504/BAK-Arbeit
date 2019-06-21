cmake_minimum_required(VERSION 3.14)
project(gmsh_c C)

set(CMAKE_C_STANDARD 11)

include_directories(/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/gmsh-4.3.0-Linux64/include)
add_executable(gmsh_c main.c)

target_link_libraries(gmsh_c /run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/gmsh-4.3.0-Linux64/lib/libgmsh.so)




