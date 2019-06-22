cmake_minimum_required(VERSION 3.14)
project(libfemlabgmsh C)

set(CMAKE_C_STANDARD 99)

include_directories(/run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/gmsh-4.3.0-Linux64/include)

add_executable(libfemlabgmsh main.c libfemlab.c libfemlab.h)

target_link_libraries(libfemlabgmsh /run/media/tobiaslafer/shared/Documents/Uni/BAK-Arbeit/gmsh-4.3.0-Linux64/lib/libgmsh.so)