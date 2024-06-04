cmake_minimum_required(VERSION 3.6)
project (mesh_imp)

file(GLOB SRC_FILES *.cpp *.h)
add_executable(mesh_imp ${SRC_FILES})
