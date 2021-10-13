cmake_minimum_required(VERSION 3.10)

project(ICELib)

# 包含目录
include_directories("D:\\Softwares\\Coding\\pybind11\\include")

# add_subdirectory("D:\\Softwares\\Coding\\pybind11")

include_directories("C:\\Users\\WenpingXie\\anaconda3\\include")

# 链接到静态库
link_directories("C:\\Users\\WenpingXie\\anaconda3\\libs")

# pybind11_add_module()

add_library(ICELib SHARED main.cpp Cylinder.cpp ArrayTable.cpp ParameterTable.cpp GasProperty.cpp Functions.cpp Valve.cpp)