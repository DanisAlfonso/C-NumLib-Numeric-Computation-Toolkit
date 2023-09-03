cmake_minimum_required(VERSION 3.12)
project(YourScientificLibrary VERSION 1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include subdirectories
add_subdirectory(include)
add_subdirectory(examples)
