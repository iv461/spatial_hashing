cmake_minimum_required(VERSION 3.3)
project(shrs)

add_compile_options(-std=c++11)

find_package(PCL REQUIRED)
find_package(Eigen3 REQUIRED)

link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})
include_directories(
 include
 include/shrs
 ${PCL_INCLUDE_DIRS}
)

add_executable(${PROJECT_NAME}_tests test/SHRSTest.cpp)
target_link_libraries(${PROJECT_NAME}_tests ${PCL_LIBRARIES})
