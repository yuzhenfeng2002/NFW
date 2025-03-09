cmake_minimum_required(VERSION 3.28)
project(BEAQuA)

set(CMAKE_CXX_STANDARD 17)

# Specify the path to Boost
set(BOOST_ROOT /opt/homebrew/opt/boost)
set(BOOST_INCLUDEDIR /opt/homebrew/opt/boost/include)

find_package(Boost REQUIRED)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
endif()

file(GLOB SOURCES "*.cpp" "*.h")
add_executable(BEAQuA ${SOURCES})
target_link_libraries(BEAQuA ${Boost_LIBRARIES})