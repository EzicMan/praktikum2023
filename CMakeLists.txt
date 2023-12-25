cmake_minimum_required(VERSION 3.16)
project(versh)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_VERBOSE_MAKEFILE 1)
set(APP_VERSION_MAJOR 1)
set(APP_VERSION_MINOR 0)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_BUILD_TYPE Debug)

add_compile_definitions(CHECK_K)
add_compile_definitions(OUT_SVTU)
#add_compile_definitions(OUT_F)
add_compile_definitions(OUT_SIGMA)
#add_compile_definitions(SMART_COMPARE)
#add_compile_definitions(OUT_BIS)

set(APP_TARGET test)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

add_executable(${APP_TARGET} main.cpp)

target_link_libraries(
    ${APP_TARGET}
    Eigen3::Eigen
)