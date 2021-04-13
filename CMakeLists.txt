#
# NOTE: This should be called by setup.py file.
# 

cmake_minimum_required(VERSION 3.12)
project(krbalancing)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)
find_package(OpenMP REQUIRED)

include(${CMAKE_SOURCE_DIR}/cmake/conan.cmake)


execute_process(COMMAND ${Python3_EXECUTABLE} setup.py --version
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE CMAKE_PROJECT_VERSION
    OUTPUT_STRIP_TRAILING_WHITESPACE)
message(STATUS " project version ${CMAKE_PROJECT_VERSION}")
add_definitions(-DCMAKE_PROJECT_VERSION=\"${CMAKE_PROJECT_VERSION}\")

# pybind11 config file has a bug that makes is impossible to use with
# conan_cmake_configure.
# https://github.com/conan-io/conan-center-index/pull/4445
conan_cmake_configure(REQUIRES 
    # pybind11/2.6.1
    eigen/3.3.9
    GENERATORS cmake_find_package)
# find_package(pybind11 REQUIRED)

conan_cmake_autodetect(settings)
conan_cmake_install(PATH_OR_REFERENCE .
    BUILD missing
    REMOTE conan-center
    SETTINGS ${settings})

find_package(Eigen3 REQUIRED)


# hack around pybind11 bug 
# https://github.com/conan-io/conan-center-index/pull/4445
conan_cmake_run(REQUIRES pybind11/2.6.1 
    BASIC_SETUP NO_OUTPUT_DIRS CMAKE_TARGETS
    BUILD missing)


pybind11_add_module(krbalancing  src/krbalancing.cpp)
target_link_libraries(krbalancing PRIVATE Eigen3::Eigen CONAN_PKG::pybind11)
target_link_libraries(krbalancing PRIVATE OpenMP::OpenMP_CXX)
