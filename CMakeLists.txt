cmake_minimum_required(VERSION 3.12)
project(krbalancing VERSION 0.0.5)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_BINARY_DIR})
list(APPEND CMAKE_PREFIX_PATH ${CMAKE_BINARY_DIR})

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

include(${CMAKE_SOURCE_DIR}/cmake/conan.cmake)

find_package(Python3 COMPONENTS Interpreter Development REQUIRED)

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


Python3_add_library(krbalancing MODULE WITH_SOABI src/krbalancing.cpp)
target_link_libraries(krbalancing PRIVATE Eigen3::Eigen CONAN_PKG::pybind11)
set_target_properties(krbalancing PROPERTIES
    INTERPROCEDURAL_OPTIMIZATION ON
    VISIBLITY_INLINES_HIDDEN ON)

configure_file(${CMAKE_SOURCE_DIR}/setup.py.in ${CMAKE_BINARY_DIR}/setup.py)

add_custom_target(wheel ALL 
    DEPENDS krbalancing
    COMMAND ${Python3_EXECUTABLE} setup.py bdist_wheel -d ${CMAKE_BINARY_DIR}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Generating wheel")
