# Copyright 2022 Carlos Reyes
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

cmake_minimum_required(VERSION 3.21)
set(CMAKE_CXX_STANDARD 20)

project(Gyrus
  VERSION 0.1.1
  DESCRIPTION "A neural network library for C++"
  LANGUAGES CXX)

include(CTest)

# -----------------------------------------------------------------------------
file(GLOB LIB_HEADER_LIST CONFIGURE_DEPENDS "${Gyrus_SOURCE_DIR}/include/brainyguy/*.hpp")
file(GLOB LIB_SOURCE_LIST CONFIGURE_DEPENDS "${Gyrus_SOURCE_DIR}/library/*.cpp")

add_library(gyrus SHARED ${LIB_HEADER_LIST} ${LIB_SOURCE_LIST})

# -----------------------------------------------------------------------------
target_include_directories(gyrus
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${Gyrus_SOURCE_DIR}/include>
)

target_compile_options(gyrus PRIVATE -Werror)
target_compile_features(gyrus PRIVATE cxx_std_20)
target_link_libraries(gyrus PUBLIC m)

# -------------------------------------------------------------------
file(GLOB TEST_SOURCE_LIST CONFIGURE_DEPENDS "${Gyrus_SOURCE_DIR}/test/*.cpp")
add_executable(gyrus_test "${TEST_SOURCE_LIST}")

target_include_directories(gyrus_test
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${Gyrus_SOURCE_DIR}/include>
)

target_link_libraries(gyrus_test PUBLIC gyrus)

# -------------------------------------------------------------------
file(GLOB SIMPLE_SAMPLE_SOURCE_LIST CONFIGURE_DEPENDS "${Gyrus_SOURCE_DIR}/sample/simple.cpp")
add_executable(simple_sample "${SIMPLE_SAMPLE_SOURCE_LIST}")

target_include_directories(simple_sample
    PUBLIC
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${Gyrus_SOURCE_DIR}/include>
)

target_compile_features(simple_sample PRIVATE cxx_std_20)
target_link_libraries(simple_sample PUBLIC gyrus)
