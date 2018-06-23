cmake_minimum_required (VERSION 3.1)
project (CleverLeaf CXX Fortran)
ENABLE_LANGUAGE(C)

set (CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake;${CMAKE_MODULE_PATH}")

find_package (SAMRAI REQUIRED)
find_package (HDF5 REQUIRED)
set (Boost_NO_BOOST_CMAKE ON)
find_package (Boost REQUIRED)
find_package (MPI REQUIRED)

include_directories (
  ${SAMRAI_INCLUDE_DIRS}
  ${HDF5_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${MPI_CXX_INCLUDE_PATH}
)

if (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model source -fp-model strict -prec-div -prec-sqrt")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ffloat-store")
endif ()

add_definitions(-DVERSION="git" -DHOST_NAME="cab")

set(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "")

add_subdirectory (src)
