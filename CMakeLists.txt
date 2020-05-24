cmake_minimum_required(VERSION 3.6)
project(npsv CXX)

include(ExternalProject)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# Set global properties
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Get external dependencies
find_package(LibLZMA REQUIRED)
find_package(ZLIB REQUIRED)
find_package(BZip2 REQUIRED)
find_package(SeqLib REQUIRED)

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

set(PYBIND11_CPP_STANDARD -std=c++14)
add_subdirectory(lib/pybind11)

# Set source directory
set(SOURCE_DIR "src/npsv")

include_directories(${SOURCE_DIR} ${SEQLIB_INCLUDE_DIRS})
set(SOURCES "${SOURCE_DIR}/aligner.cpp" "${SOURCE_DIR}/bindings.cpp")

# Add bindings here
pybind11_add_module(npsva ${SOURCES})

target_link_libraries(npsva
  PRIVATE
  ${SEQLIB_LIBRARIES}
  ${ZLIB_LIBRARIES}
  ${LIBLZMA_LIBRARIES}
  ${BZIP2_LIBRARIES}
  Threads::Threads
)