cmake_minimum_required (VERSION 3.1)
project(tools-piccante)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11)

# make sure that the default is a RELEASE
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE RELEASE)
endif()

if(MSVC)
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  add_definitions(-D_USE_MATH_DEFINES)
endif()

find_package(MPI REQUIRED)

set(CMAKE_CXX_COMPILE_FLAGS "${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS}")

include_directories(${MPI_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/src/)

add_executable(newReader        ${CMAKE_SOURCE_DIR}/src/newReader.cpp ${CMAKE_SOURCE_DIR}/src/utilities-tools.cpp)
add_executable(WangularSpectrum ${CMAKE_SOURCE_DIR}/src/WangularSpectrum.cpp ${CMAKE_SOURCE_DIR}/src/utilities-tools.cpp)
add_executable(WenergySpectrum  ${CMAKE_SOURCE_DIR}/src/WenergySpectrum.cpp ${CMAKE_SOURCE_DIR}/src/utilities-tools.cpp)
add_executable(WsliceDaBinario  ${CMAKE_SOURCE_DIR}/src/WsliceDaBinario.cpp ${CMAKE_SOURCE_DIR}/src/utilities-tools.cpp)
add_executable(WspacialSpectrum ${CMAKE_SOURCE_DIR}/src/WspacialSpectrum.cpp ${CMAKE_SOURCE_DIR}/src/utilities-tools.cpp)
add_executable(lightReader      ${CMAKE_SOURCE_DIR}/src/lightReader.cpp ${CMAKE_SOURCE_DIR}/src/utilities-tools.cpp)
add_executable(frogReader       ${CMAKE_SOURCE_DIR}/src/frogReader.cpp ${CMAKE_SOURCE_DIR}/src/utilities-tools.cpp)
add_executable(titan            ${CMAKE_SOURCE_DIR}/src/titan.cpp ${CMAKE_SOURCE_DIR}/src/utilities-tools.cpp)
add_executable(titan2           ${CMAKE_SOURCE_DIR}/src/titan2.cpp ${CMAKE_SOURCE_DIR}/src/utilities-tools.cpp)
add_executable(clusters         ${CMAKE_SOURCE_DIR}/src/clusters.cpp ${CMAKE_SOURCE_DIR}/src/utilities-tools.cpp)
add_executable(randomballs      ${CMAKE_SOURCE_DIR}/src/randomballs.cpp ${CMAKE_SOURCE_DIR}/src/utilities-tools.cpp)
add_executable(multiFrogReader  ${CMAKE_SOURCE_DIR}/src/multiFrogReader.cpp ${CMAKE_SOURCE_DIR}/src/utilities-tools.cpp)
add_executable(MPItitan         ${CMAKE_SOURCE_DIR}/src/MPItitan.cpp)
add_executable(MPIpicker        ${CMAKE_SOURCE_DIR}/src/MPIpicker.cpp)

target_link_libraries(MPItitan  ${MPI_LIBRARIES}) 
target_link_libraries(MPIpicker ${MPI_LIBRARIES}) 

install (TARGETS 
  newReader
  WangularSpectrum
  WenergySpectrum
  WsliceDaBinario
  WspacialSpectrum
  lightReader
  frogReader
  titan
  titan2
  clusters
  randomballs
  multiFrogReader
  MPItitan
  MPIpicker 
  DESTINATION ${CMAKE_SOURCE_DIR}/bin/
)

