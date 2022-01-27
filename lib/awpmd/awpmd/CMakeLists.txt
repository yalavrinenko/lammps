cmake_minimum_required(VERSION 3.8)
project(awpmd)

#set(CMAKE_CXX_STANDARD 14)

include(sources.cmake)
get_filename_component(AWPMD_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include" ABSOLUTE)
message("AWPMD include dir:" ${AWPMD_INCLUDE_DIR})

get_filename_component(IVUTILS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../ivutils/include ABSOLUTE)
message("IVUTILS include dir:" ${IVUTILS_INCLUDE_DIR})

include_directories(${AWPMD_INCLUDE_DIR})
include_directories(${IVUTILS_INCLUDE_DIR})

add_library(awpmd ${BUILDCONF} ${AWPMD_SOURCES} ${AWPMD_HEADERS})

if(WIN32)
  if(${CMAKE_GENERATOR} MATCHES "Visual Studio 15 2017")
    set(HAVE_MATH ON)
    set(HAVE_ERF ON)
    target_compile_definitions(awpmd PRIVATE -DERF_H) # don't include erf.h
  else()
    set(HAVE_MATH OFF)
    set(HAVE_ERF OFF)
  endif()
  target_compile_definitions(awpmd PRIVATE -DFU=0)

  target_compile_definitions(awpmd PRIVATE "-D_CRT_SECURE_NO_WARNINGS")
  target_compile_definitions(awpmd PRIVATE "-D_SCL_SECURE_NO_WARNINGS")
  target_compile_definitions(awpmd PRIVATE "-D_CRT_SECURE_NO_DEPRECATE")
  target_compile_definitions(awpmd PRIVATE "-D_USE_MATH_DEFINES")
    

  if(${HAVE_MATH})
    target_compile_definitions(awpmd PRIVATE -DNO_CMNMATH)
  endif()
else(NOT WIN32)
  target_compile_definitions(awpmd PRIVATE "-DUNIX=1")
  target_compile_definitions(awpmd PRIVATE "-DLINFO=1")
  target_compile_definitions(awpmd PRIVATE "-DFU=1")
endif(WIN32)
