cmake_minimum_required(VERSION 3.8)
project(awpmd)

include(sources.cmake)
get_filename_component(AWPMD_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include" ABSOLUTE)
message("AWPMD include dir:" ${AWPMD_INCLUDE_DIR})

get_filename_component(IVUTILS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../ivutils/include ABSOLUTE)
message("IVUTILS include dir:" ${IVUTILS_INCLUDE_DIR})

include_directories(${AWPMD_INCLUDE_DIR})
include_directories(${IVUTILS_INCLUDE_DIR})

add_library(awpmd ${BUILDCONF} ${AWPMD_SOURCES} ${AWPMD_HEADERS})