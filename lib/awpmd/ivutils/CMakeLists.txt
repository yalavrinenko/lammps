project(ivutils)
set(CMAKE_CXX_STANDARD 11)
include(sources.cmake)

get_filename_component(IVUTILS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include ABSOLUTE)
message("IvUtils include dir:" ${IVUTILS_INCLUDE_DIR})

if(NOT ${HAVE_ERF})   
  set(ERF_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/src/erf_namd.c" CACHE FILEPATH "Source for erf implementation")
  set(IVUTILS_SOURCES ${IVUTILS_SOURCES} ${ERF_SOURCE})
endif()

include_directories(${IVUTILS_INCLUDE_DIR})
add_library(ivutils ${BUILDCONF} ${IVUTILS_SOURCES} ${IVUTILS_HEADERS})

if(WIN32)
  if(MSVC_VERSION GREATER_EQUAL 1900)
    set(HAVE_MATH ON)
    set(HAVE_ERF ON)
    target_compile_definitions(ivutils PRIVATE -DHAVE_ERF) # don't include erf.h
  else()
    set(HAVE_MATH OFF)
    set(HAVE_ERF OFF)
  endif()
  target_compile_definitions(ivutils PRIVATE -DFU=0)

  target_compile_definitions(ivutils PRIVATE "-D_CRT_SECURE_NO_WARNINGS")
  target_compile_definitions(ivutils PRIVATE "-D_SCL_SECURE_NO_WARNINGS")
  target_compile_definitions(ivutils PRIVATE "-D_CRT_SECURE_NO_DEPRECATE")
  target_compile_definitions(ivutils PRIVATE "-D_USE_MATH_DEFINES")
    

  if(${HAVE_MATH})
    target_compile_definitions(ivutils PRIVATE -DNO_CMNMATH)
  endif()
else(NOT WIN32)
  target_compile_definitions(ivutils PRIVATE "-DUNIX=1")
  target_compile_definitions(ivutils PRIVATE "-DLINFO=1")
  target_compile_definitions(ivutils PRIVATE "-DFU=1")
endif(WIN32)

