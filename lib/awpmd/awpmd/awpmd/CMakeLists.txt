project(awpmd)

set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")


option(AWPMD_BUILD_DFT  "Enable DFT WPMD" OFF)
option(AWPMD_BUILD_DFT_GPU  "Enable DFT GPU WPMD" OFF)
option(AWPMD_USE_HOOMD_FORMAT  "Enable HOOMD i/o file support (requires wxWidgets)" OFF)

find_package(MKL REQUIRED)  # use LAPACK check instead
include_directories(${MKL_INCLUDE_DIR})

if(AWPMD_USE_HOOMD_FORMAT)
  find_package(wxWidgets REQUIRED base xml)
  add_definitions("-DwxUSE_GUI=0")
  include_directories(${wxWidgets_INCLUDE_DIRS})
endif()


include(sources.cmake)
get_filename_component(AWPMD_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/include" ABSOLUTE)
#message("AWPMD include dir:" ${AWPMD_INCLUDE_DIR})

get_filename_component(IVUTILS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../ivutils/include ABSOLUTE)
#message("IVUTILS include dir:" ${IVUTILS_INCLUDE_DIR})

include_directories(${AWPMD_INCLUDE_DIR})
include_directories(${IVUTILS_INCLUDE_DIR})

add_library(awpmd ${BUILDCONF} ${AWPMD_SOURCES} ${AWPMD_HEADERS})
target_link_libraries(awpmd ${MKL_LIBRARIES})

if(AWPMD_BUILD_DFT)
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/awpmd-dft")
  add_definitions("-DDFT_EXTENSION=1")
  set(ISCPU ON)
  add_subdirectory(awpmd-dft  awpmd-dft_bin)
  target_link_libraries(awpmd awpmd_dft)
endif()

if(AWPMD_BUILD_DFT_GPU)
  include_directories("${CMAKE_CURRENT_SOURCE_DIR}/awpmd-dft")
  add_definitions("-DDFT_EXTENSION=1")
  set(ISCUDA ON)
  add_subdirectory(awpmd-dft  awpmd-dft_bin_gpu)
  target_link_libraries(awpmd awpmd_dft_gpu)
endif()

