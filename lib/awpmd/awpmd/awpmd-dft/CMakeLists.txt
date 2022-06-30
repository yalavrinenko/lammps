cmake_minimum_required(VERSION 3.4.0)
project(wpmd_dft_unit CXX)

set(CMAKE_CXX_STANDARD 14)

#find_package(OpenMP)
#if(OPENMP_FOUND)
#    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

set(DATA_SRC DataTypes.hpp utils/ExternalTypes.hpp utils/Logger.cpp xcfunctionals/ApproxList.hpp utils/CCP4_DMap.cpp)
set(XCAPPROX_SRC xcfunctionals/LDA.cpp xcfunctionals/ModLDA.cpp xcfunctionals/LSDA.cpp xcfunctionals/TableApproximation.cpp)
set(SOURCE_FILES_CPU awpmd-dft.hpp awpmd-dft-cpu.cpp
        ${DATA_SRC}
        utils/ExternalTypes.hpp
        utils/SpaceMesh.cpp utils/DerivativesFunction.cpp
        numeric/cpu_integrator.hpp)

add_library(wpmd_dft_unit_cpu ${BUILDCONF} ${SOURCE_FILES_CPU} ${XCAPPROX_SRC})

if(WIN32)
  if(${CMAKE_GENERATOR} MATCHES "Visual Studio 15 2017")
    set(HAVE_MATH ON)
    set(HAVE_ERF ON)
    target_compile_definitions(wpmd_dft_unit_cpu PRIVATE -DERF_H) # don't include erf.h
  else()
    set(HAVE_MATH OFF)
    set(HAVE_ERF OFF)
  endif()
  target_compile_definitions(wpmd_dft_unit_cpu PRIVATE -DFU=0)

  target_compile_definitions(wpmd_dft_unit_cpu PRIVATE "-D_CRT_SECURE_NO_WARNINGS")
  target_compile_definitions(wpmd_dft_unit_cpu PRIVATE "-D_SCL_SECURE_NO_WARNINGS")
  target_compile_definitions(wpmd_dft_unit_cpu PRIVATE "-D_CRT_SECURE_NO_DEPRECATE")
  target_compile_definitions(wpmd_dft_unit_cpu PRIVATE "-D_USE_MATH_DEFINES")
    

  if(${HAVE_MATH})
    target_compile_definitions(wpmd_dft_unit_cpu PRIVATE -DNO_CMNMATH)
  endif()
else(NOT WIN32)
  target_compile_definitions(wpmd_dft_unit_cpu PRIVATE "-DUNIX=1")
  target_compile_definitions(wpmd_dft_unit_cpu PRIVATE "-DLINFO=1")
  target_compile_definitions(wpmd_dft_unit_cpu PRIVATE "-DFU=1")
endif(WIN32)
