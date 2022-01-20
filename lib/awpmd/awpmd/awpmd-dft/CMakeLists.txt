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

option(DEBUG_WAVEPACKET "Use dummy wavepacket for debug" OFF)
if (DEBUG_WAVEPACKET)
    message(WARNING "USE Dummy Wavepacket class")
    target_compile_definitions(wpmd_dft_unit_cpu PRIVATE -DSEPARATE_DEBUG)
else (DEBUG_WAVEPACKET)
    target_include_directories(wpmd_dft_unit_cpu PRIVATE "../../ivutils/include/")
endif (DEBUG_WAVEPACKET)


