cmake_minimum_required(VERSION 3.12.0)
project(wpmd_dft_unit_nvgpu CUDA CXX)

find_package(CUDA REQUIRED)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_ARCHITECTURES OFF)

set(CUDA_NVCC_FLAGS "${CMAKE_CXX_FLAGS} -rdc=true -c")

set(COMMON_SRC
    ../DataTypes.hpp
    ../utils/ExternalTypes.hpp
    ../utils/Logger.cpp
    ../xcfunctionals/ApproxList.hpp)

set(XCAPPROX_SRC
    ../xcfunctionals/LDA.cpp
    ../xcfunctionals/ModLDA.cpp
    ../xcfunctionals/LSDA.cpp
    ../xcfunctionals/TableApproximation.cpp)

set_source_files_properties(${XCAPPROX_SRC} PROPERTIES LANGUAGE CUDA)

set(SOURCE_FILES_GPU
        ../awpmd-dft.hpp
        awpmd-dft-nvgpu.hpp
        awpmd-dft-nvgpu.cu
        SpaceMesh.cu
        operators_overload.hpp
        operators_overload.cu
        xc_init_kernels.hpp
        nd4mesh.cu
        nd4mesh.hpp
        cuApproxTableObject.cu
        cuApproxTableObject.hpp
        gpu_integration.hpp
        ${DATA_SRC})


if (DFT_DISABLE_CUDA_TEXTURE)
    add_definitions(-DCUDA_NO_TEXTURE_OBJ=1)
endif()

#add_library(wpmd_dft_unit_nvgpu SHARED ${SOURCE_FILES_GPU} ${XCAPPROX_SRC})
add_library(wpmd_dft_unit_nvgpu SHARED ${SOURCE_FILES_GPU} ${XCAPPROX_SRC})

set_target_properties(wpmd_dft_unit_nvgpu PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
target_link_libraries(wpmd_dft_unit_nvgpu ${CUDA_LIBRARIES})

if(WIN32)
  target_compile_definitions(wpmd_dft_unit_nvgpu PRIVATE "-D_USE_MATH_DEFINES")
#  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
 set_target_properties(wpmd_dft_unit_nvgpu PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

