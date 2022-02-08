find_package(MKL REQUIRED)
#set(CMAKE_CXX_STANDARD 14)

if(WIN32)
    if(MSVC_VERSION GREATER_EQUAL 1900)
        set(HAVE_MATH ON)
        set(HAVE_ERF ON)
        target_compile_definitions(lammps PRIVATE -DHAVE_ERF) # don't include erf.h
    else()
        set(HAVE_MATH OFF)
        set(HAVE_ERF OFF)
    endif()
    target_compile_definitions(lammps PRIVATE -DFU=0)
    
    target_compile_definitions(lammps PRIVATE "-D_CRT_SECURE_NO_WARNINGS")
    target_compile_definitions(lammps PRIVATE "-D_SCL_SECURE_NO_WARNINGS")
    target_compile_definitions(lammps PRIVATE "-D_CRT_SECURE_NO_DEPRECATE")
    target_compile_definitions(lammps PRIVATE "-D_USE_MATH_DEFINES")
    
    if(${HAVE_MATH})
        target_compile_definitions(lammps PRIVATE -DNO_CMNMATH)
    endif()
else(NOT WIN32)
    target_compile_definitions(lammps PRIVATE "-DUNIX=1")
    target_compile_definitions(lammps PRIVATE "-DLINFO=1")
    target_compile_definitions(lammps PRIVATE "-DFU=1")
endif(WIN32)


target_compile_definitions(lammps PRIVATE "-DUSE_AWPMD=1")

target_include_directories(lammps PRIVATE "${LAMMPS_LIB_SOURCE_DIR}/awpmd/awpmd/include")
target_include_directories(lammps PRIVATE "${LAMMPS_LIB_SOURCE_DIR}/awpmd/ivutils/include")
include_directories(${MKL_INCLUDE_DIR})

add_subdirectory("${LAMMPS_LIB_SOURCE_DIR}/awpmd/awpmd" ${LAMMPS_LIB_BINARY_DIR}/awpmd)
add_subdirectory("${LAMMPS_LIB_SOURCE_DIR}/awpmd/ivutils" ${LAMMPS_LIB_BINARY_DIR}/ivutils)

if (PKG_WPMD-DFT)
    set(wpmd_dft_source_dir "${LAMMPS_SOURCE_DIR}/AWPMD/WPMD-DFT")
    file(GLOB wpmd_dft_sources ${wpmd_dft_source_dir}/[^.]*.cpp)
    file(GLOB wpmd_dft_headers ${wpmd_dft_source_dir}/[^.]*.h)
    target_sources(lammps PRIVATE ${wpmd_dft_sources})

    RegisterStyles(${wpmd_dft_source_dir})

    target_include_directories(lammps PRIVATE "${LAMMPS_LIB_SOURCE_DIR}/awpmd/awpmd/awpmd-dft" ${wpmd_dft_source_dir})
    add_subdirectory(${LAMMPS_LIB_SOURCE_DIR}/awpmd/awpmd/awpmd-dft ${CMAKE_CURRENT_BINARY_DIR}/wpmd-dft-bin)
    target_link_libraries(lammps PRIVATE wpmd_dft_unit_cpu)
endif()

if (PKG_WPMD-NVGPU-DFT)
    set(wpmd_nvgpu_dft_source_dir "${LAMMPS_SOURCE_DIR}/AWPMD/WPMD-DFT/NVGPU")
    file(GLOB wpmd_nvgpu_dft_sources ${wpmd_nvgpu_dft_source_dir}/[^.]*.cpp)
    file(GLOB wpmd_nvgpu_dft_headers ${wpmd_nvgpu_dft_source_dir}/[^.]*.h)
    target_sources(lammps PRIVATE ${wpmd_nvgpu_dft_sources})

    RegisterStyles(${wpmd_nvgpu_dft_source_dir})

    target_include_directories(lammps PRIVATE "${LAMMPS_LIB_SOURCE_DIR}/awpmd/awpmd/awpmd-dft" ${wpmd_nvgpu_dft_source_dir})
    add_subdirectory(${LAMMPS_LIB_SOURCE_DIR}/awpmd/awpmd/awpmd-dft/nvgpu ${CMAKE_CURRENT_BINARY_DIR}/wpmd-nvgpu-dft-bin)
    target_link_libraries(lammps PRIVATE wpmd_dft_unit_nvgpu)
endif()

target_link_libraries(lammps PRIVATE awpmd ivutils ${MKL_LIBRARIES})