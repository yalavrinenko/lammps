file(GLOB_RECURSE AWPMD_HEADERS RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "include/*.h")

set(AWPMD_SOURCES "src/box_hamiltonian.cpp" "src/rebuild_wp.cpp" "src/wpmd.cpp"  
                    "src/wpmd_split.cpp") 

if(AWPMD_USE_HOOMD_FORMAT)
  set(AWPMD_SOURCES ${AWPMD_SOURCES} "src/state_io.cpp")
endif()


#if(BUILD_DFT_WPMD)
#set(AWPMD_SOURCES ${AWPMD_SOURCES}
#    "src/box_hamiltonian.cpp" "src/rebuild_wp.cpp" "src/state_io.cpp" "src/wpmd.cpp"  
#                    "src/wpmd_split.cpp") 
#endif()
