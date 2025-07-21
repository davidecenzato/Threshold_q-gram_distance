## Add sdsl-lite
set(sdsl_lite_dir ${CMAKE_CURRENT_LIST_DIR}/sdsl-lite)
set(sdsl_lite_file ${sdsl_lite_dir}/CMakeLists.txt)
if(NOT EXISTS ${sdsl_lite_file})
    EXECUTE_PROCESS(COMMAND git submodule update --init external/sdsl-lite
                    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
                    OUTPUT_QUIET
            )
endif(NOT EXISTS ${sdsl_lite_file})

add_subdirectory(sdsl-lite)

## Add robin-map
set(robin_map_dir ${CMAKE_CURRENT_LIST_DIR}/robin-map)
set(robin_map_file ${robin_map_dir}/CMakeLists.txt)
if(NOT EXISTS ${robin_map_file})
    EXECUTE_PROCESS(COMMAND git submodule update --init external/robin-map
                    WORKING_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/..
                    OUTPUT_QUIET
            )
endif(NOT EXISTS ${robin_map_file})

add_subdirectory(robin-map)
