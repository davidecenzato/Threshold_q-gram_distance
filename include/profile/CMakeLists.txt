set(PROFILE_SOURCES thr_qgram_profile.hpp)

add_library(profile OBJECT ${PROFILE_SOURCES})
target_link_libraries(profile robin_map dictionary)