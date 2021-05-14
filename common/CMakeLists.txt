cmake_minimum_required(VERSION 3.12)
project(pslp_common)

if (MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wfatal-errors -pedantic)
endif()

add_library(${PROJECT_NAME}_lib INTERFACE)
target_compile_features(${PROJECT_NAME}_lib INTERFACE cxx_std_17)
target_compile_options(${PROJECT_NAME}_lib INTERFACE -DIL_STD $<$<CONFIG:Debug>:-DSCIP_DEBUG -DLOG_ACTIVE_LEVEL=LOG_LEVEL_TRACE>)
target_include_directories(${PROJECT_NAME}_lib INTERFACE include)