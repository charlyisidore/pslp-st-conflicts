cmake_minimum_required(VERSION 3.12)
project(pslp_pricers)

if (MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wfatal-errors -pedantic)
endif()

add_executable(${PROJECT_NAME}_exe src/main.cpp)
set_target_properties(${PROJECT_NAME}_exe PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_compile_features(${PROJECT_NAME}_exe PUBLIC cxx_std_17)
target_link_libraries(${PROJECT_NAME}_exe PUBLIC pslp_bp_lib)