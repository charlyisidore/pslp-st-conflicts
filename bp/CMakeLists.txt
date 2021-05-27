cmake_minimum_required(VERSION 3.12)
project(pslp_bp)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

find_package(SCIP REQUIRED)
find_package(CPLEX REQUIRED)
find_package(nlohmann_json 3.2.0 REQUIRED)

if (MSVC)
    add_compile_options(/W4)
else()
    add_compile_options(-Wall -Wextra -Wfatal-errors -pedantic)
endif()

add_library(${PROJECT_NAME}_lib STATIC
            src/branchrule_samediff.cpp
            src/cons_callback_cycle.cpp
            src/conshdlr_cycle.cpp
            src/conshdlr_samediff.cpp
            src/heur_init.cpp
            src/pricer.cpp
            src/pricer_cplex.cpp
            src/pricer_greedy.cpp
            src/pricer_cplex_lazy.cpp
            src/pricer_cplex_lo.cpp
            src/pricer_cplex_nf.cpp
            src/pricer_cplex_pa.cpp
            src/pricer_mip.cpp
            src/pricer_mip_lazy.cpp
            src/probdata.cpp
            src/reader_json.cpp)
set_target_properties(${PROJECT_NAME}_lib PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_compile_features(${PROJECT_NAME}_lib PUBLIC cxx_std_17)
target_compile_options(${PROJECT_NAME}_lib PUBLIC -DIL_STD $<$<CONFIG:Debug>:-DSCIP_DEBUG -DLOG_ACTIVE_LEVEL=LOG_LEVEL_TRACE>)
target_include_directories(${PROJECT_NAME}_lib PUBLIC include ${SCIP_INCLUDE_DIRS} ${CPLEX_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME}_lib PUBLIC pslp_common_lib ${SCIP_LIBRARIES} ${CPLEX_LIBRARIES} nlohmann_json::nlohmann_json)

add_executable(${PROJECT_NAME}_exe src/main.cpp)
set_target_properties(${PROJECT_NAME}_exe PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME}_exe PUBLIC ${PROJECT_NAME}_lib)