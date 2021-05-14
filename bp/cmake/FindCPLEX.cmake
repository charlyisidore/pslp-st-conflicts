# CPLEX_FOUND
# CPLEX_INCLUDE_DIRS
# CPLEX_LIBRARIES

find_path(CPLEX_INCLUDE_DIR
    NAMES ilcplex/cplex.h
    HINTS ${CPLEX_DIR} $ENV{CPLEX_DIR}
    PATH_SUFFIXES cplex/include/ilcplex cplex/include
    PATHS ENV CPATH
          ENV C_INCLUDE_PATH
          ENV C_PLUS_INCLUDE_PATH
          ENV INCLUDE_PATH)

find_path(CPLEX_CONCERT_INCLUDE_DIR
    NAMES ilconcert/iloenv.h
    HINTS ${CPLEX_DIR} $ENV{CPLEX_DIR}
    PATH_SUFFIXES concert/include/ilconcert concert/include
    PATHS ENV CPATH
          ENV C_INCLUDE_PATH
          ENV C_PLUS_INCLUDE_PATH
          ENV INCLUDE_PATH)

find_library(CPLEX_CONCERT_LIBRARY
    NAMES concert
    HINTS ${CPLEX_DIR} $ENV{CPLEX_DIR}
    PATH_SUFFIXES concert/lib/x86-64_linux/static_pic lib
    PATHS ENV LIBRARY_PATH
          ENV LD_LIBRARY_PATH)

find_library(CPLEX_ILOCPLEX_LIBRARY
    NAMES ilocplex
    HINTS ${CPLEX_DIR} $ENV{CPLEX_DIR}
    PATH_SUFFIXES cplex/lib/x86-64_linux/static_pic lib
    PATHS ENV LIBRARY_PATH
          ENV LD_LIBRARY_PATH)

find_library(CPLEX_LIBRARY
    NAMES cplex
    HINTS ${CPLEX_DIR} $ENV{CPLEX_DIR}
    PATH_SUFFIXES cplex/lib/x86-64_linux/static_pic lib
    PATHS ENV LIBRARY_PATH
          ENV LD_LIBRARY_PATH)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CPLEX DEFAULT_MSG CPLEX_LIBRARY CPLEX_INCLUDE_DIR CPLEX_ILOCPLEX_LIBRARY CPLEX_CONCERT_LIBRARY CPLEX_CONCERT_INCLUDE_DIR)

if(CPLEX_FOUND)
    set(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR} ${CPLEX_CONCERT_INCLUDE_DIR})
    set(CPLEX_LIBRARIES ${CPLEX_CONCERT_LIBRARY} ${CPLEX_ILOCPLEX_LIBRARY} ${CPLEX_LIBRARY} pthread dl m)
endif()