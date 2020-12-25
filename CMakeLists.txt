######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

cmake_minimum_required(VERSION 3.0)
project(overlap_analysis VERSION 0.0.0)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")

# Compiler options.
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-fsanitize=address")
#set(ASAN_OPTIONS=check_initialization_order=1)
#set(ASAN_OPTIONS=detect_leaks=1)

set(CMAKE_CXX_FLAGS "-fexceptions")

# Needed to run Boost <v1.64 on c++17, (consider upgrading boost?)
set(CMAKE_CXX_FLAGS "-fno-new-ttp-matching")

add_definitions(-ggdb3 -O0 -Wall)       # Debugging + No optimization
#add_definitions(-O3 -Wall)              # Much optimization

# Definitions needed to eliminate runtime dependency
# on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

#########################################
# ------------------------------------- #
# -------- SOURCES AND HEADERS -------- #
# ------------------------------------- #
#########################################

# Include header files
include_directories(
        "inc"
)

# Define our shared library sources. NOT test/executables.
set(SOURCES
        src/FastaReader.cpp
        src/FastaWriter.cpp
        src/SequenceElement.cpp
        src/RunlengthSequenceElement.cpp
        )

project(overlap_analysis)
add_library(overlap_analysis STATIC ${SOURCES})

# To make sure the library is named Name.so,
# get rid of the "lib" prefix.
set_target_properties(overlap_analysis PROPERTIES PREFIX "")

# Eliminate an extraneous -D during compilation.
set_target_properties(overlap_analysis PROPERTIES DEFINE_SYMBOL "")


############################################
# ---------------------------------------- #
# -------- Generating executables -------- #
# ---------------------------------------- #
############################################

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

# stupid "experimental" filesystem library for gcc <8.0
link_libraries(stdc++fs)

find_package( Boost 1.60 COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )


# -------- TESTS --------

set(TESTS

        )

foreach(FILENAME_PREFIX ${TESTS})
    add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
    set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN")
    target_link_libraries(${FILENAME_PREFIX}
            overlap_analysis
            Threads::Threads
            ${Boost_LIBRARIES}
            stdc++fs
            )

endforeach()


# -------- EXECUTABLES --------

set(EXECUTABLES
        fasta_to_RLE_fasta
        label_alignment_info
        )

foreach(FILENAME_PREFIX ${EXECUTABLES})
    add_executable(${FILENAME_PREFIX} src/executable/${FILENAME_PREFIX}.cpp)
    set_property(TARGET ${FILENAME_PREFIX} PROPERTY INSTALL_RPATH "$ORIGIN")
    target_link_libraries(${FILENAME_PREFIX}
            overlap_analysis
            Threads::Threads
            ${Boost_LIBRARIES}
            stdc++fs
            )
endforeach()

# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# The shared library goes to the bin directory.
install(TARGETS
        overlap_analysis
        DESTINATION ${BINARY_INSTALL_DIR})

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")

