######################################
# ---------------------------------- #
# -------- COMPILER VERSION -------- #
# ---------------------------------- #
######################################

cmake_minimum_required(VERSION 3.10)
project(overlap_analysis VERSION 0.0.0)

message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION: ${CMAKE_CXX_COMPILER_VERSION}")

# Compiler options.
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++17")
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS )
#set(ASAN_OPTIONS=check_initialization_order=1)
#set(ASAN_OPTIONS=detect_leaks=1)

set(CMAKE_CXX_FLAGS "-fno-new-ttp-matching -fsanitize=address -fexceptions")

add_definitions(-ggdb3 -O0 -Wall)       # Debugging + No optimization
#add_definitions(-O3 -Wall)              # Much optimization

# Definitions needed to eliminate runtime dependency
# on the boost system library.
add_definitions(-DBOOST_SYSTEM_NO_DEPRECATED)
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

# Some potentially helpful messages
message(STATUS "CMAKE_BINARY_DIR: ${CMAKE_BINARY_DIR}")
message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")
message(STATUS "CMAKE_CURRENT_BINARY_DIR: ${CMAKE_CURRENT_BINARY_DIR}")
message(STATUS "CMAKE_CURRENT_SOURCE_DIR: ${CMAKE_CURRENT_SOURCE_DIR}")


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
        src/AlignmentChain.cpp
        src/OverlapMap.cpp
        src/FastaReader.cpp
        src/FastaWriter.cpp
        src/FastqIterator.cpp
        src/graph_utils.cpp
        src/PafElement.cpp
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

##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ------------------ ogdf ------------------ #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

# Download or update library as an external project
ExternalProject_Add(project_ogdf
        GIT_REPOSITORY https://github.com/ogdf/ogdf.git
        PREFIX ${CMAKE_SOURCE_DIR}/external/ogdf
#        DOWNLOAD_COMMAND ""
#        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/ogdf/ -DCMAKE_INSTALL_LIBDIR=${CMAKE_SOURCE_DIR}/external/ogdf/lib
        STAMP_DIR ${CMAKE_SOURCE_DIR}/external/ogdf/stamp
#        SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/ogdf/src
        BINARY_DIR ${CMAKE_SOURCE_DIR}/external/ogdf/build
        INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/ogdf/
        INSTALL_COMMAND make install
        )

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_ogdf INSTALL_DIR)

# Create new library for external project (so it can be linked with main library)
add_library(ogdf STATIC IMPORTED)
set_property(TARGET ogdf
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/x86_64-linux-gnu/libOGDF.a)

# Create new library for external project (so it can be linked with main library)
add_library(coin STATIC IMPORTED)
set_property(TARGET coin
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/x86_64-linux-gnu/libCOIN.a)

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(ogdf
        project_ogdf
        coin)
add_dependencies(overlap_analysis ogdf)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(external/ogdf/include/)

# Add runtime path for main library so it can access external library
#set_property(TARGET GetBlunted PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/ogdf")

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")

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

#set(Boost_DEBUG 1)
find_package( Boost 1.60 COMPONENTS program_options REQUIRED )
include_directories( ${Boost_INCLUDE_DIR} )


# -------- TESTS --------

set(TESTS
        test_boost_interval_map
        test_chaining
        test_ogdf
#        test_networkit
        )

foreach(FILENAME_PREFIX ${TESTS})
    add_executable(${FILENAME_PREFIX} src/test/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            overlap_analysis
            Threads::Threads
            ${Boost_LIBRARIES}
            stdc++fs
            ogdf
            coin
#            networkit
#            tlx
            )

endforeach()


# -------- EXECUTABLES --------

set(EXECUTABLES
        fasta_to_RLE_fasta
        filter_chimeras_from_paf
        filter_palindromes_by_quality
        find_overlapping_alignments
        label_alignment_info
        plot_alignment_as_graph
        plot_alignment_as_double_stranded_graph
        )

foreach(FILENAME_PREFIX ${EXECUTABLES})
    add_executable(${FILENAME_PREFIX} src/executable/${FILENAME_PREFIX}.cpp)
    target_link_libraries(${FILENAME_PREFIX}
            overlap_analysis
            Threads::Threads
            ${Boost_LIBRARIES}
            stdc++fs
            ogdf
            coin
#            networkit
#            tlx
            )
endforeach()

# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# The shared library goes to the bin directory.
install(TARGETS
        overlap_analysis
        DESTINATION ${BINARY_INSTALL_DIR})


