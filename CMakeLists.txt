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

set(CMAKE_CXX_FLAGS "-fsanitize=address -fexceptions")
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
        src/PafAlignmentChain.cpp
        src/Color.cpp
        src/DagAligner.cpp
        src/DiagonalTree.cpp
        src/DiagonalIntervalTree.cpp
        src/dotplot.cpp
        src/FastaReader.cpp
        src/FastaWriter.cpp
        src/FastqIterator.cpp
        src/Graph.cpp
        src/graph_utils.cpp
        src/Kernel.cpp
        src/OverlapMap.cpp
        src/PafElement.cpp
        src/Plot.cpp
        src/SvgPlot.cpp
        src/SequenceElement.cpp
        src/SweepLineSolver.cpp
        src/wavefront.cpp
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


##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ---------------- mummer4 ----------------- #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

# Download or update library as an external project
ExternalProject_Add(project_mummer
        URL https://github.com/mummer4/mummer/releases/download/v4.0.0rc1/mummer-4.0.0rc1.tar.gz
#        PREFIX ${CMAKE_SOURCE_DIR}/external/mummer
#        DOWNLOAD_COMMAND ""
#        UPDATE_COMMAND ""
#        STAMP_DIR ${CMAKE_SOURCE_DIR}/external/mummer/stamp
        BUILD_IN_SOURCE TRUE
        SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/mummer/
        CONFIGURE_COMMAND ${CMAKE_SOURCE_DIR}/external/mummer/configure --prefix=${CMAKE_SOURCE_DIR}/external/mummer/
        INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/mummer/
        BUILD_COMMAND make
        INSTALL_COMMAND make install
        )

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_mummer INSTALL_DIR)

# Create new library for external project (so it can be linked with main library)
add_library(mummer4 SHARED IMPORTED)
set_property(TARGET mummer4
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libumdmummer.so)

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(mummer4
        project_mummer)

add_dependencies(overlap_analysis mummer4)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(external/mummer/include/)

# Add runtime path for main library so it can access external library
#set_property(TARGET GetBlunted PROPERTY INSTALL_RPATH "$ORIGIN" "${INSTALL_DIR}/bin/ogdf")

message(STATUS "INSTALL_DIR: ${INSTALL_DIR}")

##############################################
# ------------------------------------------ #
# -------- LINKING EXTERNAL LIBRARY -------- #
# ----------------- edlib ------------------ #
# ------------------------------------------ #
##############################################

# Need to explicitly enable ExternalProject functionality
include(ExternalProject)

# Download or update library as an external project
ExternalProject_Add(project_edlib
        GIT_REPOSITORY https://github.com/Martinsos/edlib.git
        PREFIX ${CMAKE_SOURCE_DIR}/external/edlib
        #        DOWNLOAD_COMMAND ""
        #        UPDATE_COMMAND ""
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_SOURCE_DIR}/external/edlib/ -DCMAKE_INSTALL_LIBDIR=${CMAKE_SOURCE_DIR}/external/edlib/lib
        STAMP_DIR ${CMAKE_SOURCE_DIR}/external/edlib/stamp
        #        SOURCE_DIR ${CMAKE_SOURCE_DIR}/external/edlib/src
        BINARY_DIR ${CMAKE_SOURCE_DIR}/external/edlib/build
        INSTALL_DIR ${CMAKE_SOURCE_DIR}/external/edlib/
        INSTALL_COMMAND make install
        )

# Define INSTALL_DIR as the install directory for external library
ExternalProject_Get_Property(project_edlib INSTALL_DIR)

# Create new library for external project (so it can be linked with main library)
add_library(edlib STATIC IMPORTED)
set_property(TARGET edlib
        PROPERTY IMPORTED_LOCATION ${INSTALL_DIR}/lib/libedlib.a)

# Define main library as dependent on the downloaded project (transitively)
add_dependencies(edlib
        project_edlib)

add_dependencies(overlap_analysis edlib)

# Ensure that main library has access to primary dependencies' and secondary dependencies' headers
include_directories(external/edlib/include/)

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
        test_boost_rtree
        test_chaining
        test_dag
        test_diagonal_interval_tree
        test_diagonal_tree
        test_edlib_runtime
        test_edlib
        test_ogdf
        test_gnuplot
        test_kernel
        test_mummer
        test_mummer_runtime
        test_svgplot
        test_sweepline
        test_wf
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
            mummer4
            edlib
            )

endforeach()


# -------- EXECUTABLES --------

set(EXECUTABLES
        chain_alignment
        dotplot
        evaluate_overlaps
        filter_chimeras_from_paf
        filter_palindromes_by_quality
        find_overlapping_alignments
        label_alignment_files
        plot_alignment_as_double_stranded_graph
        self_align
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
            mummer4
            edlib
            )

endforeach()

# -------- final steps --------

# Where to install
set(BINARY_INSTALL_DIR ${CMAKE_CURRENT_SOURCE_DIR}/bin)

# The shared library goes to the bin directory.
install(TARGETS
        overlap_analysis
        DESTINATION ${BINARY_INSTALL_DIR})


