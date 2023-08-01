cmake_minimum_required(VERSION 3.1)
project(gryphon CXX)
set(VERSION 2.0)

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

set(GRYPHON_EXTRA_SOURCES)
set(GRYPHON_EXTRA_INCLUDES)
set(GRYPHON_EXTRA_LIBRARIES)

set(CMAKE_CXX_STANDARD 14)
cmake_policy(SET CMP0048 NEW)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-std=c++14 -Wall -pedantic -ftree-vectorize -O3")
set(GRYPHON_INCLUDES)
set(GRYPHON_LIBRARIES)

# Set default build-type to release to enable performance improvements
if (NOT CMAKE_BUILD_TYPE)
    #set(CMAKE_BUILD_TYPE Release)
    set(CMAKE_BUILD_TYPE Debug)
endif(NOT CMAKE_BUILD_TYPE)

if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG")
endif (CMAKE_BUILD_TYPE STREQUAL "Debug")

message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -ffast-math")

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# PLOG (provided)
include_directories(external/plog/include)

# NamedType (provided)
include_directories(external/NamedType/include)

# C++ Threads required
#find_package(Threads REQUIRED)
#list(APPEND GRYPHON_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})

# GSL (required)
find_package(GSL REQUIRED)
list(APPEND GRYPHON_EXTRA_INCLUDES ${GSL_INCLUDE_DIR})
list(APPEND GRYPHON_EXTRA_LIBRARIES ${GSL_LIBRARIES})

# googletest (provided, see https://code.google.com/p/googletest/wiki/FAQ
#     Why is it not recommended use a pre-compiled copy of Google Test?)
option(ENABLE_TESTING "Build tests and enable test target" OFF)
if(ENABLE_TESTING)
    set(GOOGLETEST_VERSION 1.10.0)
    include_directories(external/gtest/googletest/include)
    add_subdirectory(external/gtest/googletest)
    if(APPLE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGTEST_USE_OWN_TR1_TUPLE=1")
    endif(APPLE)
endif(ENABLE_TESTING)

# Version info from Git
option(ENABLE_GIT "Embedding information about GRYPHON version from git" ON)
if(ENABLE_GIT)
    include(GetGitRevisionDescription)
    get_git_head_revision(GIT_REFSPEC GIT_SHA1)
    git_local_changes(GIT_HAS_LOCAL_CHANGES)
    set(git_revision_cpp "${CMAKE_CURRENT_BINARY_DIR}/git_revision.cpp")
    configure_file("src/utils/git_revision.cpp.in" "${git_revision_cpp}" @ONLY)
endif(ENABLE_GIT)
IF(NOT ENABLE_GIT OR (${GIT_SHA1} STREQUAL "GITDIR-NOTFOUND"))
    set(GIT_REFSPEC "")
    set(GIT_SHA1 "")
    set(GIT_DESC "${GRYPHON_RELEASE_VERSION}-no-git")
endif()
message(STATUS "GRYOPHON version: ${GIT_DESC} ${GIT_SHA1} ${GIT_REFSPEC}")


# ----------------------------------------------------------------------------
# Install
# ----------------------------------------------------------------------------
include_directories(include ${GRYPHON_EXTRA_INCLUDES})

add_library(gryphon SHARED
    # src/analyticalSolutions/beniamino.cpp
    # src/analyticalSolutions/cosmoNus.cpp
    # src/core/common.cpp
    # src/core/cosmology.cpp
    # src/core/opticalDepth.cpp
    # src/core/params.cpp
    # src/core/pid.cpp
    # src/crossSections/BreitWheeler.cpp
    # src/crossSections/KelnerAharonian2008.cpp
    # src/crossSections/PhotoPionXsecs.cpp
    # src/crossSections/PhotoDisintegrationTalysXsecs.cpp
    # src/energyLosses/AdiabaticContinuousLosses.cpp
    # src/energyLosses/BGG2006ContinuousLosses.cpp
    # src/energyLosses/PairProductionLosses.cpp
    # src/energyLosses/PhotoPionContinuousLosses.cpp
    # src/interactions/PhotoDisintegration.cpp
    # src/interactions/PhotoPionProduction.cpp
    # src/particleStacks/SingleParticleBuilder.cpp
    # src/particleStacks/SingleSourceBuilder.cpp
    # src/particleStacks/SourceEvolutionBuilder.cpp
    # src/photonFields/CmbPhotonField.cpp
    # src/photonFields/LookupTablePhotonField.cpp
    # src/photonFields/PhotonField.cpp
	src/utils/io.cpp
	src/utils/logging.cpp
    src/utils/numeric.cpp
    src/utils/progressbar.cpp
    src/utils/timer.cpp
    "${git_revision_cpp}"
    )
    target_link_libraries(gryphon ${GRYPHON_EXTRA_LIBRARIES})

# make library
#include_directories(include ${GRYPHON_INCLUDES})
#add_library(GRYPHON_LIB ${GRYPHON_SRCS} ${GRYPHON_INCLUDES})

# create the output dir
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/output)

# copy data files in build
file(GLOB GRYPHON_DATA "data/*.txt")
file(COPY ${GRYPHON_DATA} DESTINATION ${PROJECT_BINARY_DIR}/data)

# make apps
add_executable(sandbox apps/sandbox.cpp)
target_link_libraries (sandbox gryphon ${GRYPHON_EXTRA_LIBRARIES})

add_executable(run apps/run.cpp)
target_link_libraries (run gryphon ${GRYPHON_EXTRA_LIBRARIES})

# testing
# if(ENABLE_TESTING)
#     enable_testing()
#     add_executable(test_common test/testCommon.cpp)
#     target_link_libraries(test_common gryphon gtest gtest_main ${GRYPHON_EXTRA_LIBRARIES})
#     add_test(test_common test_common)

#     add_executable(test_rng test/testRng.cpp)
#     target_link_libraries(test_rng gryphon gtest gtest_main ${GRYPHON_EXTRA_LIBRARIES})
#     add_test(test_rng test_rng)

#     add_executable(test_pid test/testPid.cpp)
#     target_link_libraries(test_pid gryphon gtest gtest_main ${GRYPHON_EXTRA_LIBRARIES})
#     add_test(test_pid test_pid)

#     add_executable(test_cosmology test/testCosmology.cpp)
#     target_link_libraries(test_cosmology gryphon gtest gtest_main ${GRYPHON_EXTRA_LIBRARIES})
#     add_test(test_cosmology test_cosmology)

#     add_executable(test_photonFields test/testPhotonFields.cpp)
#     target_link_libraries(test_photonFields gryphon gtest gtest_main ${GRYPHON_EXTRA_LIBRARIES})
#     add_test(test_photonFields test_photonFields)
# endif(ENABLE_TESTING)

# Show summary.
message (STATUS "CMAKE_SYSTEM .......... = ${CMAKE_SYSTEM}")
message (STATUS "BUILD_SHARED_LIBS ..... = ${BUILD_SHARED_LIBS}")
message (STATUS "CMAKE_CXX_COMPILER .... = ${CMAKE_CXX_COMPILER}")
message (STATUS "CMAKE_CXX_FLAGS ....... = ${CMAKE_CXX_FLAGS}")
message (STATUS "CMAKE_INCLUDE_PATH .... = ${CMAKE_INCLUDE_PATH}")
message (STATUS "GRYPHON_INCLUDES ...... = ${GRYPHON_EXTRA_INCLUDES}")
message (STATUS "GRYPHON_LIBRARIES ..... = ${GRYPHON_EXTRA_LIBRARIES}")
