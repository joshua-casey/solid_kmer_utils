# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)
set (CMAKE_CXX_STANDARD 17)

# Name of project and that it is C++ only.
project (solid_kmer_utils CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

find_package (ZLIB REQUIRED)

# include the sdsl lite installed in libs/sdsl-lite/install and htslib installed in libs/htslib-1.9/install
include_directories(SYSTEM libs/htslib-1.20/install/include)
link_directories(libs/htslib-1.20/install/lib)

if (NOT TARGET kmc_api)
    add_subdirectory(libs/kmc_api EXCLUDE_FROM_ALL)
endif()

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add CXX flags
set (CMAKE_CXX_FLAGS "-pthread ${CMAKE_CXX_FLAGS} -lz -O2")

# Add executable and link
add_executable(reference_kmer_counter src/reference_kmer_counter.cpp)
target_link_libraries(reference_kmer_counter hts kmc_api ${ZLIB_LIBRARIES})
