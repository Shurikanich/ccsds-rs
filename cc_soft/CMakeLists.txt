cmake_minimum_required(VERSION 3.11)

add_library(cc_soft STATIC
    metrics.c
    tab.c
    viterbi27.c
)

# Include headers for the library
target_include_directories(cc_soft PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})


# Optional: Compiler flags
if (CMAKE_C_COMPILER_ID STREQUAL "GNU" OR CMAKE_C_COMPILER_ID STREQUAL "Clang")
    target_compile_options(cc_soft PRIVATE -Wall -Wextra -Wpedantic)
endif()