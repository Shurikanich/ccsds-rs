cmake_minimum_required(VERSION 3.10)
project(ccsds_toolkit)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define FEC and cc_soft paths
set(FEC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/fec-3.0.1)
set(CC_SOFT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/cc_soft)

# Add subdirectory for cc_soft library
add_subdirectory(${CC_SOFT_DIR})

# Include paths
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FEC_DIR}
    ${CC_SOFT_DIR}  # Add cc_soft headers
)

# FEC source files (C-only)
set(FEC_SOURCES
    ${FEC_DIR}/init_rs_char.c
    ${FEC_DIR}/encode_rs_8.c
    ${FEC_DIR}/decode_rs_8.c
    ${FEC_DIR}/encode_rs_ccsds.c
    ${FEC_DIR}/decode_rs_ccsds.c
)

# Your CCSDS C++ source files
set(CCSDS_SOURCES
    main.cc
    ccsds_rs_encoder.cc
    ccsds_rs_decoder.cc
    correlator.cc
    reed_solomon.cc
    rs_tables.cc
)

# Build FEC library
add_library(fec STATIC ${FEC_SOURCES})

# Final executable
add_executable(ccsds_main ${CCSDS_SOURCES})

# Link your code with fec and cc_soft libraries
target_link_libraries(ccsds_main fec cc_soft)