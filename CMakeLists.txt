cmake_minimum_required(VERSION 3.10)
project(ccsds_toolkit)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Define FEC source path
set(FEC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/fec-3.0.1)

# Include paths
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${FEC_DIR}
)

# FEC source files (pure C)
set(FEC_SOURCES
    ${FEC_DIR}/init_rs_char.c
    ${FEC_DIR}/encode_rs_8.c
    ${FEC_DIR}/decode_rs_8.c
    ${FEC_DIR}/encode_rs_ccsds.c
    ${FEC_DIR}/decode_rs_ccsds.c
    ${FEC_DIR}/init_rs.h
    ${FEC_DIR}/encode_rs.h
    ${FEC_DIR}/decode_rs.h
    ${FEC_DIR}/fec.h
    ${FEC_DIR}/rs-common.h
    ${FEC_DIR}/fixed.h
    ${FEC_DIR}/char.h
    ${FEC_DIR}/ccsds.h
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

# Build the FEC library
add_library(fec STATIC
    ${FEC_DIR}/init_rs_char.c
    ${FEC_DIR}/encode_rs_8.c
    ${FEC_DIR}/decode_rs_8.c
    ${FEC_DIR}/encode_rs_ccsds.c
    ${FEC_DIR}/decode_rs_ccsds.c
)

# Final executable
add_executable(ccsds_main ${CCSDS_SOURCES})

# Link your code with the fec library
target_link_libraries(ccsds_main fec)