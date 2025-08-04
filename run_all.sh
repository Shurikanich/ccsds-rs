#!/bin/bash
# run_all.sh
# This script runs the simulation executable (ccsds_main) with all configuration files in ./conf

CONFIG_DIR="./conf"
SIM_BIN="./build/ccsds_main"
OUTPUT_DIR="./res"

# If "clean" argument is passed, clear the output directory and exit
if [ "$1" == "clean" ]; then
    echo "Cleaning output directory: $OUTPUT_DIR"
    rm -f "$OUTPUT_DIR"/*
    echo "Output directory cleaned."
    exit 0
fi

# Gather all .txt configuration files in CONFIG_DIR
CONFIG_FILES=("$CONFIG_DIR"/*.txt)

# Check if any configuration files exist
if [ ${#CONFIG_FILES[@]} -eq 0 ]; then
    echo "No configuration files found in $CONFIG_DIR."
    exit 1
fi

# Run simulations for each configuration file
for cfg in "${CONFIG_FILES[@]}"; do
    if [ -f "$cfg" ]; then
        echo "Running simulation with $cfg..."
        "$SIM_BIN" "$cfg"
    else
        echo "Warning: Configuration file $cfg not found."
    fi
done

echo "All simulations complete."
