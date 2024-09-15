#!/bin/bash

# Function to check if the index files exist
function check_index {
    for ext in amb ann bwt pac sa; do
        if [[ ! -f "${1}.${ext}" ]]; then
            return 1  # Return 1 if any index file is missing
        fi
    done
    return 0  # Return 0 if all index files exist
}

# Run index only if it's missing
if ! check_index $1; then
    echo "Index files not found. Running bwa index..."
    bwa index $1
else
    echo "Index files already exist. No need to re-index."
fi