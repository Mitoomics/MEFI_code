#!/bin/bash

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <sample_directory>"
    exit 1
fi

sample_dir=$1
insert_125="fragment_study/mis_10_splitedBy125/125_splited_sam_file"
depth_dir="depth_file"

cd "${sample_dir}/${insert_125}"

mkdir -p ../${depth_dir}

fils=$(find . -type f -name '*.sam' -size +1000c)

for s in ${fils}; do
    if samtools depth -d 500000000 -a -m 0 "${s}" > "../${depth_dir}/$(basename ${s}).txt"; then
        echo "Depth information for $(basename ${s}) has been successfully generated."
    else
        echo "Failed to generate depth information for $(basename ${s})."
    fi
done
