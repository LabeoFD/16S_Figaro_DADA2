#!/bin/bash

echo "================================"
echo "FILTERING READS TO 280bp ONLY"
echo "================================"

# Create output directory for 280bp reads
mkdir -p ../fully_trimmed_reads/figaro_280bp_reads

# Process all trimmed file pairs
for r1_file in ../fully_trimmed_reads/trimmed-Sample*_R1.fastq.gz; do
  # Get corresponding R2 file
  r2_file=${r1_file/R1/R2}

  # Get just the filename without the path
  filename_r1=$(basename $r1_file)
  filename_r2=$(basename $r2_file)

  # Get sample number
  sample_num=$(echo $filename_r1 | sed 's/trimmed-Sample\([0-9]*\)_R1.fastq.gz/\1/')
  
  # Create output filenames that maintain R1/R2 pattern
  output_r1="280bp-Sample${sample_num}_R1.fastq.gz"
  output_r2="280bp-Sample${sample_num}_R2.fastq.gz"

  echo "Filtering to 280bp only: $filename_r1 and $filename_r2"

  cutadapt \
    --minimum-length 280 \
    --maximum-length 280 \
    -o ../fully_trimmed_reads/figaro_280bp_reads/$output_r1 \
    -p ../fully_trimmed_reads/figaro_280bp_reads/$output_r2 \
    $r1_file $r2_file

  echo "Completed filtering: $filename_r1 and $filename_r2"
  echo "----------------------------------------"
done

echo "All filtering complete!"
echo "280bp-only files are in the '../fully_trimmed_reads/figaro_280bp_reads' directory"
