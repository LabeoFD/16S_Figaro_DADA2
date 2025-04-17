#!/bin/bash
echo "================================"
echo "FILTERING READS TO 280bp ONLY"
echo "================================"
# Create output directory for 280bp reads
mkdir -p ../trimmed_reads/figaro_280bp_reads
# Process all trimmed file pairs
for r1_file in ../trimmed_reads/trimmed-Sample*_R1.fastq.gz; do
  # Get corresponding R2 file
  r2_file=${r1_file/R1/R2}
  # Get just the filename without the path
  filename_r1=$(basename $r1_file)
  filename_r2=$(basename $r2_file)
  # Create output filenames
  output_r1=${$filename_r1/_trimmed/_280bp/}
  output_r2=${$filename_r2/_trimmed/_280bp/}
  echo "Filtering to 280bp only: $filename_r1 and $filename_r2"
  cutadapt \
    --minimum-length 280 \
    --maximum-length 280 \
    -o ../trimmed_reads/figaro_280bp_reads/$output_r1 \
    -p ../trimmed_reads/figaro_280bp_reads/$output_r2 \
    $r1_file $r2_file
  echo "Completed filtering: $filename_r1 and $filename_r2"
  echo "----------------------------------------"
done
echo "All filtering complete!"
echo "280bp-only files are in the '../trimmed_reads/figaro_280bp_reads' directory"
