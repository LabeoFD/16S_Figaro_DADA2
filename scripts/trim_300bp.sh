#!/bin/bash
echo "================================"
echo "TRIMMING ALL READS TO 300bp FOR FIGARO"
echo "================================"
# Create output directory
mkdir -p ../fully_trimmed_reads/figaro_300bp_reads

# Process all trimmed file pairs
for r1_file in ../fully_trimmed_reads/trimmed-Sample*_R1.fastq.gz; do
  # Get corresponding R2 file
  r2_file=${r1_file/R1/R2}
  
  # Get just the sample number from the filename
  sample_num=$(echo $(basename $r1_file) | grep -o 'Sample[0-9]\+')
  
  # Create output filenames with the correct format for Figaro
  output_r1="${sample_num}_trimmed_300bp_R1.fastq.gz"
  output_r2="${sample_num}_trimmed_300bp_R2.fastq.gz"
  
  echo "Trimming to 300bp for Figaro: $(basename $r1_file) and $(basename $r2_file)"
  cutadapt \
    --length 300 \
    --minimum-length 300 \
    -o ../fully_trimmed_reads/figaro_300bp_reads/$output_r1 \
    -p ../fully_trimmed_reads/figaro_300bp_reads/$output_r2 \
    $r1_file $r2_file
    
  echo "Completed trimming: $output_r1 and $output_r2"
  echo "----------------------------------------"
done

echo "All trimming complete!"
echo "300bp fixed-length files for Figaro are in the '../fully_trimmed_reads/figaro_300bp_reads' directory"
