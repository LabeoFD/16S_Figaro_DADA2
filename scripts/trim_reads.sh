#!/bin/bash

# Go to root directory
cd .. 

# Create output directory
mkdir -p trimmed_reads

# Define parameters
ADAPTER_R1="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
ADAPTER_R2="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
TARGET_LENGTH=280

# Print parameters being used
echo "================================"
echo "ADAPTER TRIMMING PARAMETERS"
echo "--------------------------------"
echo "R1 Adapter: $ADAPTER_R1"
echo "R2 Adapter: $ADAPTER_R2"
echo "Target Length for 301bp reads: $TARGET_LENGTH bp"
echo "================================"

cd ./Samples

# Process all paired-end read files
for r1_file in *_R1.fastq.gz; do
  # Get the corresponding R2 file
  r2_file=${r1_file/_R1/_R2}

  # Check if R2 file exists
  if [ -f "$r2_file" ]; then
    echo "Processing pair: $r1_file and $r2_file"

    # First, trim adapters without length filtering
    cutadapt \
      -a $ADAPTER_R1 \
      -A $ADAPTER_R2 \
      -o ../trimmed_reads/temp_"$r1_file" \
      -p ../trimmed_reads/temp_"$r2_file" \
      "$r1_file" "$r2_file"

    # Then, trim 301bp reads to 280bp and keep only 280bp reads
    cutadapt \
      --length 301 \
      -u -21 -U -21 \
      -o ../trimmed_reads/temp2_"$r1_file" \
      -p ../trimmed_reads/temp2_"$r2_file" \
      ../trimmed_reads/temp_"$r1_file" \
      ../trimmed_reads/temp_"$r2_file"

    # Now filter to keep only exact 280bp reads
    cutadapt \
      --minimum-length $TARGET_LENGTH \
      --maximum-length $TARGET_LENGTH \
      -o ../trimmed_reads/"${r1_file%.fastq.gz}_280bp.fastq.gz" \
      -p ../trimmed_reads/"${r2_file%.fastq.gz}_280bp.fastq.gz" \
      ../trimmed_reads/temp2_"$r1_file" \
      ../trimmed_reads/temp2_"$r2_file"
    
    # Remove temporary files
    rm ../trimmed_reads/temp_"$r1_file" ../trimmed_reads/temp_"$r2_file"
    rm ../trimmed_reads/temp2_"$r1_file" ../trimmed_reads/temp2_"$r2_file"

    echo "Completed processing: $r1_file and $r2_file"
    echo "----------------------------------------"
  else
    echo "WARNING: No matching R2 file found for $r1_file"
    echo "----------------------------------------"
  fi
done

echo "All processing complete!"
echo "Trimmed files are in the '../trimmed_reads' directory"

