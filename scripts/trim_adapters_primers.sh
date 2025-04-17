#!/bin/bash

# Create output directory
mkdir -p ../fully_trimmed_reads

# Define your adapter and primer sequences
FWD_ADPT="TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG"
REV_ADPT="GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAG"
FWD_PRIMER="CCTACGGGNGGCWGCAG"
REV_PRIMER="GACTACHVGGGTATCTAATCC"

# Create the complete sequences (adapter+primer)
FWD_SEQ="${FWD_ADPT}${FWD_PRIMER}"
REV_SEQ="${REV_ADPT}${REV_PRIMER}"

echo "================================"
echo "TRIMMING ADAPTERS AND PRIMERS"
echo "--------------------------------"
echo "Forward adapter: $FWD_ADPT"
echo "Reverse adapter: $REV_ADPT"
echo "Forward primer: $FWD_PRIMER"
echo "Reverse primer: $REV_PRIMER"
echo "--------------------------------"
echo "Combined forward sequence: $FWD_SEQ"
echo "Combined reverse sequence: $REV_SEQ"
echo "================================"

# Process all sample pairs in numerical order
cd ../Samples
for num in $(ls Sample*_R1.fastq.gz | sed 's/Sample//' | sed 's/_R1.fastq.gz//' | sort -n); do
  r1_file="Sample${num}_R1.fastq.gz"
  r2_file="Sample${num}_R2.fastq.gz"
  
  # Skip if files don't exist
  if [ ! -f "$r1_file" ] || [ ! -f "$r2_file" ]; then
    continue
  fi
  
  echo "Processing Sample${num}..."
  
  # Run trimming with both adapters and primers
  cutadapt \
    -a $FWD_SEQ \
    -A $REV_SEQ \
    -o ../fully_trimmed_reads/trimmed-Sample${num}_R1.fastq.gz \
    -p ../fully_trimmed_reads/trimmed_Sample${num}_R2.fastq.gz \
    $r1_file $r2_file
    
  echo "Completed processing: Sample${num}"
  echo "----------------------------------------"
done

echo "All processing complete!"
echo "Fully trimmed files are in the '../fully_trimmed_reads' directory"
