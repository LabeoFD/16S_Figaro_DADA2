#!/bin/bash

# Go to Samples folder
cd ../Samples

# For R1 files
echo "ANALYZING R1 FILES..."
for num in $(ls *_R1*.fastq.gz | sed 's/Sample//' | sed 's/_R1.fastq.gz//' | sort -n); do
  file="Sample${num}_R1.fastq.gz"
  
  # Skip if file doesn't exist (safety check)
  if [ ! -f "$file" ]; then
    continue
  fi
  
  echo "Processing $file..."
  zcat $file | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | sort -nr | head -10
  echo "-----------------"
done

# For R2 files
echo "ANALYZING R2 FILES..."
for num in $(ls *_R2*.fastq.gz | sed 's/Sample//' | sed 's/_R2.fastq.gz//' | sort -n); do
  file="Sample${num}_R2.fastq.gz"
  
  # Skip if file doesn't exist (safety check)
  if [ ! -f "$file" ]; then
    continue
  fi
  
  echo "Processing $file..."
  zcat $file | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | sort -nr | head -10
  echo "-----------------"
done
