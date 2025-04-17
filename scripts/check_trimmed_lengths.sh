#!/bin/bash

# Check length distribution of fully trimmed reads
echo "ANALYZING READ LENGTHS AFTER ADAPTER AND PRIMER TRIMMING"
echo "=================================================="

# For R1 files
echo "ANALYZING R1 FILES..."
for num in $(ls ../fully_trimmed_reads/trimmed-Sample*_R1.fastq.gz | sed 's/.*trimmed-Sample//' | sed 's/_R1.fastq.gz//' | sort -n); do
  file="../fully_trimmed_reads/trimmed-Sample${num}_R1.fastq.gz"
  
  # Skip if file doesn't exist
  if [ ! -f "$file" ]; then
    continue
  fi
  
  echo "Processing Sample${num} R1..."
  zcat $file | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | sort -nr | head -10
  echo "-----------------"
done

# For R2 files
echo "ANALYZING R2 FILES..."
for num in $(ls ../fully_trimmed_reads/trimmed_Sample*_R2.fastq.gz | sed 's/.*trimmed_Sample//' | sed 's/_R2.fastq.gz//' | sort -n); do
  file="../fully_trimmed_reads/trimmed_Sample${num}_R2.fastq.gz"
  
  # Skip if file doesn't exist
  if [ ! -f "$file" ]; then
    continue
  fi
  
  echo "Processing Sample${num} R2..."
  zcat $file | awk 'NR%4==2 {print length($0)}' | sort -n | uniq -c | sort -nr | head -10
  echo "-----------------"
done
