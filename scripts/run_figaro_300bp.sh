#!/bin/bash
echo "================================"
echo "RUNNING FIGARO ON 300bp READS"
echo "================================"

# Create output directory for Figaro results
mkdir -p ../figaro_300bp_results

# Set the path to your 300bp reads
READ_DIR="../fully_trimmed_reads/figaro_300bp_reads"

# Run Figaro with required parameters
echo "Running Figaro on reads in $READ_DIR"
python3 ~/figaro/figaro/figaro.py \
  -i $READ_DIR \
  -o ../figaro_300bp_results \
  -n figaro_300bp_output.json \
  -a 465 \
  -f 17 \
  -r 21 \
  -m 20 \

echo "Figaro analysis complete!"
echo "Results saved to ../figaro_300bp_results/figaro_300bp_output.json"
