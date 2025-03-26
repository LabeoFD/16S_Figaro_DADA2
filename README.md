# 16S_Figaro_DADA2

# Figaro 
FIGARO is a program written by the folks at Zymo Research to take the guess work out of deciding what truncation parameters to use with the QIIME2 DADA2 plug-in. These parameters are chosen to minimize expected errors in forward and reverse reads, provide enough overlap length to merge the reads, and simultaneously maximize the number of reads passing this filtration step. A pre-print describing the program and its use is available [here](https://www.biorxiv.org/content/10.1101/610394v1). The program is available on [GitHub](https://github.com/Zymo-Research/figaro) along with instructions for installing a Docker version and a python command line version.  

# Installation



# Requirements
FIGARO requires that all sequences it is given to scan are the same length. This should be true if the sequences are directly from the sequencing facility. If, however, you have altered them in some way they may not be and FIGARO will fail. For example, if you have already trimmed the sequences of primers, they may differ in length, especially between forward and reverse reads. You can remedy that when trimming by specifying that the minimum and maximum length be the same, or after trimming by filtering all of the sequences to be the same length, i.e. the minimum length found among all of the sequences.

# Step1: Check read lengths
- Script used: `check_read_lengths.sh`
**Analysis of Read Length Distribution for Sample1**
Looking at your read length distribution for Sample1, we can observe the following:
**Forward Reads (R1)**
The forward reads show these key patterns:
Majority length: 301bp (145,545 reads)
Secondary peaks: 241bp (40,207), 240bp (21,172), and 239bp (16,207)
The steep drop-off from 301bp to other lengths suggests sequencing was done to achieve a target read length of 301bp
**Reverse Reads (R2)**
The reverse reads show a slightly different pattern:
Majority length: 300bp (143,371 reads)
Notable difference: The most common length is 300bp rather than 301bp
Secondary peaks: Similar to R1 with 241bp (37,212), 301bp (29,893), and 240bp (19,033)

NOTE: The fact that R1 reads are predominantly 301bp while R2 reads are predominantly 300bp could indicate slightly different quality profiles between forward and reverse reads. This is typical in Illumina paired-end sequencing, where reverse reads often have slightly lower quality.


# Blogs and References

1. [qiime2 Forum] https://forum.qiime2.org/t/automatically-predicting-trimming-parameters/30753/13
2. [Figaro and DADA2] (https://github.com/Zymo-Research/figaro/issues/40)
3. [Figaro tutorial](https://john-quensen.com/tutorials/figaro/)
