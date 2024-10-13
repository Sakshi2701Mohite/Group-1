# SLE777 Assessment-4

This repository contains all the files related to SLE777 of the Group-1 for the Assessment-4.

## Project Guide
Dr. Ghazhanfar Khan

## Group Members
1. Sakshi Mohite
2. Babita Acharya
3. Sachin Ghuman

Part 2

Overview
This script compares coding DNA sequences from Escherichia coli and Streptococcus thermophilus. It includes counting coding sequences, calculating nucleotide and amino acid frequencies, analyzing codon usage bias, and identifying over- and under-represented k-mers in protein sequences.

Purpose
The script analyzes and contrasts genomic features between E. coli and S. thermophilus, examining differences in coding sequences, nucleotide composition, codon preferences, and protein sequence motifs. These insights help highlight evolutionary and functional distinctions.

Inputs
The script automatically downloads and decompresses coding sequence data for both organisms from Ensembl's FTP server:
•	E. coli coding sequences.
•	Streptococcus thermophilus coding sequences.

Outputs
The script generates:
1.	Coding sequence counts for each organism.
2.	Total coding DNA length for each organism.
3.	Boxplots of coding sequence length distributions.
4.	Nucleotide and amino acid frequency bar plots.
5.	Codon usage bias plot (RSCU).
6.	K-mer analysis plots for over- and under-represented 3-mers.

Code Sections
1.	Loading Libraries:
Loads seqinr, R.utils, and ggplot2 libraries for sequence manipulation, file decompression, and plotting, respectively.
2.	Data Download and Decompression:
Downloads the coding sequence files for each organism from specified URLs and decompresses them. This allows the script to read and analyze the sequences.
3.	Loading Sequences into R:
Loads the decompressed FASTA files into R as lists using read.fasta, which enables processing of the sequence data.
4.	Counting Coding Sequences:
Counts the total number of coding sequences for each organism by calculating the length of the sequence list. Results are stored in a table for comparison.
5.	Calculating Total Coding DNA Length:
Computes the total length of all coding DNA sequences for each organism by summing individual sequence lengths, then displays the results in a table.
6.	Coding Sequence Length Analysis:
Calculates mean and median lengths of the coding sequences. A boxplot is generated to visualize the length distributions for both organisms, providing insights into genome structure.
7.	Nucleotide Frequency Analysis:
Combines all sequences into a single string for each organism, counts occurrences of each nucleotide, and then plots these frequencies to show the nucleotide composition differences.
8.	Amino Acid Frequency Analysis:
Translates DNA sequences into proteins, counts occurrences of each amino acid, and plots these frequencies to compare protein composition between the organisms.
9.	Codon Usage Bias (RSCU):
Calculates Relative Synonymous Codon Usage (RSCU) for each organism to measure codon usage preferences, then plots these values to highlight biases in codon usage.
10.	K-mer Analysis in Protein Sequences:
Computes 3-mer frequencies in protein sequences and identifies the top 10 over- and under-represented motifs for each organism. Plots illustrate the prevalence of specific motifs and highlight sequence differences.

Instructions for Running the Script
1.	Ensure R and required libraries (seqinr, R.utils, and ggplot2) are installed on your system.
2.	Run the script in R or RStudio. The script will automatically download the data, perform analyses, and generate visualizations.
3.	Review console outputs for tables and the plot window for graphs, which provide insights into the comparative genomic analysis.

