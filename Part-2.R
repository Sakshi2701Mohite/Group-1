# Load necessary libraries
suppressPackageStartupMessages({
  library(seqinr)
  library(R.utils)
  library(ggplot2)
})

# Define URLs for sequence data
ecoli_url <- "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_0_collection/escherichia_coli_str_k_12_substr_mg1655_gca_000005845/cds/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.cds.all.fa.gz"
streptococcus_url <- "https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/release-59/fasta/bacteria_121_collection/streptococcus_thermophilus_cag_236_gca_000434755/cds/Streptococcus_thermophilus_cag_236_gca_000434755.MGS236.cds.all.fa.gz"

# Download and decompress sequence files
download.file(ecoli_url, destfile = "ecoli_cds.fa.gz")
gunzip("ecoli_cds.fa.gz")
download.file(streptococcus_url, destfile = "strepto_cds.fa.gz")
gunzip("strepto_cds.fa.gz")

# Load sequences into R
ecoli_seqs <- read.fasta("ecoli_cds.fa")
strepto_seqs <- read.fasta("strepto_cds.fa")

# 1. Number of Coding Sequences and Generating Table
coding_counts <- data.frame(
  Organism = c("E. coli", "Streptococcus thermophilus"),
  Coding_Sequences = c(length(ecoli_seqs), length(strepto_seqs))
)
cat("Coding Sequence Counts Table:\n")
print(coding_counts)
# Differences Between the Organisms:
# E. coli has significantly more coding sequences (4239) compared to Streptococcus thermophilus (1741). This suggests that E. coli may have a larger genome with more genes, or more fragmented coding regions compared to S. thermophilus.

# Calculate total coding DNA lengths
ecoli_total_length <- sum(sapply(ecoli_seqs, length))
strepto_total_length <- sum(sapply(strepto_seqs, length))

# Create a table for total coding DNA length
total_coding_dna_table <- data.frame(
  Organism = c("E. coli", "Streptococcus thermophilus"),
  Total_Coding_DNA_Length_bp = c(ecoli_total_length, strepto_total_length)
)

# Print the table
cat("Total Coding DNA Length Table:\n")
print(total_coding_dna_table)
# Differences Between the Organisms:
# E. coli has a substantially larger total coding DNA length (3,978,528 bp) compared to Streptococcus thermophilus (1,374,033 bp). This indicates that E. coli may have a larger genome and potentially a more complex set of genes than S. thermophilus.

# Calculate sequence lengths for both organisms
ecoli_lengths <- sapply(ecoli_seqs, length)
strepto_lengths <- sapply(strepto_seqs, length)

# Create a data frame for plotting
lengths_df <- data.frame(
  Organism = rep(c("E. coli", "Streptococcus thermophilus"), 
                 times = c(length(ecoli_lengths), length(strepto_lengths))),
  Sequence_Length = c(ecoli_lengths, strepto_lengths)
)
# The boxplot shows that both E. coli and Streptococcus thermophilus have similar distributions of coding sequence lengths, with the majority falling around 900-1000 bp.
# Both organisms exhibit a similar range with a few longer outliers, indicating comparable gene structures in terms of sequence length. However, S. thermophilus shows slightly more variability in longer sequences. 

# Calculate mean and median lengths
mean_ecoli <- mean(ecoli_lengths)
median_ecoli <- median(ecoli_lengths)
mean_strepto <- mean(strepto_lengths)
median_strepto <- median(strepto_lengths)

# Output mean and median lengths
cat("Mean and Median Coding Sequence Lengths:\n")
cat("E. coli - Mean:", mean_ecoli, "bp, Median:", median_ecoli, "bp\n")
cat("Streptococcus thermophilus - Mean:", mean_strepto, "bp, Median:", median_strepto, "bp\n")

# Boxplot of coding sequence lengths
boxplot(Sequence_Length ~ Organism, data = lengths_df,
        main = "Coding Sequence Length Distribution",
        xlab = "Organism", ylab = "Sequence Length (bp)",
        col = c("blue", "red"))

# Define standard nucleotides
standard_nucleotides <- c("A", "T", "C", "G")

# Combine all sequences into one long sequence for each organism
ecoli_combined_seq <- paste(unlist(ecoli_seqs), collapse = "")
strepto_combined_seq <- paste(unlist(strepto_seqs), collapse = "")

# Convert combined sequences to character vectors for counting
ecoli_nuc_seq <- s2c(ecoli_combined_seq)
strepto_nuc_seq <- s2c(strepto_combined_seq)

# Calculate nucleotide frequencies and ensure all standard nucleotides are represented
ecoli_nuc_freq <- table(factor(ecoli_nuc_seq, levels = standard_nucleotides))
strepto_nuc_freq <- table(factor(strepto_nuc_seq, levels = standard_nucleotides))

# Create a data frame for nucleotide frequencies
nucleotide_df <- data.frame(
  Nucleotide = rep(standard_nucleotides, 2),
  Frequency = c(as.numeric(ecoli_nuc_freq), as.numeric(strepto_nuc_freq)),
  Organism = rep(c("E. coli", "Streptococcus thermophilus"), each = length(standard_nucleotides))
)

# Plot nucleotide frequencies
ggplot(nucleotide_df, aes(x = Nucleotide, y = Frequency, fill = Organism)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Nucleotide Frequency in Coding Sequences", x = "Nucleotide", y = "Frequency")

# Translate sequences to protein sequences
ecoli_protein_seqs <- lapply(ecoli_seqs, translate)
strepto_protein_seqs <- lapply(strepto_seqs, translate)

# Combine all protein sequences into one for each organism and calculate amino acid frequencies
ecoli_combined_protein <- unlist(ecoli_protein_seqs)
strepto_combined_protein <- unlist(strepto_protein_seqs)

# Get unique amino acids across both organisms
all_amino_acids <- union(names(table(ecoli_combined_protein)), names(table(strepto_combined_protein)))

# Calculate amino acid frequencies and ensure all amino acids are represented
ecoli_aa_freq <- table(factor(ecoli_combined_protein, levels = all_amino_acids))
strepto_aa_freq <- table(factor(strepto_combined_protein, levels = all_amino_acids))

# Create a data frame for amino acid frequencies
aa_df <- data.frame(
  Amino_Acid = rep(all_amino_acids, 2),
  Frequency = c(as.numeric(ecoli_aa_freq), as.numeric(strepto_aa_freq)),
  Organism = rep(c("E. coli", "Streptococcus thermophilus"), each = length(all_amino_acids))
)

# Plot amino acid frequencies
ggplot(aa_df, aes(x = Amino_Acid, y = Frequency, fill = Organism)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Amino Acid Frequency in Coding Sequences", x = "Amino Acid", y = "Frequency")
# E. coli has a wider distribution of amino acid frequencies, reflecting a potentially larger and more diverse set of proteins.
# Specific amino acids, such as L (Leucine), appear significantly more in E. coli, suggesting that E. coli may favor certain amino acids more heavily than S. thermophilus.
# Define a function to calculate RSCU for all sequences in a list

calculate_rscu <- function(sequence_list) {
  codon_usage <- do.call(rbind, lapply(sequence_list, function(seq) uco(seq, index = "rscu", as.data.frame = TRUE)))
  aggregated_rscu <- aggregate(RSCU ~ codon, data = codon_usage, mean)
  return(aggregated_rscu)
}

# Calculate RSCU for each organism
ecoli_rscu <- calculate_rscu(ecoli_seqs)
strepto_rscu <- calculate_rscu(strepto_seqs)

# Add organism labels
ecoli_rscu$Organism <- "E. coli"
strepto_rscu$Organism <- "Streptococcus thermophilus"

# Combine data for plotting
combined_rscu <- rbind(ecoli_rscu, strepto_rscu)

# Plot the RSCU values to show codon usage bias
ggplot(combined_rscu, aes(x = codon, y = RSCU, fill = Organism)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Relative Synonymous Codon Usage (RSCU)", x = "Codon", y = "RSCU") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Define a custom function to calculate k-mer frequencies
calculate_kmer_frequencies <- function(sequences, k) {
  # Collapse all sequences into one single character string
  combined_sequence <- paste(unlist(sequences), collapse = "")
  
  # Initialize an empty list to store k-mers
  kmer_counts <- list()
  
  # Loop over the sequence to extract k-mers
  for (i in 1:(nchar(combined_sequence) - k + 1)) {
    kmer <- substr(combined_sequence, i, i + k - 1)
    
    if (kmer %in% names(kmer_counts)) {
      kmer_counts[[kmer]] <- kmer_counts[[kmer]] + 1
    } else {
      kmer_counts[[kmer]] <- 1
    }
  }
  
  # Convert kmer_counts list to a data frame
  kmer_df <- data.frame(
    Kmer = names(kmer_counts),
    Frequency = unlist(kmer_counts)
  )
  
  return(kmer_df)
}

# Calculate k-mer counts for k = 3 for Streptococcus thermophilus and E. coli protein sequences
strepto_kmer_3 <- calculate_kmer_frequencies(strepto_protein_seqs, 3)
ecoli_kmer_3 <- calculate_kmer_frequencies(ecoli_protein_seqs, 3)

# Define a function to identify top 10 over- and under-represented k-mers
identify_top_kmers <- function(kmer_df) {
  top_over <- head(kmer_df[order(-kmer_df$Frequency), ], 10)
  top_under <- head(kmer_df[order(kmer_df$Frequency), ], 10)
  return(list(over = top_over, under = top_under))
}

# Identify top over- and under-represented k-mers for Streptococcus thermophilus and E. coli
strepto_top_kmer_3 <- identify_top_kmers(strepto_kmer_3)
ecoli_top_kmer_3 <- identify_top_kmers(ecoli_kmer_3)

# Combine data for visualization
combined_top_kmers <- rbind(
  data.frame(Kmer = strepto_top_kmer_3$over$Kmer, Frequency = strepto_top_kmer_3$over$Frequency, Type = "Over-represented", Organism = "Streptococcus thermophilus"),
  data.frame(Kmer = strepto_top_kmer_3$under$Kmer, Frequency = strepto_top_kmer_3$under$Frequency, Type = "Under-represented", Organism = "Streptococcus thermophilus"),
  data.frame(Kmer = ecoli_top_kmer_3$over$Kmer, Frequency = ecoli_top_kmer_3$over$Frequency, Type = "Over-represented", Organism = "E. coli"),
  data.frame(Kmer = ecoli_top_kmer_3$under$Kmer, Frequency = ecoli_top_kmer_3$under$Frequency, Type = "Under-represented", Organism = "E. coli")
)

# Plot over- and under-represented k-mers
ggplot(combined_top_kmers, aes(x = Kmer, y = Frequency, fill = Organism)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~Type, scales = "free_y") +
  labs(title = "Over- and Under-represented 3-mers in Protein Sequences", x = "3-mer", y = "Frequency") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))