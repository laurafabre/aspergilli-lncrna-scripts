# 2_predict_cis_targets.R
library(GenomicRanges)
library(dplyr)

# 2 ways to define nearby genes, or distance of 10 kb or 3 genes up and down of the lncrna

# Function to read BED file and create a data frame of gene coordinates
read_bed_file <- function(bed_file) {
  bed_df <- read.table(bed_file, header = FALSE, stringsAsFactors = FALSE)
  colnames(bed_df) <- c("chrom", "start", "end", "gene_id", "strand")
  return(bed_df)
}

# Paths to the BED files
bed_file_afla <- "../orthology_synteny/bed_files/afla_lncrna_intergenic.bed"
bed_file_afum <- "../orthology_synteny/bed_files/afum_lncrna_intergenic.bed"
bed_file_anid <- "../orthology_synteny/bed_files/anid_lncrna_intergenic.bed"
bed_file_anig <- "../orthology_synteny/bed_files/anig_lncrna_intergenic.bed"
#bed_phi_afum <- "../paper/PHI_afum_nearbygenes_correlated.bed"

# Read the BED files
afla_bed <- read_bed_file(bed_file_afla)
afum_bed <- read_bed_file(bed_file_afum)
anid_bed <- read_bed_file(bed_file_anid)
anig_bed <- read_bed_file(bed_file_anig)
#phi_afum <- read_bed_file(bed_phi_afum)

# Add species information for genes
afla_bed$Species <- "A.flavus"
afum_bed$Species <- "A.fumigatus"
anid_bed$Species <- "A.nidulans"
anig_bed$Species <- "A.niger"

# Paths to the gene BED files
bed_file_genes_afla <- "bed_files/afla_mRNA.bed"
bed_file_genes_afum <- "bed_files/afum_mRNA.bed"
bed_file_genes_anid <- "bed_files/anid_mRNA.bed"
bed_file_genes_anig <- "bed_files/anig_mRNA.bed"

# Read the gene BED files
afla_genes_bed <- read_bed_file(bed_file_genes_afla)
afum_genes_bed <- read_bed_file(bed_file_genes_afum)
anid_genes_bed <- read_bed_file(bed_file_genes_anid)
anig_genes_bed <- read_bed_file(bed_file_genes_anig)

# Add species information for genes
afla_genes_bed$species <- "afla"
afum_genes_bed$species <- "afum"
anid_genes_bed$species <- "anid"
anig_genes_bed$species <- "anig"

# genomic range
find_nearby_genes_for_species <- function(lncrna_bed, gene_bed, distance_threshold = 100000) {
  lncrna_bed$strand <- factor(lncrna_bed$strand, levels = c("+", "-", "*"))
  gene_bed$strand <- factor(gene_bed$strand, levels = c("+", "-", "*"))
  
  # Create GRanges objects
  lncrna_gr <- GRanges(
    seqnames = lncrna_bed$chrom,
    ranges = IRanges(start = lncrna_bed$start, end = lncrna_bed$end),
    strand = lncrna_bed$strand
  )
  
  genes_gr <- GRanges(
    seqnames = gene_bed$chrom,
    ranges = IRanges(start = gene_bed$start, end = gene_bed$end),
    strand = gene_bed$strand
  )
  
  # Find nearby genes and calculate distances
  nearby_genes <- findOverlaps(lncrna_gr, genes_gr, maxgap = distance_threshold, ignore.strand = TRUE)
  distances <- distance(lncrna_gr[queryHits(nearby_genes)], genes_gr[subjectHits(nearby_genes)], ignore.strand = TRUE)
  
  # Create data frame with distances and add strand information
  nearby_genes_df <- data.frame(
    lncRNA_ID = lncrna_bed$gene_id[queryHits(nearby_genes)],
    gene_ID = gene_bed$gene_id[subjectHits(nearby_genes)],
    lncRNA_strand = as.character(lncrna_bed$strand[queryHits(nearby_genes)]),
    gene_strand = as.character(gene_bed$strand[subjectHits(nearby_genes)]),
    species = gene_bed$species[queryHits(nearby_genes)],
    distance = distances
  )
  
  return(nearby_genes_df)
}


find_nearby_genes_10kb_each_side <- function(lncrna_bed, gene_bed) {
  # Create GRanges objects
  lncrna_gr <- GRanges(
    seqnames = lncrna_bed$chrom,
    ranges = IRanges(start = lncrna_bed$start, end = lncrna_bed$end),
    strand = lncrna_bed$strand,
    gene_id = lncrna_bed$gene_id
  )
  
  genes_gr <- GRanges(
    seqnames = gene_bed$chrom,
    ranges = IRanges(start = gene_bed$start, end = gene_bed$end),
    strand = gene_bed$strand,
    gene_id = gene_bed$gene_id,
    species = gene_bed$species
  )
  
  # Expand lncRNA regions by 10kb on each side
  lncrna_expanded <- promoters(lncrna_gr, upstream = 100000, downstream = 100000 + width(lncrna_gr))
  
  # Find overlaps (genes within 10kb upstream or downstream)
  overlaps <- findOverlaps(lncrna_expanded, genes_gr, ignore.strand = TRUE)
  
  # Calculate precise distances (negative for upstream, positive for downstream)
  distance_df <- data.frame(
    lncRNA_ID = lncrna_gr$gene_id[queryHits(overlaps)],
    gene_ID = genes_gr$gene_id[subjectHits(overlaps)],
    lncRNA_strand = as.character(strand(lncrna_gr[queryHits(overlaps)])),
    gene_strand = as.character(strand(genes_gr[subjectHits(overlaps)])),
    species = genes_gr$species[subjectHits(overlaps)],
    distance = distance(lncrna_gr[queryHits(overlaps)], genes_gr[subjectHits(overlaps)], ignore.strand = TRUE)
  )
  
  # Adjust distance sign based on relative position and strand
  distance_df$signed_distance <- apply(distance_df, 1, function(row) {
    lnc_pos <- ifelse(row["lncRNA_strand"] == "-", 
                      end(lncrna_gr[lncrna_gr$gene_id == row["lncRNA_ID"]]),
                      start(lncrna_gr[lncrna_gr$gene_id == row["lncRNA_ID"]]))
    gene_pos <- ifelse(row["gene_strand"] == "-",
                       end(genes_gr[genes_gr$gene_id == row["gene_ID"]]),
                       start(genes_gr[genes_gr$gene_id == row["gene_ID"]]))
    
    if (row["lncRNA_strand"] == "-") {
      gene_pos - lnc_pos  # For negative strand, upstream is positive direction
    } else {
      lnc_pos - gene_pos  # For positive strand, upstream is negative direction
    }
  })
  
  # Classify as upstream or downstream
  distance_df$position <- ifelse(distance_df$signed_distance > 0, "upstream", "downstream")
  
  return(distance_df)
}
nearby_genes_afum2 <- find_nearby_genes_10kb_each_side(afum_bed, afum_genes_bed)
nearby_genes_afum <- find_nearby_genes_for_species(afum_bed, afum_genes_bed)
nearby_genes_afla <- find_nearby_genes_for_species(afla_bed, afla_genes_bed)
nearby_genes_anid <- find_nearby_genes_for_species(anid_bed, anid_genes_bed)
nearby_genes_anig <- find_nearby_genes_for_species(anig_bed, anig_genes_bed)

# Combine all the gene data into one data frame
cis_all <- rbind(nearby_genes_afum, nearby_genes_afla, nearby_genes_anid, nearby_genes_anig)

# Define the distance threshold for filtering (for example, 1000 bp)
distance_threshold_utr <- 100  # Adjust the distance threshold as needed

# Filter out genes in a putative UTR region only if they are on the same strand
filtered_nearby_genes <- cis_all %>%
  filter(
    # Retain genes on the opposite strand regardless of distance
    gene_strand != lncRNA_strand |
      # Keep genes on the same strand only if they are outside the UTR region threshold
      distance > distance_threshold_utr
  )

# Identify lncRNAs that appear at <500 bp at least once
lncRNAs_to_remove <- cis_all %>%
  filter(gene_strand == lncRNA_strand & distance <= 100) %>%
  pull(lncRNA_ID) %>%
  unique() 

filtered_cis_all <-cis_all %>%
  filter(!lncRNA_ID %in% lncRNAs_to_remove)

save(lncRNAs_to_remove, filtered_nearby_genes, file = "../Aspergilli_paper/data/processed/correlation_function/lncrna_to_remove_100bpdistance_UTR_and_nearbygenes.RData")
write.csv(cis_all, "../Aspergilli_paper/data/processed/correlation_function/cis_targets_all_species.csv", row.names = FALSE)
write.csv(filtered_cis_all, "../Aspergilli_paper/data/processed/correlation_function/cis_targets_all_species_filtered_100bp.csv", row.names = FALSE)

#Filter if Needed
# filter for specific strand relationships (e.g., same strand), specific genes, distances <10kb, etc.
cis_all_10kb <- subset(filtered_cis_all, distance <= 10000)