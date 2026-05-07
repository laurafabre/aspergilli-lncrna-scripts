library(GENESPACE)

wd <- "data_analysis/genespace"
path2mcscanx = "softwares/MCScanX"
genomeRepo <- "data_analysis/genespace"

# Parse
# list of folders of the species: copy and paste the protein fasta and gff3 downloaded from FungiDB (ref strains)
list_ref <- c("afla", "afum", "anid", "anig")

parsedPaths3 <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = c(list_ref),
  gffString = "gff3",
  gffIdColumn='ID',
  gffStripText=';',
  faString = "fa",
  headerSep = " | ",
  headerStripText = "gene=",
  headerEntryIndex = 5,
  presets = "none",
  genespaceWd = wd,
  troubleShoot = TRUE,
  overwrite = TRUE)

# checking
gpar <- init_genespace(wd = wd,
  path2mcscanx = path2mcscanx)
  
# run 
out <- run_genespace(gpar, overwrite = T)

#out <- run_genespace(gpar, overwrite = T, overwriteInBlkOF = T)  #'run_orthofinderInBlk' optionally re-runs orthofinder within each syntenic block, returning phylogenetically hierarchical orthogroups (HOGs)
