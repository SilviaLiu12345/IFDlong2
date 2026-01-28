

suppressPackageStartupMessages({
  install_and_load <- function(pkgs) {
    # Set CRAN mirror if not already set
    if (is.null(getOption("repos")) || getOption("repos")["CRAN"] == "@CRAN@") {
      options(repos = c(CRAN = "https://cloud.r-project.org"))
    }
    
    for (p in pkgs) {
      if (!requireNamespace(p, quietly = TRUE)) {
        install.packages(p)
      }
      library(p, character.only = TRUE)
    }
  }
  
  install_and_load(c("data.table", "parallel"))
})

#### parameter
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 5)

mainPath <- args[1]
sampleName <- args[2]
Aligner <- args[3]
refEXON <- args[4]
ncores <- as.integer(args[5])

# Load reference CDS
cat("Loading reference exon file:", refEXON, "\n")
cds <- fread(refEXON, col.names = c("chr", "start", "end", "name", "score", "strand"))
cds[, isoform := tstrsplit(name, "__")[[1]]]

splitDir <- file.path(mainPath, Aligner, "split_example")
baseFiles <- list.files(splitDir, pattern = "example_part[0-9]+\\.bed$", full.names = TRUE)

process_pair <- function(bedFile) {
  interFile <- gsub("\\.bed$", "_mapped_woSecond_intersectS.bed", bedFile)
  
  # output name
  outFile <- gsub("\\.bed$", "_woSecond_intersectS_EXONuncover.bed", bedFile)
  
  
  # Reading files
  inter <- fread(interFile, header = FALSE)
  setnames(inter, c("chr", "start", "end", "sampleID", "score", "strand", 
                    "CDS_chr", "CDS_start", "CDS_end", "CDS_name", 
                    "CDS_score", "CDS_strand", "n_base"))
  
  # Identify blocks not overlapping exons (CDS_name == ".")
  uncovered_keys <- inter[CDS_name == ".", paste(chr, start, end, sampleID, sep = ":")]
  
  if(length(uncovered_keys) == 0) return(NULL)
  
  bed <- fread(bedFile, header = FALSE)
  bed_keys <- bed[, paste(V1, V2, V3, V4, sep = ":")]
  
  # Match and extract
  uncovered_bed <- bed[bed_keys %in% uncovered_keys]
  
  # Write part output
  fwrite(uncovered_bed, file = outFile, sep = "\t", col.names = FALSE)
  return(outFile)
}


results <- mclapply(baseFiles, process_pair, mc.cores = ncores)

