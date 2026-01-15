
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
  
  install_and_load(c("data.table", "parallel", "stringr", "rlist", "dplyr", "purrr"))
})


#### parameter
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) >= 12)

PATH <- args[1]
sampleName <- args[2]
Aligner <- args[3]
buffer <- as.integer(args[4])
anchorLen <- as.integer(args[5])
refEXON <- args[6] #"./refData/allexon_NO.bed"
refAA <- args[7] #"./refData/isoformAA_Ref.txt"
refPseudo <- args[8] #"./refData/pseudogenes.rds"
refRoot <- args[9] #"./refData/rootNames.txt"
refHMmatch <- args[10] #"./refData/Hm_Mm_match.rds"
species <- args[11] #hg38 or mm10
ncores <- as.integer(args[12])


### ref files
allCDS <- fread(refEXON,
                col.names = c("chr","start","end","name","score","strand"))
allCDS[, isoform := tstrsplit(name, "__", keep = 1)]
setkey(allCDS, chr, isoform)

isoformAA <- fread(refAA)
isoformAA[, isoformID := tstrsplit(isoformID, "__", keep = 1)]
setkey(isoformAA, isoformID)

read_if_exists <- function(path, fun) if (file.exists(path)) fun(path) else NULL
pseudogenes <- read_if_exists(refPseudo, readRDS)
rootNames   <- read_if_exists(refRoot, fread)
Hm_Mm_match <- read_if_exists(refHMmatch, readRDS)

#### input files
interSbedfile=paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS.bed")
intergenebedfile=paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_geneTol500intersectS.bed")
coverSReOut=paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_Rep.csv")


uncoverFilter <- function(interbedfile, intergenebedfile) {
  intersect_tab <- read.table(interbedfile, stringsAsFactors = FALSE)
  colnames(intersect_tab) <- c("chr", "start", "end", "SampleID", "score", "strand", 
                               "CDS_chr", "CDS_start", "CDS_end", "CDS_name", 
                               "CDS_score", "CDS_strand", "n_base")
  
  match_intersect <- subset(intersect_tab, CDS_name != ".")
  CDS_info <- strsplit(match_intersect$CDS_name, "__")
  
  match_info <- match_intersect[, c(1:4, 6:9, 12, 13)]
  match_info$gene    <- sapply(CDS_info, `[[`, 6)
  match_info$isoform <- sapply(CDS_info, `[[`, 1)
  match_info$order   <- sapply(CDS_info, `[[`, 7)
  
  if (inherits(try(read.table(intergenebedfile), silent = TRUE), "try-error")) {
    read_info <- match_info
  } else {
    uncover_info <- read.table(intergenebedfile, stringsAsFactors = FALSE)
    colnames(uncover_info) <- c("chr", "start", "end", "SampleID", "score", "strand", 
                                "gene_chr", "gene_start", "gene_end", "gene_name", 
                                "gene_score", "gene_strand", "n_base")
    
    unmatch_info <- uncover_info[, c(1:4, 6:9, 12, 13)]
    unmatch_info$gene <- ifelse(uncover_info$gene_name != ".", 
                                sapply(strsplit(uncover_info$gene_name[uncover_info$gene_name != "."], "__"), `[[`, 1), 
                                "undefined")
    
    unmatch_info$isoform <- "undefined"
    unmatch_info$order   <- "undefined"
    colnames(unmatch_info) <- c("chr", "start", "end", "SampleID", "strand", 
                                "CDS_chr", "CDS_start", "CDS_end", "CDS_strand", 
                                "n_base", "gene", "isoform", "order")
    
    read_info <- rbind(match_info, unmatch_info)
  }
  
  read_info$CDS_strand[read_info$CDS_strand == "."] <- "undefined"
  read_info$CDS_chr[read_info$CDS_chr == "."]       <- "undefined"
  
  message(" Extract the CDS Covered Alignments Done!")
  return(as.data.frame(read_info))
}

match_info <- uncoverFilter(interSbedfile, intergenebedfile)
setDT(match_info)

match_info[, gene    := as.character(gene)]
match_info[, isoform := as.character(isoform)]
match_info[, order   := as.character(order)]

print(dim(match_info))

###### split df1 and df2
gene_counts <- match_info[gene != "undefined",
                          .(n_gene = uniqueN(gene)),
                          by = SampleID]

match_info <- match_info[gene != "undefined"]
match_info <- gene_counts[match_info, on = "SampleID"]

df1 <- match_info[n_gene == 1][, n_gene := NULL]
df2 <- match_info[n_gene > 1][, n_gene := NULL]

print(dim(df1))
print(dim(df2))

### iso ref
isoform_summary <- allCDS[, .(
  nExon_isof = .N,
  length_isof = sum(end - start)
), by = isoform]
setkey(isoform_summary, isoform)

######### df1 
Sys.time()
df1_summary <- df1[, .(
  gene = first(gene),
  gene_strand = first(strand),
  isoform = paste(unique(isoform), collapse="||"),
  nblock = .N,
  NO.Exon = paste(sort(unique(as.numeric(order))), collapse="-"),
  position = paste(CDS_chr,CDS_start,CDS_end,CDS_strand,
                   sep=":", collapse=";"),
  fusion = "N",
  continuous = all(diff(sort(as.numeric(order))) == 1)
), by = SampleID]

df1_summary[, `:=`(
  note = fifelse(continuous,
                 "continuous CDS and edge-matching",
                 "discontinuous CDS"),
  type = fifelse(continuous,"normal","novel with deletion"))]
Sys.time()

### func
add_iso_stats <- function(dt) {
  dt[, iso_list := strsplit(isoform, "\\|\\|")]
  dt[, nExon_isof := sapply(iso_list,
                            function(x) paste(isoform_summary[x, nExon_isof], collapse="||"))]
  dt[, length_isof := sapply(iso_list,
                             function(x) paste(isoform_summary[x, length_isof], collapse="||"))]
  dt[, iso_list := NULL]
  dt
}

df1_summary <- add_iso_stats(df1_summary)
Sys.time()


###### df2 fusion
block_match <- function(chr_val, start_val, end_val, isoform_val, ref_df, buffer=9) {
  ref_sub <- ref_df[chr == chr_val & isoform == isoform_val]
  if (nrow(ref_sub) == 0) return(FALSE)
  any((start_val >= (ref_sub$start - buffer) & start_val <= (ref_sub$end + buffer)) &
        (end_val >= (ref_sub$start - buffer) & end_val <= (ref_sub$end + buffer)))
}


Sys.time()
if (nrow(df2) == 0) {
  message("No fusion reads detected.") 
} else {
  
  Sys.time()
  # ######## the correct one, but slowerrrr
  df2_summary <- df2 %>%
    group_by(SampleID, gene) %>%
    group_modify(~{
      
      x <- .x
      
      # --- build isoform-specific info ---
      iso_info <- x %>%
        group_split(isoform) %>%
        map_df(function(g){
          
          ord <- as.numeric(g$order)
          
          tibble(
            isoform     = unique(g$isoform),
            cds_len     = sum(g$CDS_end - g$CDS_start + 1),
            continuous  = all(diff(sort(ord)) == 1),
            nblock      = length(ord),
            exon_order  = paste(sort(ord), collapse = "-"),
            pos_sig     = paste(
              g$CDS_chr, g$CDS_start, g$CDS_end,
              g$CDS_strand, sep = ":", collapse = ";"
            )
          )
        })
      
      # --- isoform selection rules ---
      if (any(iso_info$continuous)) {
        
        iso_sub <- iso_info %>% filter(continuous)
        
        # longest CDS among continuous
        max_len <- max(iso_sub$cds_len)
        iso_sub <- iso_sub %>% filter(cds_len == max_len)
        
        # same long & same pos → keep them all
        chosen_pos <- iso_sub$pos_sig[1]
        chosen_isoforms <- iso_sub %>% filter(pos_sig == chosen_pos) %>% pull(isoform)
        
      } else {
        # fallback to longest discontinuous
        max_len <- max(iso_info$cds_len)
        chosen_isoforms <- iso_info %>% filter(cds_len == max_len) %>% pull(isoform)
      }
      
      chosen_rows <- x %>% filter(isoform %in% chosen_isoforms)
      
      tibble(
        SampleID   = unique(x$SampleID),
        gene       = unique(x$gene),
        gene_strand= first(x$strand),
        
        isoform    = paste(chosen_isoforms, collapse="||"),
        
        nblock     = n_distinct(chosen_rows$order),
        
        NO.Exon    = paste(
          sort(unique(as.numeric(chosen_rows$order))),
          collapse = "-"
        ),
        
        position   = paste(
          paste(x$CDS_chr, x$CDS_start, x$CDS_end, x$CDS_strand, sep = ":"),
          collapse = ";"
        ),
        
        nExon_isof = NA,
        length_isof = NA,
        fusion = "Y",   # df2 is fusion → changed to Y automatically
        
        note = if (any(iso_info$continuous)) {
          "continuous CDS and edge-matching"
        } else {
          "discontinuous CDS"
        },
        
        type = if (any(iso_info$continuous)) {
          "normal"
        } else {
          "novel with deletion"
        }
      )
    }) %>%
    ungroup()
  Sys.time()
  
  dim(df2_summary)
  
  df2_summary_cont <- df2_summary%>%filter(note == "continuous CDS and edge-matching")
  dim(df2_summary_cont)
  #head(df2_summary_cont)
  
  
  Sys.time()
  df2_summary_cont <- df2_summary_cont%>%
    rowwise() %>%
    mutate(
      # Split position string into individual blocks as a list
      blocks = list(str_split(position, ";")[[1]]),
      # Check each block against reference
      matched = all(sapply(blocks[[1]], function(b) {
        parts <- str_split(b, ":")[[1]]
        chr <- parts[1]
        start <- as.integer(parts[2])
        end <- as.integer(parts[3])
        block_match(chr, start, end, isoform, allCDS)
      })),
      # Update note if any block does not match
      note = ifelse(matched, note, "continuous CDS and edge-unmatching")
    ) %>%
    ungroup() %>%
    select(-blocks, -matched)
  Sys.time()
  
  df2_summary[df2_summary$note == "continuous CDS and edge-matching", ] <- df2_summary_cont
  
  table(df2_summary$note)
  
  
  
  #### add nExon_isof and length_isof cols
  df2_summary <- df2_summary %>%
    rowwise() %>%
    mutate(
      nExon_isof = if_else(
        isoform == "undefined", 
        NA_character_,
        paste(
          str_split(isoform, "\\|\\|")[[1]] %>%
            sapply(function(x) {
              val <- isoform_summary$nExon_isof[isoform_summary$isoform == x]
              if(length(val) == 0) NA else val
            }),
          collapse = "||"
        )
      ),
      length_isof = if_else(
        isoform == "undefined",
        NA_character_,
        paste(
          str_split(isoform, "\\|\\|")[[1]] %>%
            sapply(function(x) {
              val <- isoform_summary$length_isof[isoform_summary$isoform == x]
              if(length(val) == 0) NA else val
            }),
          collapse = "||"
        )
      )
    ) %>%
    ungroup()
  
  head(df2_summary)
  
  
  ########### v2
  head(df2_summary)
  
  df_counts <- df2_summary%>%
    group_by(SampleID)%>%
    mutate(n_rows = n())
  
  dim(df_counts)
  
  # df1: groups where row count > 2
  df_counts21 <- df_counts%>%
    filter(n_rows > 2) %>%
    select(-n_rows)
  
  dim(df_counts21)
  
  combine_all_cols <- function(df) {
    # Ensure all columns have names
    if (any(names(df) == "")) {
      names(df)[names(df) == ""] <- paste0("V", seq_len(sum(names(df) == "")))
    }
    
    # Group by position
    pos_groups <- df %>%
      group_by(position) %>%
      summarise(across(everything(), ~list(.x)), .groups = "drop")
    
    # If only 1 unique position, combine all rows with &
    if (nrow(pos_groups) == 1) {
      combined <- pos_groups %>%
        mutate(across(everything(), ~paste(.x[[1]], collapse = "&")))
      return(combined)
    }
    
    # All pairwise position combinations
    idx <- combn(nrow(pos_groups), 2)
    
    # Initialize a named list to hold final column values
    final_combined <- setNames(vector("list", length = ncol(pos_groups)), names(pos_groups))
    final_combined <- map(final_combined, ~character(0))
    
    for (k in seq_len(ncol(idx))) {
      i <- idx[, k]
      row1 <- pos_groups[i[1], ]
      row2 <- pos_groups[i[2], ]
      
      for (colname in names(pos_groups)) {
        cross <- expand.grid(row1[[colname]][[1]], row2[[colname]][[1]], stringsAsFactors = FALSE)
        combined <- paste0(cross$Var1, "&", cross$Var2)
        final_combined[[colname]] <- c(final_combined[[colname]], combined)
      }
    }
    
    # Collapse each column by #
    final_combined <- map(final_combined, ~paste(.x, collapse = "#"))
    
    # Return as flat tibble
    tibble::as_tibble(final_combined)
  }
  
  # Apply per SampleID
  if (nrow(df_counts21) == 0) {
    message(" ")
  } else {
    # Continue only when df_counts21 is not empty
    df2_summary_grouped1 <- df_counts21 %>%
      group_by(SampleID) %>%
      filter(n() >= 3) %>%
      group_modify(~ combine_all_cols(.x))
  }
  
  
  ###################### 
  # df2: groups where row count == 2
  df_counts22 <- df_counts%>%
    filter(n_rows == 2)%>%
    select(-n_rows)
  
  dim(df_counts22)
  
  ######### fusion 1x1 
  Sys.time()
  df2_summary_grouped <- df_counts22%>%
    group_by(SampleID) %>%
    summarise(
      # Check if positions are all identical
      same_position = n_distinct(position) == 1,
      
      # Merge each column accordingly
      gene = if(same_position) paste(gene, collapse = "#") else paste(gene, collapse = "&"),
      gene_strand = if(same_position) paste(gene_strand, collapse = "#") else paste(gene_strand, collapse = "&"),
      isoform = if(same_position) paste(isoform, collapse = "#") else paste(isoform, collapse = "&"),
      position = if(same_position) paste(position, collapse = "#") else paste(position, collapse = "&"),
      nblock = if(same_position) paste(nblock, collapse = "#") else paste(nblock, collapse = "&"),
      NO.Exon = if(same_position) paste(NO.Exon, collapse = "#") else paste(NO.Exon, collapse = "&"),
      nExon_isof = if(same_position) paste(nExon_isof, collapse = "#") else paste(nExon_isof, collapse = "&"),
      length_isof = if(same_position) paste(length_isof, collapse = "#") else paste(length_isof, collapse = "&"),
      note = if(same_position) paste(note, collapse = "#") else paste(note, collapse = "&"),
      type = if(same_position) paste(type, collapse = "#") else paste(type, collapse = "&"),
      
      # Fusion logic
      fusion = if(same_position) "N" else "Y",
      .groups = "drop"
    )
  
  head(df2_summary_grouped)
  table(df2_summary_grouped$fusion)
  
}

########## merge
collist <- c("SampleID","gene","gene_strand","isoform","position","nblock","NO.Exon","nExon_isof","length_isof","fusion","note","type")

dfs_to_bind <- list(df1_summary)

if (exists("df2_summary_grouped")) {
  dfs_to_bind <- c(dfs_to_bind, list(df2_summary_grouped))
}

if (exists("df2_summary_grouped1")) {
  dfs_to_bind <- c(dfs_to_bind, list(df2_summary_grouped1))
}

# Apply column selection and rbind
final <- do.call(
  rbind,
  lapply(dfs_to_bind, function(x) {
    # If it's a data.table, use the .. prefix to select columns via variable
    if (inherits(x, "data.table")) {
      return(x[, ..collist])
    } else {
      # If it's a data.frame or tibble, standard selection works
      return(x[, collist])
    }
  })
)

#head(com)
dim(final)

write.table(final, coverSReOut, row.names = FALSE, sep = ',')


#### AA part
AA_lookup <- isoformAA[, .(isoformID, AAseq, note)]
setkey(AA_lookup, isoformID)

annot_AA <- function(isof) {
  ids <- unlist(strsplit(isof, "[|&]"))
  hit <- AA_lookup[J(ids)]
  c(
    paste(fifelse(is.na(hit$AAseq),"noncoding",hit$AAseq), collapse="||"),
    paste(fifelse(is.na(hit$note),"noncoding",hit$note), collapse="||")
  )
}

AAres <- t(sapply(final$isoform, annot_AA))
final[, `:=`(AAseq = AAres[,1], AAnote = AAres[,2])]

fullRepOut <- paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_fullRep.csv")

fwrite(final, fullRepOut)                            
                             
fusionfiltReOut <- paste0(PATH,"/",Aligner,"/",sampleName,"_mapped_woSecond_intersectS_buffer",buffer,"bp_fusionRep_anchor",anchorLen,"bp.filt.csv")

partlength=function(parts) {
  partAs=sapply(parts,"[[",1)
  partBs=sapply(parts,"[[",2)
  blocksAs=strsplit(partAs,";")
  blocksBs=strsplit(partBs,";")
  partAslen=sapply(blocksAs,function(p) {
    blockAs=strsplit(p,":")
    blockAstart=as.numeric(sapply(blockAs,"[[",2))
    blockAend=as.numeric(sapply(blockAs,"[[",3))
    blocklen=abs(blockAend-blockAstart)
    partlen=sum(blocklen)
    return(partlen)
  })
  
  partBslen=sapply(blocksBs,function(p) {
    blockBs=strsplit(p,":")
    blockBstart=as.numeric(sapply(blockBs,"[[",2))
    blockBend=as.numeric(sapply(blockBs,"[[",3))
    blocklen=abs(blockBend-blockBstart)
    partlen=sum(blocklen)
    return(partlen)
  })
  partlen=cbind(partAslen,partBslen)
  colnames(partlen)=c("partA_len","partB_len")
  return(partlen)
}
      
filteredRep=function(reportPath,fusionfiltPath,min.len=10,pseudogenes,rootNames,species="hg38",Hm_Mm_match) {
  report=read.csv(reportPath)
  allannot=strsplit(report$position, "#")
  nannot=sapply(allannot,length)
  #print(table(nannot))
  allparts=strsplit(sapply(allannot,"[[",1), "&")
  nparts=sapply(allparts,length)
  #print(which(nparts>2))
  #print(table(nparts))
  
  #criteria I
  if (sum(nparts > 1) > 0) {
    
    if (sum(nparts==2)>0) {
      fusionPart=partlength(allparts[nparts==2])
      Rep_2p=report[nparts==2,]
      report$fusionlen="pass"
      report$fusionlen[(report$SampleID %in% Rep_2p$SampleID[apply(fusionPart,1,min)<=min.len])]="failed"
    }
    if (sum(nparts > 2)>0) {
      partslen.list=lapply(allparts[nparts > 2],function(x) {
        np=length(x)
        partslen=c()
        for (i in 1:np) {
          blocks=strsplit(x[i],";")
          partlen=sapply(blocks,function(p) {
            blocklist=strsplit(p,":")
            blockstart=as.numeric(sapply(blocklist,"[[",2))
            blockend=as.numeric(sapply(blocklist,"[[",3))
            blocklen=abs(blockend-blockstart)
            partlen=sum(blocklen)
            return(partlen)
          })
          partslen[i]=partlen
        }
        return(partslen)
      })
      Rep_ge3p=report[nparts > 2,]
      if ("fusionlen" %in% colnames(report)) {
        report$fusionlen[(report$SampleID %in% Rep_ge3p$SampleID[sapply(partslen.list,min)<=min.len])]="failed"
        report$fusionlen[!(report$SampleID %in% Rep_2p$SampleID|report$SampleID %in% Rep_ge3p$SampleID)]=NA
      }
      else {
        report$fusionlen="pass"
        report$fusionlen[(report$SampleID %in% Rep_ge3p$SampleID[sapply(partslen.list,min)<=min.len])]="failed"
        report$fusionlen[!(report$SampleID %in% Rep_ge3p$SampleID)]=NA
      }
    }
  }
  else {
    report$fusionlen=NA
  }
  
  #criteria II
  if (!is.null(pseudogenes)) {
    pseudoID=sapply(1:length(report$gene),function(i) {
      x=report$gene[i]
      annolist=strsplit(x,"#")
      annopseudo=sapply(annolist,function(g) {
        genes=unlist(strsplit(g,"&"))
        #npseudo=sum(genes %in% unique(pseudoName))
        npseudo=sum(genes %in% unique(pseudogenes))
      })
      return(all(annopseudo!=0))
    })
    report$pseudogene="pass"
    report$pseudogene[pseudoID]="failed" 
  }
  else {
    report$pseudogene="pass"
  }
  
  #criteria III
  if (!is.null(rootNames)) {
    if (sum(report$fusion!="N")>0) {
      if (species=="mm10") {
        familyID=sapply(1:sum(report$fusion!="N"),function(i) {
          x=report$gene[report$fusion!="N"][i]
          annolist=strsplit(x,"#")
          annofam=sapply(annolist,function(g) {
            genes=unlist(strsplit(g,"&"))
            Hm_match=keep(Hm_Mm_match,function(x) x$mouse %in% genes)
            if (length(Hm_match)==0) return(0)
            FAMlist=sapply(rootNames$Common.root.gene.symbol,function(g) {
              return(sum(sapply(Hm_match,function(i) length(str_which(unlist(i),paste0("^",g,"[^a-zA-Z]+"))))))
            })
            FAMid=sapply(FAMlist,function(j) {
              return(j==length(genes))
            })
            nfam=sum(unlist(FAMid))
            return(nfam)
          })
          return(all(annofam!=0))
        })
      }
      else {
        familyID=sapply(1:sum(report$fusion!="N"),function(i) {
          x=report$gene[report$fusion!="N"][i]
          annolist=strsplit(x,"#")
          annofam=sapply(annolist,function(g) {
            genes=unlist(strsplit(g,"&"))
            FAMlist=sapply(rootNames$Common.root.gene.symbol,function(g) grep(paste0("^",g,"[^a-zA-Z]+"),genes))
            FAMid=sapply(FAMlist,function(j) {
              return(length(j)==length(genes))
            })
            nfam=sum(unlist(FAMid))
          })
          return(all(annofam!=0))
        })
      }
      report$FamGene="pass"
      report$FamGene[report$fusion=="N"]=NA
      report$FamGene[report$fusion!="N"][familyID]="failed" #1481
    }
    
    else {
      report$FamGene=NA
    }
  }
  else {
    report$FamGene="pass"
    report$FamGene[report$fusion=="N"]=NA
  }
  
  write.csv(report,reportPath,row.names = F)
  write.csv(report[report$fusionlen!="failed" & report$pseudogene!="failed" & report$FamGene!="failed" & report$fusion!="N" ,],fusionfiltPath,row.names = F)
  invisible(report)
}


filteredRep(
  reportPath = fullRepOut,
  fusionfiltPath = fusionfiltReOut,
  min.len = anchorLen,
  pseudogenes = pseudogenes,
  rootNames = rootNames,
  species = species,
  Hm_Mm_match = Hm_Mm_match
)

print("finish!!!")

