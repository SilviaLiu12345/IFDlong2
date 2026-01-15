
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
  
  inter <- fread(interbedfile, header = FALSE)
  setnames(inter, c("chr","start","end","SampleID","score","strand",
                    "CDS_chr","CDS_start","CDS_end","CDS_name",
                    "CDS_score","CDS_strand","n_base"))
  
  inter <- inter[CDS_name != "."]
  cds_split <- tstrsplit(inter$CDS_name, "__")
  
  inter[, `:=`(
    gene    = cds_split[[6]],
    isoform = cds_split[[1]],
    order   = cds_split[[7]])]
  
  if (!file.exists(intergenebedfile)) {
    return(inter[, .(chr,start,end,SampleID,strand,
                     CDS_chr,CDS_start,CDS_end,CDS_strand,
                     n_base,gene,isoform,order)])
  }
  
  uncov <- fread(intergenebedfile, header = FALSE)
  setnames(uncov, c("chr","start","end","SampleID","score","strand",
                    "gene_chr","gene_start","gene_end","gene_name",
                    "gene_score","gene_strand","n_base"))
  
  uncov[, gene := ifelse(gene_name == ".", "undefined",
                         tstrsplit(gene_name, "__", keep = 1))]
  uncov[, `:=`(
    isoform="undefined", order="undefined",
    CDS_chr="undefined", CDS_start=NA, CDS_end=NA,
    CDS_strand="undefined")]
  
  rbind(
    inter[, .(chr,start,end,SampleID,strand,
              CDS_chr,CDS_start,CDS_end,CDS_strand,
              n_base,gene,isoform,order)],
    uncov[, .(chr,start,end,SampleID,strand,
              CDS_chr,CDS_start,CDS_end,CDS_strand,
              n_base,gene,isoform,order)],
    fill = TRUE
  )
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


#### df2 fusion
Sys.time()
df2_summary <- df2[, .(
  gene = paste(unique(gene), collapse="&"),
  gene_strand = paste(unique(strand), collapse="&"),
  isoform = paste(unique(isoform), collapse="||"),
  position = paste(CDS_chr,CDS_start,CDS_end,CDS_strand,
                   sep=":", collapse=";"),
  nblock = .N,
  NO.Exon = paste(sort(unique(as.numeric(order))), collapse="-"),
  fusion = "Y",
  continuous = all(diff(sort(as.numeric(order))) == 1)
), by = SampleID]

df2_summary[, `:=`(
  note = fifelse(continuous,
                 "continuous CDS and edge-matching",
                 "discontinuous CDS"),
  type = fifelse(continuous,"normal","novel with deletion"))]
Sys.time()

Sys.time()
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
df2_summary <- add_iso_stats(df2_summary)
Sys.time()

### merge
final <- rbindlist(list(df1_summary, df2_summary), fill = TRUE)

collist <- c("SampleID","gene","gene_strand","isoform","position",
             "nblock","NO.Exon","nExon_isof","length_isof",
             "fusion","note","type")

final <- final[, ..collist]
fwrite(final, coverSReOut)
                             
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

fullRepOut <- paste0(PATH,"/",Aligner,"/",sampleName,
                     "_mapped_woSecond_intersectS_buffer",
                     buffer,"bp_fullRep.csv")

fwrite(final, fullRepOut)                            
                             
fusionfiltReOut <- paste0(PATH,"/",Aligner,"/",sampleName,
                          "_mapped_woSecond_intersectS_buffer",
                          buffer,"bp_fusionRep_anchor",
                          anchorLen,"bp.filt.csv")

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
  return(report)
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

