#!/bin/bash
set -uo pipefail

echo "Command: $0 $@"

### Functions ###
usage() {
    echo "IFDlong perform isoform-level annotation of long-read RNA-seq data, detect gene fusions, and quantify both fusions and isoforms."
    printf "\n"
    echo "Usage: bash IFDlong.sh -o output_directory -n sample_name -i input_file -l 'self_align' -g 'hg38' -t 9 -a 10 -c 1"
    echo "Options:"
    echo "  -h, --help        Check the usage."
    echo "  -o, --outDir      The directory to save the output."
    echo "  -n, --name        The sample name."
    echo "  -i, --inFile      Input FASTQ file of long read RNAseq data OR BAM file after alignment and indexing."
    echo "  -l, --aligner     The aligner used to generate the BAM file; set it to self_align if the input format is BAM."
    echo "  -g, --ghc         Human (hg38), mouse (mm10) or other self-defined species (the same value as -g in refDataSetup.sh), hg38 by default"
    echo "  -t, --bufferLen   The buffer length for novel isoform identification, 9 by default."
    echo "  -a, --anchorLen   The anthor length for fusion filtering, 10 by default."
    echo "  -c, --ncores      How many cores are assigned to run the pipeline in parallel. Use 4 core by default"
    printf "\n"
    printf "\n"
    echo "Questions or issues? Contact: Silvia (shl96[at].pitt.edu)"
    #echo "Modified date: 08 Aug. 2025"
    #echo "Modified date: 21 Nov. 2025"
    echo "Modified date: 05 Feb. 2026"
}


symlink_path () {
    SD="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
    echo $SD
}


align () {
    echo "### The process the input ###"
    if [[ -z "$inFile" ]]; then
      echo "Usage: $0 <input_file.fq|input_file.bam>"
      exit 1
    fi
    
    # Get file extension
    ext="${inFile##*.}"
    
    filename=$(basename "$inFile")
    if [[ "$filename" == *.fq.gz ]]; then
        ext="fq.gz"
    elif [[ "$filename" == *.fastq.gz ]]; then
        ext="fastq.gz"
    else
        ext="${inFile##*.}"
    fi
    
    
    if [[ "$ext" == "fastq" || "$ext" == "fq" || "$ext" == "fa" || "$ext" == "fasta" || "$ext" == "fq.gz" || "$ext" == "fastq.gz" ]]; then 
      echo "Detected FASTQ file: $inFile"
      
      echo Begin Alignment by minimap2 $(date '+%Y-%m-%d %H:%M:%S')
      Aligner="MINIMAP2"
      outPath=$mainPath/$Aligner
      BAMfile=$outPath/$sample.bam

      echo "Output Path: $outPath"
      echo "Aligner:     $Aligner"

      mkdir -p "$outPath"
    
      $minimap2 -ax splice "$genomeDir" "$inFile" -t "$ncores" | samtools view -Sb | samtools sort -o $BAMfile
      if [[ $? -ne 0 ]]; then
          echo "[ERROR] Alignment pipeline failed (out of memory)!" >&2
          exit 1
      fi
      
      echo Alignment Done $(date '+%Y-%m-%d %H:%M:%S')
      echo samtools index
      $samtools index $BAMfile 
      echo samtools index Done $(date '+%Y-%m-%d %H:%M:%S')
    
    elif [[ "$ext" == "bam" ]]; then
      echo "Detected BAM file: $inFile"
      Aligner="self_align"
      outPath=$mainPath/$Aligner

      mkdir -p "$outPath"

      BAMfile=$inFile
    
    else
      echo "[ERROR] Unsupported file format: $ext"
      exit 1
    fi

}


filter () {
    echo filter out unmapped and multiple alignment for $sample $(date '+%Y-%m-%d %H:%M:%S')
    $samtools view -b -F 4 "$BAMfile" | \
    $samtools view -b -F 256 - > "$outPath/${sample}_mapped_woSecond.bam"
    # rm "$outPath/$sample"_mapped.bam"
    echo Finish Filtering for $sample $(date '+%Y-%m-%d %H:%M:%S')

    echo Generate BED intersect for $sample $(date '+%Y-%m-%d %H:%M:%S')
    $bedtools bamtobed -i "$outPath/${sample}_mapped_woSecond.bam" -split -cigar > $outPath/${sample}"_mapped_woSecond.bed"
    echo Finish convert the bam to bed file $(date '+%Y-%m-%d %H:%M:%S')

}

split_bed_by_readname () {
    local bedfile="$outPath/${sample}_mapped_woSecond.bed"
    local ncores="$ncores"
    local splitDir="$outPath/split_${sample}"

    echo "Splitting BED into $ncores parts"
    mkdir -p "$splitDir"

    awk -v T="$ncores" -v OUT="$splitDir/${sample}_part" '
    {
        read = $4
        if (!(read in idx)) {
            idx[read] = (count % T)
            count++
        }
        print $0 >> (OUT idx[read] ".bed")
    }
    ' "$bedfile"

    echo "Splitting done. Files:"
    ls "$splitDir"/*.bed
}

process_split_beds () {
    local splitDir="$outPath/split_${sample}"

    for bed in "$splitDir"/*.bed; do
        base=$(basename "$bed" .bed)

        echo "Processing $base"

        $bedtools intersect \
            -a "$bed" \
            -b "$refFile" \
            -f 0.90 -wao \
            > "$splitDir/${base}_mapped_woSecond_intersectS.bed"
    done
    echo Intersecting Finish for $sample $(date '+%Y-%m-%d %H:%M:%S')
}


blocks () {
    echo Begin EXON-uncovered blocks generating $(date '+%Y-%m-%d %H:%M:%S')
    $Rscript $EXONuncover $mainPath $sample $Aligner $refFile $ncores
    echo EXON-uncovered blocks generated in bed file Done!
    
    echo Begin gene range covering
    local splitDir="$outPath/split_${sample}"
    for bed in "$splitDir"/*_woSecond_intersectS_EXONuncover.bed; do
        filename=$(basename "$bed" .bed)
        base=${filename%%_woSecond_intersectS_EXONuncover}
    
        echo "Processing $base"
    
        # Run bedtools intersect
        bedtools intersect \
            -a "$bed" \
            -b "$refFiletol" \
            -f 0.90 -wao \
            > "$splitDir/${base}_mapped_woSecond_geneTol500intersectS.bed"
    done
    echo Gene range covering done!
}


anno () {
    echo Begin isoform annotation $(date '+%Y-%m-%d %H:%M:%S')
    $Rscript $report $mainPath $sample $Aligner $buffer $anchorLen $refFile $refAAFile $refPseudoFile $refRootFile $hmmatchFile $ghc $ncores
    echo Isoform annotation done! $(date '+%Y-%m-%d %H:%M:%S')
}


quant () {
    echo Begin isoform quantification $(date '+%Y-%m-%d %H:%M:%S')
    $Rscript $quant $mainPath $sample $Aligner $buffer $anchorLen $refGTFFile $ncores
    echo Isoform quantification done!$(date '+%Y-%m-%d %H:%M:%S')

    local splitDir="$outPath/split_${sample}"

    rm -r $splitDir
}





### Initialization ###
# Check if no arguments were provided
if [ $# -eq 0 ]; then
    usage
    exit 0

fi

######## default parameter
#Aligner="self_align"
sample=""
inFile=""
ghc="hg38"
ncores=1
mainPath="output"
codeBase="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

#echo $mainPath 
#echo $codeBase

### Argument Parsing ###
while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)
      usage
      exit 0
      ;;

    -o|--outDir)
      mainPath="$2"
      shift 2
      ;;

    -n|--name)
      sample="$2"
      shift 2
      ;;

    -i|--inFile)
      inFile="$2"
      shift 2
      ;;

    -l|--aligner)
      Aligner="$2"
      shift 2
      ;;

    -g|--ghc)
      ghc="$2"
      shift 2
      ;;

    -t|--bufferLen)
      bufferLen="$2"
      shift 2
      ;;

    -a|--anchorLen)
      anchorLen="$2"
      shift 2
      ;;

    -c|--ncores)
      ncores="$2"
      shift 2
      ;;

    *)
      echo "Invalid option: $1"
      usage
      ;;
  esac
done

### Required Argument Validation ###
if [ -z "$inFile" ]; then
  echo "Missing required option: -i|--inFile (input file)"
  missing_arg=true
fi

if [ -z "${sample:-}" ]; then
  sample=$(basename "$inFile")
  sample="${sample%%.*}"
  echo "Sample name not provided. Using '$sample' from input file."
fi

# Ensure output directory is set
if [ -z "${mainPath:-}" ]; then
  mainPath="output"
  echo "Output directory not provided. Using default: '$mainPath'."
fi

echo "Initialization complete."


##### check the dependencies
#Rscript=$(which Rscript)
#samtools=$(which samtools)
#bedtools=$(which bedtools)
#minimap2=$(which minimap2)

check_tool() {
    tool_name="$1"

    # Get all paths
    mapfile -t tool_paths < <(which -a "$tool_name" 2>/dev/null)

    # If none found
    if [ ${#tool_paths[@]} -eq 0 ]; then
        echo "Error: $tool_name not found in PATH."
        return 1
    fi

    # Choose the last path (latest one in PATH order)
    latest_path="${tool_paths[0]}"
    latest_version="0"

    for path in "${tool_paths[@]}"; do
        version=$("$path" --version 2>&1 | grep -oP '\d+(\.\d+)+' | head -n1)

        # Compare versions
        if [[ $(printf '%s\n%s\n' "$version" "$latest_version" | sort -V | tail -n1) == "$version" ]]; then
            latest_version="$version"
            latest_path="$path"
        fi
    done

    echo "[OK] $tool_name found at: $latest_path (version $latest_version)"
    eval "${tool_name}='$latest_path'"
}


# List of required tools
check_tool "Rscript" || exit 1
check_tool "samtools" || exit 1
check_tool "bedtools" || exit 1
check_tool "minimap2" || exit 1

### Check required files
if [ -f "$codeBase/refData/$ghc/genome.fa.gz" ] && [ ! -f "$codeBase/refData/$ghc/genome.fa" ]; then
    gunzip "$codeBase/refData/$ghc/genome.fa.gz"
fi

files=(
  "$codeBase/refData/$ghc/genome.fa"
  "$codeBase/refData/$ghc/genes.gtf"
  "$codeBase/refData/$ghc/isoformAA.txt"
  "$codeBase/refData/$ghc/pseudogenes.rds"
  "$codeBase/refData/$ghc/rootName.txt"
  "$codeBase/refData/$ghc/Hm_Mm_match.rds"
  "$codeBase/refData/$ghc/allexon_NO.bed"
  "$codeBase/refData/$ghc/gene_range_tol500.bed"
)

all_valid=true

min_size=1024
echo "Checking Reference Files"
for f in "${files[@]}"; do
    # Use -e to check if the path (or symlink target) exists
    if [ ! -e "$f" ]; then
        echo "[MISSING] $f"
        all_valid=false
    else
        # Get size of target file
        fsize=$(stat -L -c%s "$f" 2>/dev/null)
        
        # Check if fsize is empty (happens with broken links) or too small
        if [ -z "$fsize" ] || [ "$fsize" -lt "$min_size" ]; then
            echo "[INVALID] $(basename "$f") is too small or broken ($fsize bytes)."
            all_valid=false
        else
            echo "[OK] $(basename "$f") ($((fsize/1024)) KB)"
        fi
    fi
done

if [ "$all_valid" = true ]; then
    echo "[Success] All reference files exist and are valid."
else
    echo "[Error] Some files are missing, broken links, or invalid LFS pointers."
    echo "If these are 134-byte files, run 'git lfs pull'."
    exit 1
fi


### Main ###
echo "Pipeline Begin $(date '+%Y-%m-%d %H:%M:%S')"


######## ref file path
refFile=$codeBase/refData/$ghc/allexon_NO.bed
refFiletol=$codeBase/refData/$ghc/gene_range_tol500.bed
genomeDir=$(ls "$codeBase/refData/$ghc/"*.fa | head -n 1)
refGTFFile=$(ls "$codeBase/refData/$ghc/"*.gtf | head -n 1)

buffer="${bufferLen:-9}"
anchorLen="${anchorLen:-10}"
refAAFile=$codeBase/refData/$ghc/isoformAA.txt
refPseudoFile=$codeBase/refData/$ghc/pseudogenes.rds
refRootFile=$codeBase/refData/$ghc/rootName.txt
hmmatchFile=$codeBase/refData/$ghc/Hm_Mm_match.rds

## Path to scripts used by the IFDlong pipeline
#echo $codeBase
EXONuncover="${codeBase}/scripts/EXONuncover.R"
report="${codeBase}/scripts/Report.r"
quant="${codeBase}/scripts/quant.R"


symlink_path
align
filter
split_bed_by_readname
process_split_beds
blocks
anno
quant

if [ -f "$outPath/${sample}_mapped_woSecond_intersectS_buffer9bp_Isof_quant.csv" ]; then
    cp "$outPath/${sample}_mapped_woSecond_intersectS_buffer9bp_Isof_quant.csv" "$mainPath/${sample}_Isof_quant.csv"
else
    echo "[Warning] Isof quant results not found. If you encounter any issues, feel free to report them on GitHub."
fi

# Copy Fusion quant file if it exists
if [ -f "$outPath/${sample}_mapped_woSecond_intersectS_buffer9bp_Fusion_quant_anchor10bp.csv" ]; then
    cp "$outPath/${sample}_mapped_woSecond_intersectS_buffer9bp_Fusion_quant_anchor10bp.csv" "$mainPath/${sample}_Fusion_quant.csv"
else
    echo "[Warning] Fusion quant results not found. If you encounter any issues, feel free to report them on GitHub."
fi




echo "Completed!!!!!"

echo LOG END $(date '+%Y-%m-%d %H:%M:%S')

echo Used memory $(free |grep Mem|awk '{print $3}') 
echo Used memory percentage $(free |grep Mem|awk '{print $3/$2 * 100.0}')


exit 0



