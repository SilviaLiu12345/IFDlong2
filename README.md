# IFDlong #

## About ##

IFDlong is a bioinformatics pipeline that can perform long-read RNA-seq annotation at isoform levels, fusion detection, as well as fusion and isoform quantification.


## Installation ##
### Dependencies ###
R (≥ 4.4.0) along with a compatible version of gcc (12.2.0).  
The following tools are required for running IFDlong: minimap2 (≥ 2.30), bedtools (≥ 2.31), and samtools (≥ 1.17)  
If these tools are not already installed, they will be installed automatically during the install.sh step.
Required R packages include: rlist, parallel, stringr, dplyr, Rcpp, seqRFLP, BiocManager, rtracklayer, and Biostrings.
These packages will also be installed automatically during processing.

### Installation pipeline ###
Method 1: Install via conda
```
git clone https://github.com/SilviaLiu12345/IFDlong2.git
cd IFDlong2
conda env create -f IFDlong.yaml
conda activate IFDlong
```

Method 2: Install maunally  
It will automatically install the dependent tools and put their paths to the `tools.path` file.  
```bash
git clone https://github.com/SilviaLiu12345/IFDlong2.git
cd IFDlong2
bash Install.sh
```


## Reference Database ##
Our tools include built-in support for human (hg38) and mouse (mm10) reference datasets.
The corresponding genome.fa and genes.gtf files from the UCSC Genome Browser [illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html). 

Method 1: Install via git lfs  
Before installation, make sure **Git LFS** is [installed](https://docs.github.com/en/repositories/working-with-files/managing-large-files/installing-git-large-file-storage).  
```
cd IFDlong2
git lfs install
git lfs pull
```
Method 2: Install maunally  
From the GitHub repository `refData/(hg38 or mm10)`, copy the following files into the corresponding folder.
```
cd IFDlong2/refData/(hg38 or mm10)
- Copy `genome.fa.gz`
- Copy `genes.gtf`
```

## Usage ##
### Quick Start
```bash
bash IFDlong.sh -o output_directory -n sample_name -i input_file -l "self_align" -g "hg38" -t 9 -a 10
```

#### Demo
To verify your installation and run a example.

```
bash IFDlong.sh -o out -n example -i example/demo.fq.gz -g "hg38" -t 9 -a 10
```

### Run IFDlong Pipeline
Our tool accepts both FASTQ and BAM files as input. If using FASTQ, it will be aligned using Minimap2. If using BAM, please ensure the files are already aligned and indexed.

Running with default settings:
```bash
bash IFDlong.sh -o output_directory -n sample_name -i input_file -l "self_align" -g "hg38" -t 9 -a 10
```

Required options:
```
-h, --help        Check the usage.
-o, --outDir      Output directory to save results
-n, --name        Sample name prefix
-i, --inFile      Input file (supported formats: .fq, .fq.gz, .fastq, .fastq.gz, .fa, .fa.gz, .fasta, .fasta.gz, .bam)
-l, --aligner     The aligner used to generated the bam file. Set to be self_align if missing.
-g, --ghc         Human (hg38), mouse (mm10) or other self-defined species (the same value as -g in refDataSetup.sh), hg38 by default
-t, --bufferLen   The buffer length for novel isoform identification, 9 by default
-a, --anchorLen   The anthor length for fusion filtering, 10 by default
-c, --ncores      How many cores are assigned to run the pipeline in parallel. Use 4 core by default

```

## OUTPUT FORMAT ##
Two files, **\[sampleID]\_Isof\_quant.csv** and **\[sampleID]\_Fusion\_quant.csv**, will be generated in the output folder.
1. Isof_quant.csv contains the following columns:
isoform, gene, group, prop, count
2. Fusion_quant.csv contains the following columns:
fusion, gene, group, prop, count, fusion_counts, gene_counts



## Citation ##
The study describing the IFDlong method can be found in: 

## License ##
The software is under the MIT License. Please see the LICENSE file for details.

