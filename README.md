# IFDlong #

## About ##

IFDlong is a probabilistic framework and software suite for detecting isoform and fusion transcripts from bulk or single-cell long-RNA-seq data. IFDlong annotates each long read, identifies novel isoforms, quantifies expression via an expectation-maximization algorithm, and profiles fusion transcripts. 


## Installation ##
### Dependencies ###
R (version ≥ 4.4.0) is required, along with a compatible version of gcc (12.2.0).

The following external tools are required to run IFDlong: minimap2 (≥ 2.24), bedtools (≥ 2.31), and samtools (≥ 1.17). If these tools are not already installed, they will be installed automatically during the Install.sh step.

Required R packages include data.table, parallel, stringr, rlist, dplyr, purrr, tidyr. These packages will also be installed automatically during the installation or processing steps.


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
source path/to/IFDlong2/tools.path
```


## Reference Database ##
The tool includes built-in support for human (hg38) and mouse (mm10) reference datasets. Corresponding genome.fa and genes.gtf files are downloaded from the UCSC Genome Browser [illumina iGenomes](https://support.illumina.com/sequencing/sequencing_software/igenome.html). 

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
- Copy `isoformAA.txt`
```

## Usage ##
### Demo
To verify your installation and run an example.

```
bash IFDlong.sh -o out -n example -i example/demo.fq.gz -g "hg38" -t 9 -a 10
```

### Run IFDlong Pipeline
The tool accepts both FASTQ and BAM files as input. FASTQ files will be aligned using Minimap2. For BAM input, please ensure that the files are already aligned and indexed.

Running with default settings:
```bash
bash IFDlong.sh -o output_directory -n sample_name -i input_file -l "self_align" -g "hg38" -t 9 -a 10 -c 1
```

Required options:
```
-h, --help        Check the usage.
-o, --outDir      Output directory to save results
-n, --name        Sample name prefix for output files
-i, --inFile      Input file (supported formats: .fq, .fq.gz, .fastq, .fastq.gz, .fa, .fa.gz, .fasta, .fasta.gz, .bam)
-l, --aligner     The aligner used to generate the bam file. Set to be 'self_align' if missing.
-g, --ghc         Human (hg38), mouse (mm10), hg38 by default
-t, --bufferLen   The buffer length in base pairs for novel isoform identification, 9 by default
-a, --anchorLen   The anchor length in base pairs for fusion identification, 10 by default
-c, --ncores      Number of cores are assigned to run the pipeline in parallel. Use 1 core by default

```

## OUTPUT FORMAT ##
Two files, **\[sampleID]\_Isof\_quant.csv** and **\[sampleID]\_Fusion\_quant.csv**, will be generated in the output folder.
1. Isof_quant.csv contains the following columns:
isoform, gene, group, prop, count
2. Fusion_quant.csv contains the following columns:
fusion, gene, group, prop, count, isoform_counts, gene_counts



## Citation ##
The study describing the IFDlong method can be found in:  

Wang W, Li Y, Ko S, Feng N, Zhang M, Liu JJ, Zheng S, Ren B, Yu YP, Luo JH, Tseng GC, Liu S. IFDlong: an isoform and fusion detector for accurate annotation and quantification of long-read RNA-seq data. bioRxiv [Preprint]. 2024 May 14:2024.05.11.593690. doi: 10.1101/2024.05.11.593690. PMID: 38798496; PMCID: PMC11118288.

## License ##
The software is under the MIT License. Please see the LICENSE file for details.

