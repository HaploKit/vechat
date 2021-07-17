
<!-- ![Aaron Swartz](https://raw.githubusercontent.com/xiaoluo91/vechat/master/logo2.png) -->

<img src=https://raw.githubusercontent.com/xiaoluo91/vechat/master/logo2.png width=30% />

Haplotype-aware Error Correction for Noisy Long Reads Using Variation Graphs

## Description

Error self-correction is usually the first step and thus crucial in long-read sequencing data analysis. Current methods generally perform self-correction by computing a consensus sequence for each read, which may lose haplotype-specific variations in  heterogeneous genomes such as polyploid genomes and metagenomes. Here, we present vechat, a novel approach to perform haplotype-aware error correction for noisy long reads using variation graphs. Noisy variation graphs are firstly constructed from raw long reads, preserving true [v]()ariations and sequencing [e]()rrors. These graphs could be then pruned based on the correlations([chat]()) between nodes (i.e. variations or errors). The pruned variation graphs are subsequently used to perform haplotype-aware error correction.

## Installation and dependencies
vechat relies on the following dependencies:
- [racon](https://github.com/lbcb-sci/racon)
- [minimap2](https://github.com/lh3/minimap2)
- [fpa](https://github.com/natir/fpa)

To install vechat, firstly, it is recommended to intall the dependencies through [Conda](https://docs.conda.io/en/latest/):
```
conda create -n phasebook python=3.7
conda activate phasebook
conda install -c bioconda whatshap=0.18 minimap2 longshot samtools bcftools racon fpa=0.5
```

Subsequently, pull down the code to the directory where you want to install, and compile the code:
```
git clone https://github.com/phasebook/phasebook.git
cd phasebook
sh install.sh
```

## Running and options

The input read file is only required and the format should be FASTA or FASTQ. Other parameters are optional.
Please run `python phasebook.py -h` to get details of optional parameters setting. 
The final polished haplotype aware contigs are included in the `contigs.fa` file under output directory.

Before running phasebook, please read through the following basic parameter settings, 
which may be helpful to obtain better assemblies. Note that the option `-x` indicates 
using preset parameters for assembly, which is recommended.
```
positional arguments:
  sequences             input file in FASTA/FASTQ format (can be compressed
                        with gzip) containing sequences used for correction
  target_sequences      input file in FASTA/FASTQ format (can be compressed
                        with gzip) containing sequences which will be
                        corrected

optional arguments:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        output file (default: reads.corrected.fa)
  --platform PLATFORM   sequencing platform: pb/ont (default: pb)
  --split               split target sequences into chunks (recommend for
                        FASTQ > 20G or FASTA > 10G) (default: False)
  --split-size SPLIT_SIZE
                        split target sequences into chunks of desired size in
                        lines, only valid when using --split (default:
                        1000000)
  -d MIN_CONFIDENCE, --min-confidence MIN_CONFIDENCE
                        minimum confidence for keeping edges in the graph
                        (default: 0.2)
  -s MIN_SUPPORT, --min-support MIN_SUPPORT
                        minimum support for keeping edges in the graph
                        (default: 0.2)
  -t THREADS, --threads THREADS
                        number of threads (default: 1)
```

## Examples
One can test the program using the small PacBio HiFi reads file `example/reads.fq`.
`-g` is used to set the running mode for small or large genomes. If set `-g large`, 
it will utilize more efficient approaches for read overlap calculation and filtering,
 as well as sequencing error correction, but may at the cost of assembly performance.
 In general,

For small genomes or genomic regions assembly:
- PacBio CLR reads
```
cd example
vechat reads.fq.gz reads.fq.gz -t 8 --platform pb -o reads.corrected.fa 
```
- ONT reads
```
vechat reads.fq.gz reads.fq.gz -t 8 --platform ont -o reads.corrected.fa 
```

For running large datasets on a single local machine, one could add `--split` to reduce memory usage:
```
vechat reads.fq reads.fq -t 48 --platform pb --split -o reads.corrected.fa 
```


For running large datasets on HPC, one could refer to `scripts/vechat_hpc.sh` or `scripts/vechat_hpc.fast.sh` for how to split and submit jobs.


## Citation
