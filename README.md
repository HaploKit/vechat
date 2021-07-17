
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
- gcc 4.8+ or clang 3.4+
- cmake 3.2+
- zlib

#### 1.Install from [Conda](https://docs.conda.io/en/latest/)
This is easy and recommended:
```
conda create -n vechat
conda activate vechat
conda install -c bioconda vechat
```

#### 2.Install from source code
To install vechat, firstly, it is recommended to intall the dependencies through [Conda](https://docs.conda.io/en/latest/):
```
conda create -n vechat
conda activate vechat
conda install -c bioconda minimap2 fpa=0.5
```

Subsequently, pull down the code to the directory where you want to install, and compile the code:
```
git clone https://github.com/xiaoluo91/vechat.git
cd vechat
mkdir build;cd build;cmake -DCMAKE_BUILD_TYPE=Release ..;make
```

Add the following code to your `$HOME/.bashrc`, please change the `/path/to/vechat` as your own path:
```
export PATH=/path/to/vechat/scripts/:$PATH
```

## Running and options

The input read file is only required and the format should be FASTA/FASTQ (can be compressed with gzip). Other parameters are optional.
Please run `vechat -h` to get details of optional arguments. 

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
One can test the program using the small sequencing read file `example/reads.fq.gz`:

- Correcting PacBio CLR reads
```
cd example
vechat reads.fq.gz reads.fq.gz -t 8 --platform pb -o reads.corrected.fa 
```
- Correcting ONT reads
```
vechat reads.fq.gz reads.fq.gz -t 8 --platform ont -o reads.corrected.fa 
```

- For running large datasets on a single local machine, one could add `--split` to reduce memory usage:
```
vechat reads.fq.gz reads.fq.gz -t 48 --platform pb --split -o reads.corrected.fa 
```


- For running large datasets on HPC to speed up, one could refer to `scripts/vechat_hpc.sh` or `scripts/vechat_hpc.fast.sh` for how to split and submit jobs.


## Citation
