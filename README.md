# vechat
Haplotype-aware error correction for noisy long reads using variation graphs

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
-i INFILE, --infile INFILE
                        input file in FASTA/FASTQ format (default: None)
-o OUTDIR, --outdir OUTDIR
                        output directory (default: .)
-t THREADS, --threads THREADS
                        number of threads (default: 1)
-p PLATFORM, --platform PLATFORM
                        sequencing platform(PacBio CLR/PacBio HiFi/Oxford
                        Nanopore): [pb/hifi/ont] (default: pb)
-x PRESET, --preset PRESET
                        use preset parameters
-g GENOMESIZE, --genomesize GENOMESIZE
                        genome size: small/large (default: small)
--overlaps OVERLAPS   input file in PAF format (default: None)
--min_cov MIN_COV     min coverage for trimming consensus (default: 4.0)
--min_identity MIN_IDENTITY
                        min identity for filtering overlaps (default: 0.75)
--min_read_len MIN_READ_LEN
                        min read length for processing (default: 1000)
--min_sread_len MIN_SREAD_LEN
                        min seed read length (default: 1000)
--min_ovlp_len MIN_OVLP_LEN
                        min overlap length for super reads construction
                        (default: 1000)
--n_correct N_CORRECT
                        times for self error correction of raw reads (default:
                        0)
--n_polish N_POLISH   times for super reads polishing (default: 2)
--sp_min_identity SP_MIN_IDENTITY
                        super reads min identity for filtering overlaps
                        (default: 0.98)
--min_cluster_size MIN_CLUSTER_SIZE
                        min size of read clusters (default: 4)
--trim_ends TRIM_ENDS
                        trim the erroneous bases in both ends, should be
                        either True or False (default: False)
--ctg_asm CTG_ASM    method to assemble super reads: [rb/naive], rb is time consuming, 
                        which is only recommended for small genomes (default: rb)
--correct_mode CORRECT_MODE
                        method to correct raw reads: [msa/hybrid], msa is much
                        faster than hybrid, which is recommended for large
                        genomes (default: msa)
--max_het_snps MAX_HET_SNPS
                        maximum number of heterozygous SNPs to determine the
                        contig overlap is from the identical haplotype or not
                        (default: 0)
--min_allele_cov MIN_ALLELE_COV
                        number of observations of each allele (default: 4)
--n_final_polish N_FINAL_POLISH
                        polish times for final contigs (default: 1)
```

## Examples
One can test the program using the small PacBio HiFi reads file `example/reads.fa`.
`-g` is used to set the running mode for small or large genomes. If set `-g large`, 
it will utilize more efficient approaches for read overlap calculation and filtering,
 as well as sequencing error correction, but may at the cost of assembly performance.
 In general,

For small genomes or genomic regions assembly (roughly, size < 50Mbp):
- PacBio HiFi reads
```
    cd example
    python ../scripts/phasebook.py -i reads.fa -t 8 -p hifi -g small -x 
```
- PacBio CLR reads
```
    python phasebook.py -i reads.fa -t 8 -p pb -g small -x 
```

- ONT reads
```
    python phasebook.py -i reads.fa -t 8 -p ont -g small -x 
```

For large genomes or genomic regions assembly:
- PacBio HiFi reads
```
    python phasebook.py -i reads.fa -t 8 -p hifi -g large -x 
```
- PacBio CLR reads
```
    python phasebook.py -i reads.fa -t 8 -p pb -g large -x 
```
- ONT reads
```
    python phasebook.py -i reads.fa -t 8 -p ont -g large -x 
```

## Possible issues during installation (optional)

- If `g++` version of the system is not satisfied, one could try this to install:
```
conda install -c conda-forge gxx_linux-64=7.3.0
# replace the /path/to/ with your own path
ln -s /path/to/miniconda3/envs/phasebook/bin/x86_64-conda-cos6-linux-gnu-g++ /path/to/miniconda3/envs/phasebook/bin/g++
ln -s /path/to/miniconda3/envs/phasebook/bin/x86_64-conda-cos6-linux-gnu-gcc /path/to/miniconda3/envs/phasebook/bin/gcc
```
- If `boost` or `zlib` library is not installed, one could try this to install:
```
conda install -c conda-forge boost zlib
# set envionment variables
export LD_LIBRARY_PATH=/path/to/miniconda3/envs/phasebook/lib/:$LD_LIBRARY_PATH
export CPATH=/path/to/miniconda3/envs/phasebook/include/:$CPATH
```

- If compile error occurs something like `/path/to/miniconda3/envs/phasebook/x86_64-conda_cos6-linux-gnu/bin/ld: cannot find -lboost_timer `
or `cannot find -lz`, 
 which means it fails to link `boost` or `zlib` library, one could try this to solve:
```
ln -s /path/to/miniconda3/envs/phasebook/lib/libboost_* /path/to/miniconda3/envs/phasebook/x86_64-conda_cos6-linux-gnu/lib/
ln -s /path/to/miniconda3/envs/phasebook/lib/libz.* /path/to/miniconda3/envs/phasebook/x86_64-conda_cos6-linux-gnu/lib/
# then re-complile and install
sh install.sh
```

## Citation
