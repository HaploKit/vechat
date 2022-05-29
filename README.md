
![logo](logo.png)

Correcting Errors in Noisy Long Reads Using Variation Graphs

## Description

Error correction is the canonical first step in long-read sequencing data analysis. The current standard is to make use of a consensus sequence as a template. However, in mixed samples, such as metagenomes or organisms of higher ploidy, consensus induced biases can mask true variants affecting haplotypes of lower frequencies, because they are mistaken as errors.

The novelty presented here is to use graph based, instead of sequence based consensus as a template for identifying errors. The advantage is that graph based reference systems also capture variants of lower frequencies, so do not mistakenly mask them as errors. We present VeChat, as a novel approach to implement this idea: VeChat distinguishes errors from haplotype-specific true variants based on variation graphs, which reflect a popular type of data structure for pangenome reference systems. Upon initial construction of an ad-hoc variation graph from the raw input reads, nodes and edges that are due to errors are pruned from that graph by way of an iterative procedure that is based on principles from frequent itemset mining. Upon termination, the graph exclusively contains nodes and edges reflecting true sequential phenomena. Final re-alignments of the raw reads indicate where and how reads need to be corrected. VeChat is implemented with C++ and Python3.


## Installation and dependencies
VeChat relies on the following dependencies:
- [racon](https://github.com/lbcb-sci/racon)
- [minimap2](https://github.com/lh3/minimap2)
- [fpa](https://github.com/natir/fpa) and [yacrd](https://github.com/natir/yacrd)
- gcc 4.8+ or clang 3.4+
- cmake 3.2+
- zlib

#### 1.Install from [Conda](https://docs.conda.io/en/latest/) (TODO)
This is easy and recommended:
```
conda create -n vechat
conda activate vechat
conda install -c bioconda vechat
```

#### 2.Install from source code
To install VeChat, firstly, it is recommended to intall the dependencies through [Conda](https://docs.conda.io/en/latest/):
```
conda create -n vechat
conda activate vechat
conda install -c bioconda minimap2 yacrd fpa=0.5
```

## Running and options
If you pulled the source repo; to run vechat you just need to invoke the script with python intepreter in created conda enviornment from the repo directory
```bash
git clone https://github.com/HaploKit/vechat.git
cd vechat
python ./scripts/vechat -h
```
Modified `racon` dependency will be built automatically.

The input read file is only required and the format should be FASTA/FASTQ (can be compressed with gzip). Other parameters are optional.
Please run `vechat -h` to get details of optional arguments. 

```
positional arguments:
  sequences             input file in FASTA/FASTQ format (can be compressed
                        with gzip) containing sequences used for correction

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
  --scrub               scrub chimeric reads (default: False)
  --min-identity-cns MIN_IDENTITY_CNS
                        minimum sequence identity between read overlaps in the
                        consensus round (default: 0.99)
  -t THREADS, --threads THREADS
                        number of threads (default: 1)
```

## Examples
One can test the program using the small sequencing read file `example/reads.fq.gz`:

- Correcting PacBio CLR reads
```
cd example
vechat reads.fq.gz -t 8 --platform pb -o reads.corrected.fa 
```
- Correcting ONT reads
```
vechat reads.fq.gz -t 8 --platform ont -o reads.corrected.fa 
```

- For running large datasets on a single local machine, one could add `--split` to reduce memory usage:
```
vechat reads.fq.gz -t 48 --platform pb --split -o reads.corrected.fa 
```


- For running large datasets on HPC to speed up, one could refer to `scripts/vechat_hpc.fast.sh` for how to split and submit jobs.


## Citation
