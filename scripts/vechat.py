#!/usr/bin/env python

from __future__ import print_function
import os
import sys
import time
import argparse
import subprocess
import gzip


def run_error_correction(
        sequences, chunk_target_sequence,
        include_unpolished, platform,split,
        fragment_correction, linear, min_confidence,
        min_support, corrected_file,
        window_length, quality_threshold,
        error_threshold, match, mismatch, gap, threads,
        cudaaligner_batches, cudapoa_batches, cuda_banded_alignment):

    racon_path = os.path.dirname(
        os.path.abspath(__file__))+"/../build/bin/racon"
    overlap='overlap.paf'
    try:
        # compute overlap and filter
        os.system("minimap2 -x ava-{} --dual=yes {} {} -t {} 2>/dev/null|awk '$11>=500'|\
            fpa drop --same-name --internalmatch  - >{}"
                  .format(platform, chunk_target_sequence, sequences, threads,overlap))
    except:
        raise Exception("Unable to compute overlaps!")
    
    #extract only query sequences which have overlaps with target sequences to reduce memory
    if split:
        sub_query_sequences = extract_sub_sequences(sequences,overlap,chunk_target_sequence)
    else:
        sub_query_sequences = sequences

    try:
        if not linear:
            # perform haplotype aware error correction
            print("perform variation graph based (haplotype-aware) error correction")
            os.system("{} -f -p -d {} -s {} -t {} {} overlap.paf {} >{}"
                      .format(racon_path, min_confidence, min_support, threads, sub_query_sequences,
                              chunk_target_sequence, corrected_file))
        else:
            # perform normal error correction, which is the same with orginal 'racon'
            print("perform linear sequence based error correction")
            os.system("{} -f  -t {} {} overlap.paf {} >{}"
                      .format(racon_path, threads, sub_query_sequences,
                              chunk_target_sequence, corrected_file))
    except:
        raise Exception("Unable to run error correction!")
    return corrected_file

def extract_sub_sequences(sequences,overlap,chunk_target_sequence):
    mode=fq_or_fa(chunk_target_sequence)
    target_reads={}
    query_reads={}
    if mode == 'fq':
        i=0
        with open(chunk_target_sequence) as fr:
            for line in fr:
                i+=1
                if(i%4==1):
                    target_reads[line[1:].rstrip().split()[0]]=1
                else:
                    continue
    else:
        i=0
        with open(chunk_target_sequence) as fr:
            for line in fr:
                i+=1
                if(i%2==1):
                    target_reads[line[1:].rstrip().split()[0]]=1
                else:
                    continue
    with open(overlap) as fr:
        for line in fr:
            a=line.split()
            if a[0] in target_reads or a[5] in target_reads:
                query_reads[a[0]]=1
                query_reads[a[5]]=1

    sub_query_records=[]
    if_extract=False
    if mode == 'fq':
        i=0
        with open(sequences) as fr:
            for line in fr:
                if(i%4==0):
                    if line[1:].rstrip().split()[0] in query_reads:
                        if_extract=True 
                        sub_query_records.append(line)
                    else:
                        if_extract=False
                elif(i%4==1):
                    if if_extract:
                        sub_query_records.append(line)
                elif(i%4==2):
                    if if_extract:
                        sub_query_records.append(line)
                else:
                    if if_extract:
                        sub_query_records.append(line)
                i+=1 

    if mode == 'fa':
        i=0
        with open(sequences) as fr:
            for line in fr:
                if(i%2==0):
                    if line[1:].rstrip().split()[0] in query_reads:
                        if_extract=True 
                        sub_query_records.append(line)
                    else:
                        if_extract=False
                else:
                    if if_extract:
                        sub_query_records.append(line)
                i+=1   

    sub_query_sequences='query_sequences.tmp.{}'.format(mode)
    with open(sub_query_sequences,'w') as fw:
        fw.write(''.join(sub_query_records))
    return sub_query_sequences

def fq_or_fa(file):
    if file.endswith('.gz'):
        fr = gzip.open(file, 'rt')
    else:
        fr = open(file, 'r')
    s = fr.readline()[0]
    mode = ''
    if s == '>':
        mode = 'fa'
    elif s == '@':
        mode = 'fq'
    else:
        raise Exception(
            "invalid input file, must be FASTA/FASTQ format.", file)
    return mode


# *******************************************************************************


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='''Variation Graph based Long Read Error Correction''',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('sequences', help='''input file in FASTA/FASTQ format
        (can be compressed with gzip) containing sequences used for correction''')
    # parser.add_argument('overlaps', help='''input file in MHAP/PAF/SAM format
    #     (can be compressed with gzip) containing overlaps between sequences and
    #     target sequences''')
    parser.add_argument('target_sequences', help='''input file in FASTA/FASTQ
        format (can be compressed with gzip) containing sequences which will be
        corrected''')

    parser.add_argument('-o', '--outfile', required=False, default='reads.corrected.fa',
                        help="output file")
    parser.add_argument('--platform', default='pb',
                        help='''sequencing platform: pb/ont''')
    parser.add_argument('--split', action='store_true',
                        help='''split target sequences into chunks (recommend for FASTQ > 20G or FASTA > 10G)''')
    parser.add_argument('--split-size', default=1000000,type=int, help='''split target sequences into chunks of
        desired size in lines, only valid when using --split''')

    parser.add_argument('-u', '--include-unpolished', action='store_true',
                        help='''output unpolished target sequences''')
    parser.add_argument('-f', '--fragment-correction', action='store_false',
                        help='''perform fragment correction instead of contig polishing
        (overlaps file should contain dual/self overlaps!)''')

    parser.add_argument('--linear', action='store_true',
                        help='''perform linear based fragment correction rather than variation graph based
                        fragment correction''')
    parser.add_argument('-d', '--min-confidence', default=0.2,
                        help='''minimum confidence for keeping edges in the graph''')
    parser.add_argument('-s', '--min-support', default=0.2,
                        help='''minimum support for keeping edges in the graph''')
    # parser.add_argument('-k', '--num-prune', default=3,
    #                     help='''number of iterations for pruning the graph''')
    parser.add_argument('--iter', default=1,type=int,
                        help='''number of iterations for error correction''')
    parser.add_argument('-w', '--window-length', default=500, help='''size of
        window on which POA is performed''')
    parser.add_argument('-q', '--quality-threshold', default=10.0,
                        help='''threshold for average base quality of windows used in POA''')
    parser.add_argument('-e', '--error-threshold', default=0.3, help='''maximum
        allowed error rate used for filtering overlaps''')
    parser.add_argument('-t', '--threads', default=1,
                        help='''number of threads''')

    parser.add_argument('-m', '--match', default=5, help='''score for matching
        bases''')
    parser.add_argument('-x', '--mismatch', default=-4, help='''score for
        mismatching bases''')
    parser.add_argument('-g', '--gap', default=-8, help='''gap penalty (must be
        negative)''')
    parser.add_argument('--cudaaligner-batches', default=0,
                        help='''number of batches for CUDA accelerated alignment''')
    parser.add_argument('-c', '--cudapoa-batches', default=0,
                        help='''number of batches for CUDA accelerated polishing''')
    parser.add_argument('-b', '--cuda-banded-alignment', action='store_true',
                        help='''use banding approximation for polishing on GPU. Only applicable when -c is used.''')

    args = parser.parse_args()

    if(args.split):
        try:
            suffix = fq_or_fa(args.target_sequences)

            query_sequences_file = ''
            target_sequences_file = ''
            corrected_file = ''
            for i in range(1, args.iter+1):
                print("Performing the {} iteration for error correction...".format(i))

                if i == 1:
                    query_sequences_file = args.sequences
                    target_sequences_file = args.target_sequences
                    #need gsplit if on MacOS
                    os.system("split -l {} -d --additional-suffix .{} {} reads_chunk".format(
                        args.split_size, suffix, args.target_sequences))
                    chunk_target_sequences = os.popen(
                        "ls reads_chunk*").read().strip().split()
                elif i > 1:
                    new_split_size = args.split_size if suffix == 'fa' else int(
                        args.split_size/2)
                    query_sequences_file = "reads.corrected.tmp{}.fa".format(
                        i-1)
                    target_sequences_file = "reads.corrected.tmp{}.fa".format(
                        i-1)

                    os.system("split -l {} -d --additional-suffix .{} {} reads_chunk".format(
                        new_split_size, 'fa', target_sequences_file))
                    chunk_target_sequences = os.popen(
                        "ls reads_chunk*.fa").read().strip().split()
                else:
                    raise ValueError("Invalid iteration: {}".format(i))

                corrected_file = "reads.corrected.tmp{}.fa".format(i)
                corrected_chunk_files = []
                j = 0
                for chunk_target_sequence in chunk_target_sequences:
                    j += 1
                    print("processing chunk {}...".format(j))
                    corrected_chunk_file = "reads.corrected.tmp.chunk{}.fa".format(
                        j)
                    run_error_correction(
                        query_sequences_file, chunk_target_sequence,
                        args.include_unpolished, args.platform,args.split,
                        args.fragment_correction, args.linear, args.min_confidence,
                        args.min_support, corrected_chunk_file,
                        args.window_length, args.quality_threshold,
                        args.error_threshold, args.match, args.mismatch, args.gap, args.threads,
                        args.cudaaligner_batches, args.cudapoa_batches, args.cuda_banded_alignment)
                    corrected_chunk_files.append(corrected_chunk_file)

                # os.system("echo -n >{}".format(corrected_file))
                open(corrected_file,'w').close()
                for corrected_chunk_file in corrected_chunk_files:
                    os.system("cat {} >>{}".format(
                        corrected_chunk_file, corrected_file))
                    os.system("rm -f {} ".format(corrected_chunk_file))

            os.system("mv {} {}".format(corrected_file, args.outfile))
            os.system("rm -f reads.corrected.tmp*.fa")
            os.system("rm -f reads_chunk*")

        except OSError:
            eprint('[RaconWrapper::run] error: unable to run racon!')
            sys.exit(1)
    else:
        try:
            query_sequences_file = ''
            target_sequences_file = ''
            corrected_file = ''
            for i in range(1, args.iter+1):
                print("Performing the {} iteration for error correction...".format(i))
                if i == 1:
                    query_sequences_file = args.sequences
                    target_sequences_file = args.target_sequences
                elif i > 1:
                    query_sequences_file = "reads.corrected.tmp{}.fa".format(
                        i-1)
                    target_sequences_file = "reads.corrected.tmp{}.fa".format(
                        i-1)
                else:
                    raise ValueError("Invalid iteration: {}".format(i))
                corrected_file = "reads.corrected.tmp{}.fa".format(i)
                run_error_correction(
                    query_sequences_file, target_sequences_file,
                    args.include_unpolished, args.platform,args.split,
                    args.fragment_correction, args.linear, args.min_confidence,
                    args.min_support, corrected_file,
                    args.window_length, args.quality_threshold,
                    args.error_threshold, args.match, args.mismatch, args.gap, args.threads,
                    args.cudaaligner_batches, args.cudapoa_batches, args.cuda_banded_alignment)
            os.system("mv {} {}".format(corrected_file, args.outfile))
            os.system("rm -f reads.corrected.tmp*.fa")
        except OSError:
            eprint('[RaconWrapper::run] error: unable to run racon!')
            sys.exit(1)
