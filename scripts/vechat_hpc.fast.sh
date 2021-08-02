#!/bin/bash

raw_read=$PWD/reads.fq
platform=pb
threads=48
binpath='/prj/whatshap-denovo/software/miniconda3/bin' #env path
outdir=$PWD/out
# n_lines=4000
# n_lines=200000 #1Gb fastq
n_lines=40000 #200M fastq

min_confidence=0.2
min_support=0.2
min_corrected_len=1000 
# min_identity_cns=0.99 #default
min_identity_cns=0.98 #high compleixty metagenome 20x
#######################################
# SCRIPTDIR=`dirname $0`
SCRIPTDIR='/prj/whatshap-denovo/github/vechat/scripts'
hapracon=$SCRIPTDIR/../build/bin/racon
corrected_read=$outdir/reads.corrected.fa

# num_lines=`cat $raw_read|wc -l`
# n_lines=`expr $num_lines / $n_split`

rm -rf $outdir 
mkdir  $outdir 

#split reads, TODO:fa/fq because minimap2 only accept .fa/.fq
split -l $n_lines -d --additional-suffix .tmp.fq $raw_read $outdir/reads_chunk 

#run the first round 
for target_read in $outdir/reads_chunk*.tmp.fq
do
    overlap=$target_read.paf
    echo "$binpath/minimap2 -x ava-$platform --dual=yes  $target_read $raw_read -t $threads |awk '\$11>=500' |$binpath/fpa drop --same-name --internalmatch  - > $overlap; "
done >run_overlap1.sh 

for target_read in $outdir/reads_chunk*.tmp.fq
do
    overlap=$target_read.paf
    # query_read=$raw_read
    query_read=$target_read.query.fq
    echo -n "perl -e 'my%h;open A,\$ARGV[1] or die; while(<A>){my@a=split;\$h{\$a[0]}=1; \$h{\$a[5]}=1;}close A; open A,\$ARGV[0] or die; open O,\">\$ARGV[2]\" or die; my\$flag=0;while(<A>){if(\$.%4==1){chomp;s/^@//;if(exists \$h{\$_}){\$flag=1;print O \"@\".\"\$_\\n\";}else{\$flag=0;} }elsif(\$flag){print O \$_;} }close A;close O;' $raw_read  $overlap  $query_read; "
    echo -n "$hapracon -f -p -d $min_confidence -s $min_support -t $threads  $query_read  $overlap  $target_read >$target_read.corrected.tmp.fa;"
    echo -n "$SCRIPTDIR/filter_fa $target_read.corrected.tmp.fa $min_corrected_len >$target_read.corrected.fa; "
    # echo "rm -f $overlap $target_read $target_read.corrected.tmp.fa;"
    # echo "rm -f $overlap $target_read.corrected.tmp.fa;"
    echo ""
done >run_round1.sh 
exit 

#submit to HPC 

split -l 1 -d  run_overlap1.sh sub-veovlp
for i in `ls sub-veovlp*`;do qsub -cwd -P fair_share -S /bin/bash -l arch=lx-amd64 -l h_rt=100000:00:00,h_vmem=80G,vf=80G $i;done

split -l 1 -d  run_round1.sh sub-1r
for i in `ls sub-1r*`;do qsub -cwd -P fair_share -S /bin/bash -l arch=lx-amd64 -l h_rt=100000:00:00,h_vmem=200G,vf=200G $i;done

#check if finished successfully manually !!!

#####################################
######## run the second round #######
#####################################

#Note that in the second round, it will lost many overlaps if splitting reads when running minimap2.
for i in $outdir/reads_chunk*corrected.fa
do
    cat $i 
done >$outdir/reads.round1.fa  


for target_read in $outdir/reads_chunk*corrected.fa
do
    overlap=$target_read.paf
    echo "$binpath/minimap2 -x ava-$platform --dual=yes  $target_read $outdir/reads.round1.fa -t $threads |awk '\$11>=1000 && \$10/\$11>=$min_identity_cns' |cut -f 1-12|$binpath/fpa drop --same-name --internalmatch  - > $overlap; "
done >run_overlap2.sh 

for target_read in $outdir/reads_chunk*corrected.fa
do
    overlap=$target_read.paf
    # query_read=$raw_read
    query_read=$target_read.query.fq
    echo -n "perl -e 'my%h;open A,\$ARGV[1] or die; while(<A>){my@a=split;\$h{\$a[0]}=1; \$h{\$a[5]}=1;}close A; open A,\$ARGV[0] or die; open O,\">\$ARGV[2]\" or die; my\$flag=0;while(<A>){if(\$.%2==1){chomp;s/^>//;my@a=split; if(exists \$h{\$a[0]}){\$flag=1;print O \">\".\"\$a[0]\\n\";}else{\$flag=0;} }elsif(\$flag){print O \$_;} }close A;close O;' $outdir/reads.round1.fa  $overlap  $query_read; "
    # echo -n "$hapracon -f -u -t $threads  $outdir/reads.round1.fa   $overlap  $target_read >$target_read.corrected.tmp;"
    echo -n "/prj/whatshap-denovo/software/miniconda3/bin/racon -f -u -t $threads  $outdir/reads.round1.fa   $overlap  $target_read >$target_read.corrected.tmp;" #use racon installed by conda to run on HPC

    echo -n "$SCRIPTDIR/filter_fa $target_read.corrected.tmp $min_corrected_len >$target_read.corrected2.fa; "
    # echo "rm -f $overlap $target_read.corrected.tmp;"
    echo;
done >run_round2.sh 

#submit to HPC 

split -l 1 -d  run_overlap2.sh sub-ovlp
for i in `ls sub-ovlp*`;do qsub -cwd -P fair_share -S /bin/bash -l arch=lx-amd64 -l h_rt=100000:00:00,h_vmem=100G,vf=100G $i;done

split -l 1 -d  run_round2.sh sub-2r
for i in `ls sub-2r*`;do qsub -cwd -P fair_share -S /bin/bash -l arch=lx-amd64 -l h_rt=100000:00:00,h_vmem=200G,vf=200G $i;done


#check if finished successfully manually !!!

#merge
for i in $outdir/reads_chunk*corrected2.fa
do
    cat $i 
done >$corrected_read 

rm -f $outdir/reads_chunk*

