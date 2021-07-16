
raw_read=''
platform=pb
threads=48
binpath='/prj/whatshap-denovo/software/miniconda3/bin' #env path
outdir='./out'
n_split=100 
min_confidence=0.2
min_support=0.2
min_corrected_len=1000 
#######################################
SCRIPTDIR=`dirname $0`
hapracon=$SCRIPTDIR/../build/bin/racon
corrected_read=$outdir/reads.corrected.fa

num_lines=`cat $raw_read|wc -l`
n_lines=`expr $num_lines / $n_split`

rm -rf $outdir 
mkdir  $outdir 

#split reads
split -l $n_lines -d --additional-suffix .tmp $raw_read $outdir/reads_chunk

#run the first round 
for target_read in $outdir/reads_chunk*.tmp
do
    overlap=$target_read.paf
    echo -n $binpath/minimap2 -x ava-$platform --dual=yes $raw_read $target_read -t $threads |awk '$11>=500' |$binpath/fpa drop --same-name --internalmatch  - > $overlap; 
    echo -n hapracon -f -p -d $min_confidence -s $min_support -t $threads  $raw_read  $overlap  $target_read >$target_read.corrected.tmp.fa;
    echo $SCRIPTDIR/filter_fa $target_read.corrected.tmp.fa $min_corrected_len >$target_read.corrected.fa; rm -f $overlap $target_read $target_read.corrected.tmp.fa;
done >run_round1.sh 

#submit to HPC 


#run the second round 
for i in $outdir/reads_chunk*corrected.fa
do
    cat $i 
done >$outdir/reads.round1.fa  

for target_read in $outdir/reads_chunk*corrected.fa
do
    overlap=$target_read.paf
    echo -n $binpath/minimap2 -x ava-$platform --dual=yes $outdir/reads.round1.fa $target_read -t $threads |awk '$11>=1000 && $10/$11>=0.99' |cut -f 1-12|$binpath/fpa drop --same-name --internalmatch  - > $overlap; 
    echo -n hapracon -f  -t $threads  $raw_read  $overlap  $target_read >$target_read.corrected.tmp;
    echo $SCRIPTDIR/filter_fa $target_read.corrected.tmp $min_corrected_len >$target_read.corrected2.fa; rm -f $overlap $target_read $target_read.corrected.tmp;
done 

#submit to HPC 

#merge
for i in $outdir/reads_chunk*corrected2.fa
do
    cat $i 
done >$corrected_read 

rm -f $outdir/reads_chunk*corrected2.fa

