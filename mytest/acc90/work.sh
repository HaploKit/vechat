
#check the error rate(defined by Racon) of overlap
cat overlaps.paf |perl -lne 'my@a=split;my$i=($a[3]-$a[2]);my$j=($a[8]-$a[7]);my$min=$i;my$max=$j; if($i>$j){$min=$j;$max=$i;};my$r=1-$min/$max;print "Err: $r";'|sort -k2nr|le

#filter overlaps
scp vincent@dolores:/export/scratch3/vincent/project/vglr/vglr_result/clr/3-HIV/acc90/overlaps.paf .

./racon -t 4 -f reads.fa overlaps.paf reads.fa --include-unpolished --no-trimming >reads.racon.notrim.fa

./racon -t 4 -f reads.fa overlaps.paf reads.fa > reads.racon.fa 2> racon.log

./hapracon -t 4 -f  --haplotype reads.fa overlaps.paf reads.fa > reads.hapracon.fa 2> hapracon.log

./hapracon -t 4 -f --haplotype reads.fa overlaps.paf reads.fa --include-unpolished --no-trimming > reads.hapracon.notrim.fa 2>/dev/null

