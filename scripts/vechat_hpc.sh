
platform=pb
iter=1
size=100000
#size=1000000  #default

#************************************#
file_mode=fq
threads=56

#######################################
raw_read=reads.raw.$file_mode
ln -fs ../../data/reads.$file_mode $raw_read

corrected_read=reads.corrected.fa

python /prj/whatshap-denovo/github/HapRacon/scripts/vechat.py $raw_read $raw_read --split   --split-size $size -t $threads --iter $iter --platform $platform -o tmp.fa
#python /prj/whatshap-denovo/github/HapRaconRobin/scripts/vechat.py $raw_read $raw_read --split   --split-size $size -t $threads --iter $iter --platform $platform -o tmp.fa

#filter very short reads <1000bp
/prj/whatshap-denovo/bin/filter_fa tmp.fa 1000 > $corrected_read

#rm -f tmp.fa
