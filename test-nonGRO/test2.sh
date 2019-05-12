softwaredir=../
sh ${softwaredir}/TSD.tophat.sh -i /data/data1/liuc/reference_genome/hg19-bowtie-1.1.2-index/hg19 -x test_1.fastq -y test_2.fastq -p 10 -w ${softwaredir} 
python ${softwaredir}/TSD.py -i test_config2.txt -o test-filted
