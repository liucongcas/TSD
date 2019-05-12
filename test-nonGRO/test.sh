softwaredir=PATH-to-TSD
sh ${softwaredir}/TSD.tophat.sh -i $index -x test_1.fastq -y test_2.fastq -p 10 -w ${softwaredir} 
python ${softwaredir}/TSD.py -i test_config.txt -o test-filted
