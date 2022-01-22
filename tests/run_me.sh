echo simulated_genome.fa > reads.fof
/usr/bin/time -v -o time_installation_test.txt LINKS -f simulated_contigs.fa -s reads.fof -d 1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,12000,14000,16000,18000,20000 -b installation_test 
