### RAN: ./LINKSrecipe_athaliana_raw.sh 21 5 5 0.3 at unix prompt - will require ~84GB RAM, performed on machine with 264GB RAM for manuscript (Fig S7 and S8 in manuscript)
### http://1001genomes.org/data/MPI/MPISchneeberger2011/releases/current/Ler-1/Assemblies/Allpaths_LG/Ler-1.allpaths_lg.final.assembly.fasta
### pacbio.fof | ftp://qb.cshl.edu/schatz/ectools/arabidopsis/Pacbio.fasta.gz
### k = 21, l = 5, a = 0.3
../LINKS -f Ler-1.allpaths_lg.final.assembly.fasta -s pacbio.fof -b rawATlinks1 -d 5000 -t 20 -k $1 -l $3 -a $4 -o 10
../LINKS -f rawATlinks1.scaffolds.fa -s pacbio.fof -b rawATlinks2 -d 10000 -t $2 -k $1 -l $3 -a $4 -r rawATlinks1.bloom -o 10
../LINKS -f rawATlinks2.scaffolds.fa -s pacbio.fof -b rawATlinks3 -d 15000 -t $2 -k $1 -l $3 -a $4 -r rawATlinks1.bloom -o 10
../LINKS -f rawATlinks3.scaffolds.fa -s pacbio.fof -b rawATlinks4 -d 20000 -t $2 -k $1 -l $3 -a $4 -r rawATlinks1.bloom -o 10
