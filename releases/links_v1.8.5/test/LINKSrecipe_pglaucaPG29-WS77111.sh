### RAN: ./LINKSrecipe_pglaucaPG29-WS77111.sh 26 200 5 0.3 at unix prompt - will require ~132GB RAM, performed on machine with 264GB RAM for manuscript (Fig 4 in manuscript) -- LINKS v1.1
### run6_scaffold500.fa is PG29-V3 (Genbank:ALWZ030000000)
### ws77111sealed1_500.fa is WS77111-V1 (Genbank:JZKD010000000)
### k = 26, l = 5, a = 0.3
../LINKS -f run6_scaffold500.fa -s ws77111sealed1_500.fa -b links1 -d 5000 -t 200 -k $1 -l $3 -a $4
../LINKS -f links1.scaffolds.fa -s ws77111sealed1_500.fa -b links2 -d 7500 -t 200 -k $1 -l $3 -a $4
../LINKS -f links2.scaffolds.fa -s ws77111sealed1_500.fa -b links3 -d 10000 -t 150 -k $1 -l $3 -a $4
../LINKS -f links3.scaffolds.fa -s ws77111sealed1_500.fa -b links4 -d 12500 -t 150 -k $1 -l $3 -a $4
../LINKS -f links4.scaffolds.fa -s ws77111sealed1_500.fa -b links5 -d 15000 -t 150 -k $1 -l $3 -a $4
../LINKS -f links5.scaffolds.fa -s ws77111sealed1_500.fa -b links6 -d 20000 -t 100 -k $1 -l $3 -a $4
../LINKS -f links6.scaffolds.fa -s ws77111sealed1_500.fa -b links7 -d 30000 -t 75 -k $1 -l $3 -a $4
../LINKS -f links7.scaffolds.fa -s ws77111sealed1_500.fa -b links8 -d 40000 -t 50 -k $1 -l $3 -a $4
../LINKS -f links8.scaffolds.fa -s ws77111sealed1_500.fa -b links9 -d 50000 -t 50 -k $1 -l $3 -a $4
../LINKS -f links9.scaffolds.fa -s ws77111sealed1_500.fa -b links10 -d 60000 -t 50 -k $1 -l $3 -a $4
../LINKS -f links10.scaffolds.fa -s ws77111sealed1_500.fa -b links11 -d 70000 -t 50 -k $1 -l $3 -a $4
../LINKS -f links11.scaffolds.fa -s ws77111sealed1_500.fa -b links12 -d 80000 -t 50 -k $1 -l $3 -a $4
../LINKS -f links12.scaffolds.fa -s ws77111sealed1_500.fa -b links13 -d 90000 -t 50 -k $1 -l $3 -a $4
../LINKS -f links13.scaffolds.fa -s ws77111sealed1_500.fa -b links14 -d 100000 -t 50 -k $1 -l $3 -a $4
