../LINKS -f StyphiH58illumina-baseline.fa -s STH58_2D.fof -b STH58links1 -d 500 -t $2 -k $1 -l $3 -a $4
../LINKS -f STH58links1.scaffolds.fa -s STH58_2D.fof -b STH58links2 -r STH58links1.bloom -d 750 -t $2 -k $1 -l $3 -a $4
../LINKS -f STH58links2.scaffolds.fa -s STH58_2D.fof -b STH58links3 -r STH58links1.bloom -d 1000 -t $2 -k $1 -l $3 -a $4
../LINKS -f STH58links3.scaffolds.fa -s STH58_2D.fof -b STH58links4 -r STH58links1.bloom -d 1250 -t $2 -k $1 -l $3 -a $4
../LINKS -f STH58links4.scaffolds.fa -s STH58_2D.fof -b STH58links5 -r STH58links1.bloom -d 1500 -t $2 -k $1 -l $3 -a $4
../LINKS -f STH58links5.scaffolds.fa -s STH58_2D.fof -b STH58links6 -r STH58links1.bloom -d 1750 -t $2 -k $1 -l $3 -a $4
../LINKS -f STH58links6.scaffolds.fa -s STH58_2D.fof -b STH58links7 -r STH58links1.bloom -d 3000 -t $2 -k $1 -l $3 -a $4
../LINKS -f STH58links7.scaffolds.fa -s STH58_2D.fof -b STH58links8 -r STH58links1.bloom -d 3250 -t $2 -k $1 -l $3 -a $4
../LINKS -f STH58links8.scaffolds.fa -s STH58_2D.fof -b STH58links9 -r STH58links1.bloom -d 3500 -t $2 -k $1 -l $3 -a $4
../LINKS -f STH58links9.scaffolds.fa -s STH58_2D.fof -b STH58links10 -r STH58links1.bloom -d 3750 -t $2 -k $1 -l $3 -a $4
../LINKS -f STH58links10.scaffolds.fa -s STH58_2D.fof -b STH58links11 -r STH58links1.bloom -d 4000 -t $2 -k $1 -l $3 -a $4
