d=`date`
echo ---------------------------------------------------------------------------------------------
echo Running LINKS on ABySS E.coli assembly with Full 2D -High Quality- Oxford Nanopore long reads
echo --------------------------------------------------------------------------------------------- 
echo WITH CURRENT PARAMETERS, YOU NEED A SYSTEM WITH AT LEAST 8GB RAM. INCREASE -t IF YOU DO NOT
echo Table 1C and 1D in manuscript 
echo ---------------------------------------------------------------------------------------------
wget ftp://ftp.bcgsc.ca/supplementary/LINKS/K12_full2dONT_longread.fa
wget ftp://ftp.bcgsc.ca/supplementary/LINKS/ecoliK12_abyss_illumina_contig_baseline.fa
wget ftp://ftp.bcgsc.ca/supplementary/LINKS/ecoliK12_abyss_illumina_scaffold_baseline.fa
echo ---------------------------------------------------------------------------------------------
echo done. Initiating LINKS scaffolding ETA 1-2 min depending on system...
echo ---------------------------------------------------------------------------------------------
echo K12_full2dONT_longread.fa > K12_F2D.fof
/usr/bin/time -v -o timeLINKS_ECK12singleTIG.txt ../LINKS -f ecoliK12_abyss_illumina_contig_baseline.fa -s K12_F2D.fof -b ecoliK12-ONT_linksSingleIterationTIG
/usr/bin/time -v -o timeLINKS_ECK12singleSCAFF.txt ../LINKS -f ecoliK12_abyss_illumina_scaffold_baseline.fa -s K12_F2D.fof -b ecoliK12-ONT_linksSingleIterationSCAFF
echo ---------------------------------------------------------------------------------------------
echo done. Initiating LINKS scaffolding with iterative distances ETA 3 min depending on system...
echo ---------------------------------------------------------------------------------------------
/usr/bin/time -v -o timeLINKS_ECK12singleSCAFFrange.txt ../LINKS -f ecoliK12_abyss_illumina_scaffold_baseline.fa -s K12_F2D.fof -d 1000,2000,3000,4000,5000,6000,7000,8000,9000,10000,12000,14000,16000,18000,20000 -b ecoliK12-ONT_linksSingleIterationSCAFFrange
echo ---------------------------------------------------------------------------------------------
