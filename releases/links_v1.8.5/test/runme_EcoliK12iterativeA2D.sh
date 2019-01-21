d=`date`
echo ---------------------------------------------------------------------------------------------
echo Running LINKS on baseline E. coli K12 MG1655 assembly with 2D Oxford Nanopore long reads
echo Table 1F in manuscript 
echo ---------------------------------------------------------------------------------------------
echo *Must be running on a computer with at least 16 GB RAM, with current parameters
echo ---------------------------------------------------------------------------------------------
echo Downloading baseline assembly and ONT long reads... be patient.
wget http://www.bcgsc.ca/downloads/supplementary/LINKS/K12_all2dONT_longread.fa.gz
wget http://www.bcgsc.ca/downloads/supplementary/LINKS/K12illumina-baseline.fa
gunzip K12_all2dONT_longread.fa.gz
echo ---------------------------------------------------------------------------------------------
echo Running LINKS iteratively a few times...be patient!
echo K12_all2dONT_longread.fa > K12_A2D.fof
/usr/bin/time -v -o timeLINKS_ECK12A2D.txt ./runIterativeLINKS_ECK12A2D.sh 15 2
