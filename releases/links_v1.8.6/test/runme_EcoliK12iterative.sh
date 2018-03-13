d=`date`
echo ---------------------------------------------------------------------------------------------
echo Running LINKS on baseline E. coli K12 MG1655 assembly with 2D Oxford Nanopore long reads
echo Table 1F in manuscript 
echo ---------------------------------------------------------------------------------------------
echo *Must be running on a computer with at least 16 GB RAM, with current parameters
echo ---------------------------------------------------------------------------------------------
echo Downloading baseline assembly and ONT long reads... be patient.
wget ftp://ftp.bcgsc.ca/supplementary/LINKS/K12_full2dONT_longread.fa
wget ftp://ftp.bcgsc.ca/supplementary/LINKS/K12illumina-baseline.fa
echo ---------------------------------------------------------------------------------------------
echo Running LINKS iteratively a few times...be patient!
echo K12_full2dONT_longread.fa > K12_F2D.fof
/usr/bin/time -v -o timeLINKS_ECK12.txt ./runIterativeLINKS_ECK12.sh 15 2
