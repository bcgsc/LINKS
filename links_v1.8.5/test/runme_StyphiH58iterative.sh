d=`date`
echo ---------------------------------------------------------------------------------------------
echo Running LINKS on baseline E. coli K12 MG1655 assembly with 2D Oxford Nanopore long reads
echo Figure 1 in manuscript - black diamonds
echo ---------------------------------------------------------------------------------------------
echo *Must be running on a computer with at least 16 GB RAM, with current parameters
echo ---------------------------------------------------------------------------------------------
echo Downloading baseline assembly and ONT long reads... be patient.
wget ftp://ftp.bcgsc.ca/supplementary/LINKS/StyphiH58_2dONT_longread.fa
wget ftp://ftp.bcgsc.ca/supplementary/LINKS/StyphiH58illumina-baseline.fa
echo ---------------------------------------------------------------------------------------------
echo StyphiH58_2dONT_longread.fa > STH58_2D.fof
echo Running LINKS iteratively a few times...be patient!
/usr/bin/time -v -o timeLINKS_STH58.txt ./runIterativeLINKS_STH58.sh 15 1 5 0.1
