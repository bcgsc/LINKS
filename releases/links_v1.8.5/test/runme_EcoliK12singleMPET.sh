d=`date`
echo ---------------------------------------------------------------------------------------------
echo Running LINKS on ABySS E.coli assembly with MPET data 
echo Downloading data...
echo --------------------------------------------------------------------------------------------- 
wget http://www.bcgsc.ca/downloads/supplementary/LINKS/EcMG1_S7_L001_R1_001.fastq.gz
wget http://www.bcgsc.ca/downloads/supplementary/LINKS/EcMG1_S7_L001_R2_001.fastq.gz
wget http://www.bcgsc.ca/downloads/supplementary/LINKS/set1_R1.mp.fastq.gz
wget http://www.bcgsc.ca/downloads/supplementary/LINKS/set1_R2.mp.fastq.gz
wget http://www.bcgsc.ca/downloads/supplementary/LINKS/ecoliK12_abyss_illumina_scaffold_baseline.fa
echo ---------------------------------------------------------------------------------------------
echo done. Prepping files...
echo ---------------------------------------------------------------------------------------------
echo Making fasta from MPET fastq
gunzip -c EcMG1_S7_L001_R1_001.fastq.gz | perl -ne '$ct++;if($ct>4){$ct=1;}print if($ct<3);' > mpet4k_1.fa
gunzip -c EcMG1_S7_L001_R2_001.fastq.gz | perl -ne '$ct++;if($ct>4){$ct=1;}print if($ct<3);' > mpet4k_2.fa
gunzip -c set1_R1.mp.fastq.gz | perl -ne '$ct++;if($ct>4){$ct=1;}print if($ct<3);' > trimmedmpet4k_1.fa
gunzip -c set1_R2.mp.fastq.gz | perl -ne '$ct++;if($ct>4){$ct=1;}print if($ct<3);' > trimmedmpet4k_2.fa
echo Generate the paired input:
../tools/makeMPETOutput2EQUALfiles.pl mpet4k_1.fa mpet4k_2.fa
###the parameter 1 in the command below instruct the script to re-orient the -><- input to <-->
../tools/makeMPETOutput2EQUALfiles.pl trimmedmpet4k_1.fa trimmedmpet4k_2.fa 1
echo mpet4k_1.fa_paired.fa > mpet.fof
echo trimmedmpet4k_1.fa_paired.fa > trimmedmpet.fof
echo ---------------------------------------------------------------------------------------------
echo done. Running LINKS sequentially on trimmed MPET and raw MPET..
echo ---------------------------------------------------------------------------------------------
/usr/bin/time -v -o timeLINKS_ECK12singleSCAFFtrimmedmpet.txt ../LINKS -k 15 -d 4000 -f ecoliK12_abyss_illumina_scaffold_baseline.fa -s trimmedmpet.fof -b ecoliK12-ONT_linksSingleIterationSCAFF-trimmedMPET -m 1
/usr/bin/time -v -o timeLINKS_ECK12singleSCAFFmpet.txt ../LINKS -k 15 -d 4000 -f ecoliK12_abyss_illumina_scaffold_baseline.fa -s mpet.fof -b ecoliK12-ONT_linksSingleIterationSCAFF-MPET -m 1
echo ---------------------------------------------------------------------------------------------
