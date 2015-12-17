# LINKS
Long Interval Nucleotide K-mer Scaffolder
LINKS v1.5.2 Rene L. Warren, 2014-2015
email: rwarren at bcgsc.ca

Name
----

LINKS: Long Interval Nucleotide K-mer Scaffolder


Description
-----------

LINKS is a genomics application for scaffolding genome assemblies with imperfect, long reads, such as those produced by Oxford Nanopore Technologies Ltd.

What's new in v1.5.2 ?
---------------------

LINKS outputs a scaffold graph in gv format, with highlighted merges and edge attributes


What's new in v1.5.1 ?
---------------------

Fixed a bug that prevented the creation of Bloom filters with a different false positive rate (FPR) than default.
Using lower FPR does not influence scaffolding itself, only run time.
For large genomes (>1Gbp), using a higher FPR is recommended when compute memory (RAM) is limiting.


What's new in v1.5 ?
---------------------

LINKS uses a Bloom filter to limit hashed paired k-mers to only those found in the sequence file to re-scaffold.
This feature decreases RAM usage by over 60%, while the run time is nearly unchanged.
When ran iteratively, users can re-use Bloom filters with the -r options, which results in faster run times up to half compared to v1.3 and earlier.


What's new in v1.3 ?
---------------------

Added support for fastq files.
Added support for multiple long-reads files. With v1.3, the reads file is not supplied directly through -s, but with a file-of-filenameinstead, which is a text file listing the fullpath/FASTA or FASTQ on your system. The file-of-filenames supplied through the -s option could include a mixture of FASTA and FASTQ files.


What's new in v1.2 ?
---------------------

Fixed bug that prevented reading traditional FASTA sequences (where a sequence is represented as a series of linestypically no longer than 120 characters)


What's new in v1.1 ?
---------------------

Included offset option (-o option) - Enable LINKS to explore a wider k-mer space range when running iteratively
Minor fixes: IUPAC codes are now preserved


Implementation and requirements
-------------------------------

LINKS is implemented in PERL and runs on any OS where PERL is installed.
In v1.5, there is a single dependency to Bloom::Faster - an extension for the c library libbloom.


Install
-------

Download the tar ball, gunzip and extract the files on your system using:

gunzip links_v1-5-2.tar.gz
tar -xvf links_v1-5-2.tar

The Bloom::Faster PERL library was built against various version of PERL

./lib/Bloom-Faster-1.7/
bloom5-10-0
bloom5-16-0
bloom5-16-3
bloom5-18-1

To re-build against YOUR PERL version, follow these easy steps:

1. cd ./lib/Bloom-Faster-1.7
2. /gsc/software/linux-x86_64/perl-5.10.0/bin/perl Makefile.PL PREFIX=./bloom5-10-0
(perl at YOUR favourite location)
3. make
4. make install
5. change line 30 in LINKS to reflect the location of the Bloom library
5. test: /gsc/software/linux-x86_64/perl-5.10.0/bin/perl ../../LINKS

Or change the shebang line of LINKS to point to the version of perl installed on your system and you're good to go.


Documentation
-------------

Refer to the LINKS.readme file on how to run LINKS and the LINKS web site for information about the software and its performance 
www.bcgsc.ca/bioinfo/software/links

Questions or comments?  We would love to hear from you!
rwarren at bcgsc.ca


Citing LINKS
------------

Thank you for using, developing and promoting this free software.

René L. Warren*, Benjamin P. Vandervalk, Steven JM Jones and Inanç Birol
LINKS: Scaffolding genome assemblies with kilobase-long nanopore reads
doi: http://dx.doi.org/10.1101/016519
*rwarren@bcgsc.ca
http://biorxiv.org/content/early/2015/03/20/016519


Running LINKS
-------------

e.g. /usr/bin/time -v -o timeLINKS_ECK12singleTIG.txt ../LINKS -f ecoliK12_abyss_illumina_contig_baseline.fa -s K12_F2D.fof -b ecoliK12-ONT_linksSingleIterationTIG

Usage: ../LINKS [v1.5.2]

-f  sequences to scaffold (Multi-FASTA format, required)

-s  file-of-filenames, full path to long sequence reads (Multi-FASTA/fastq format, required)

-d  distance between k-mer pairs (eg. target distance to re-scaffold on. default -d 4000, optional)

-k  k-mer value (default -k 15, optional)

-t  step of sliding window when extracting k-mer pairs from long reads (default -t 2, optional)

-o  offset position for extracting k-mer pairs (default -o 0, optional)

-e  error (%) allowed on -d distance   e.g. -e 0.1  == distance +/- 10% (default -e 0.1, optional)

-l  minimum number of links (k-mer pairs) to compute scaffold (default -l 5, optional)

-a  maximum link ratio between two best contig pairs (default -a 0.3, optional)
    *higher values lead to least accurate scaffolding*

-b  base name for your output files (optional)

-r  Bloom filter input file for sequences supplied in -s (optional, if none provided will output to .bloom)
    NOTE: BLOOM FILTER MUST BE DERIVED FROM THE SAME FILE SUPPLIED IN -f WITH SAME -k VALUE
    IF YOU DO NOT SUPPLY A BLOOM FILTER, ONE WILL BE CREATED (.bloom)

-p  Bloom filter false positive rate (default -p 0.0001, optional; increase to prevent memory allocation errors)

-x  Turn off Bloom filter functionality (-x 1 = yes, default = no, optional)

-v  Runs in verbose mode (-v 1 = yes, default = no, optional)

Notes:

-s K12_F2D.fof specifies a file-of-filenames (text file) listing: K12_full2dONT_longread.fa (see ./test folder)

-x When turned on (-x 1), LINKS will run with a behaviour equivalent to v1.3 (no Bloom filters).  

This may be useful for large genome assembly drafts and when long reads are extremely high quality.


Tips to minimize memory usage
-----------------------------

The most important parameters for decreasing RAM usage are -t and -d.
The largest dataset used for scaffolding by our group was a draft assembly of the white spruce genome (20 Gb)* - For this, a large sliding window, -t (200) was used and was decreased as the k-mer distance -d increased. 
*refer to LINKSrecipe_pglaucaPG29-WS77111.sh in the ./test folder

Because you want want to start with a low -d for scaffolding, you have to estimate how many minimum links (-l) would fit in a -d window +/- error -e given sliding window -t.
For instance, it may not make sense to use -t 200, -d 500 at low coverages BUT if you have at least 10-fold coverage it might since, in principle, you should be able to derive sufficient k-mer pairs within same locus if there's no bias in genome sequencing.

For re-scaffolding white spruce, only 1X coverage was available (since the re-scaffolding used a draft assembly instead of long reads), but even -t 200 -d 5000 (1st iteration) did merge scaffolds even though, in theory, the -e parameter will play an important role limiting linkages outside of the target range -d (+/-) -e %. 

On the data side of things, reducing the coverage (using less long reads), and limiting to only the highest quality reads would help decrease RAM usage.

In v1.5, LINKS builds a Bloom filter that comprises all k-mer of a supplied (-f) genome draft and uses it to only hash k-mer pairs from longreads having an equivalent in the Bloom filter. When LINKS runs iteratively, the bloom filter built at the first iteration is re-used thus saving execution time.

These are the best tips I can offer at the moment, until we address it further programatically using even more efficient data structures & code.


Test data
---------

Go to ./test

run:
-------------------------------------
./runme_EcoliK12single.sh

-rwxr-xr-x 1 rwarren users 1.2K Feb 26 07:41 runme_ScerevisiaeW303iterative.sh

The script will download the baseline E. coli abyss scaffold assembly and full 2D ONT reads (Quick et al 2014) and used the latter to re-scaffold the former, with default parameters (Table 1D in paper). 
NEED ~8GB RAM WITH CURRENT PARAMETERS. Increase (-t) to use less RAM.

-------------------------------------
./runme_EcoliK12iterative.sh

The script will download the baseline E. coli abyss scaffold assembly and full 2D ONT reads (Quick et al 2014) and used the latter to re-scaffold the former, iteratively 30 times increasing the distance between k-mer pairs at each iteration (Table 1F in paper). 
NEED ~16GB RAM WITH CURRENT PARAMETERS. Increase (-t) to use less RAM.

-------------------------------------
./runme_ScerevisiaeW303iterative.sh

This script will download the S. cerevisiae W303 raw ONT long reads and used them iteratively to scaffold a baseline IlluminaMiSeq assembly (both data from http://schatzlab.cshl.edu/data/nanocorr/).  You will need a computer with at least 132GB RAM.  This process was clocked at 6:08:21 (h:mm:ss wall clock) and used 118GB RAM on a Intel(R) Xeon(R) CPU E5-2699 v3 @ 2.30GHz 16 dualcore (but running on a single CPU). (Fig 1 in the main ms, FigS8 in preprint). 
NEED <132GB RAM WITH CURRENT PARAMETERS. Increase (-t) to use less RAM.

-------------------------------------
./runme_ScerevisiaeS288citerative.sh

This script will download the S. cerevisiae W303 raw ONT long reads and used
them iteratively to scaffold a baseline ABySS assembly of Illumina data.  You will need a computer with at least 132GB RAM.
(Fig 1 in the main ms, FigS8 in preprint).
NEED <132GB RAM WITH CURRENT PARAMETERS. Increase (-t) to use less RAM.

-------------------------------------
./runme_StyphiH58iterative.sh

This script will download the S. typhi H58 2D ONT long reads and used
them iteratively to scaffold a baseline assembly of Illumina data (both from
Ashton,P.M. 2015. Nat.Biotechnol.33,296–300). You will need a computer with at least 132GB RAM.
(Fig 1 in the main ms, FigS8 in preprint).

-------------------------------------
or 
./runall.sh

This script will run ALL of the above examples.


Additional info:

The file:

LINKSrecipe_pglaucaPG29-WS77111.sh

is provided to show the re-scaffolding recipe used to produce a re-scaffolded white spruce (P. glauca) genome assembly.

The file:

LINKSrecipe_athaliana_ectools.sh

and

LINKSrecipe_athaliana_raw.sh

are provided to show the re-scaffolding of the A. thaliana high-quality genome draft using ECTools-corrected or raw Pacific Biosciences reads.


How it works
------------

Process: nanopore reads are supplied as input (-s option, fasta format) and k-mer pairs are extracted using user-defined k-mer length (-k) and distance between the 5’-end of each pairs (-d) over a sliding window (-t). Unique k-mer pairs at set distance are hashed. Fasta sequences to scaffold are sup-plied as input (-f), and are shredded to k-mers on both strands, tracking the [contig] sequence of origin, k-mer positions and frequencies of observation. 

Algorithm: LINKS has two main stages, a contig pairing and a scaffold layout phase. Cycling through each k-mer pair, k-mers that are uniquely placed on contigs are identified, and putative contig pairs are formed if k-mers are placed on different contigs. Contig pairs are only considered if the calculated distances between them satisfy the mean distance provided (-d) while allowing for a deviation (-e). Contig pairs having a valid gap or over-lap are allowed to proceed to the scaffolding stage. Contigs in pairs may be ambiguous: a given contig may link to multiple contigs. To mitigate, the number of spanning k-mer pairs (links) between any given contig pair is recorded, along with a mean putative gap or overlap. Once pairing between contigs is complete, the scaffolds are built using contigs as seeds. Contigs are used in turn until all have been incorporated into a scaffold. Scaffolding is controlled by merging sequences only when a minimum number of links (-l) join two contig pairs, and when links are dominant compared to that of another possible pairing (-a). The predecessor of LINKS is the unpublished scaffolding engine in the widely-used SSAKE assembler (Warren et al. 2007), and foundation of the SSPACE-LongRead scaffolder (Boetzer and Pirovano, 2014).

Output: A summary of the scaffold layout is provided (.scaffold) as a text file and captures the linking information of successful scaffolds. A fasta file (.scaffold.fa) is generated using that information, placing sized N-pads for gaps and a single “n” in cases of overlaps between contigs. A log summary of k-mer pairing in the assembly is provided (.log) along with a text file describing possible issues in pairing (.pairing_issues) and pairing distribu-tion (.pairing_distribution.csv).


Consider the following contig pairs (AB, AC and rAD):

    A         B
========= ======== 
  ->       <-
   ->        <-
    ->      <-
       ->       <-

    A       C
========= ======
  ->        <-
    ->        <-

   rA        D           equivalent to rDA, in this order
========= =======
      ->   <-
     ->   <-
       ->   <-

Two parameters control scaffolding (-l and -a).  The -l option specifies the minimum number of links (read pairs) a valid contig pair MUST have to be considered.  The -a option specifies the maximum ratio between the best two contig pairs for a given seed/contig being extended.  For example, contig A shares 4 links with B and 2 links with C, in this orientation.  contig rA (reverse) also shares 3 links with D.   When it's time to extend contig A (with the options -l and -a set to 2 and 0.7, respectively), both contig pairs AB and AC are considered.  Since C (second-best) has 2 links and B (best) has 4 (2/4) = 0.5 below the maximum ratio of 0.7, A will be linked with B in the scaffold and C will be kept for another extension. If AC had 3 links the resulting ratio (0.75), above the user-defined maximum 0.7 would have caused the extension to terminate at A, with both B and C considered for a different scaffold.  A maximum links ratio of 1 (not recommended) means that the best two candidate contig pairs have the same number of links -- LINKS will accept the first one since both have a valid gap/overlap.  When a scaffold extension is terminated on one side, the scaffold is extended on the "left", by looking for contig pairs that involve the reverse of the seed (in this example, rAD).  With AB and AC having 4 and 2 links, respectively and rAD being the only pair on the left, the final scaffolds outputted by LINKS would be:

1) rD-A-B
2) C 

LINKS outputs a .scaffolds file with linkage information between contigs (see "Understanding the .scaffolds csv file" below)
Accurate scaffolding depends on many factors.  Number and nature of repeats in your target sequence, optimum adjustments of distance (-d), deviation on the distance (-e), kmer sizes (-k), Minimum number of links (-l) and link ratio (-a) and data quality will all affect LINKS's ability to build scaffolds.

NOTE: IT IS ADVISED TO RUN LINKS WITH SMALLER DISTANCES (-d) FIRST, ESPECIALLY WHEN ASSEMBLIES ARE VERY FRAGMENTED.

.log                      :: text file; Logs execution time / errors / pairing stats
.pairing_distribution.csv :: comma-separated file; 1st column is the calculated distance for each pair (template) with reads that assembled logically within the same contig.  2nd column is the number of pairs at that distance
.pairing_issues           :: text file; Lists all pairing issues encountered between contig pairs and illogical/out-of-bounds pairing
.scaffolds                :: comma-separated file; see below
.scaffolds.fa             :: fasta file of the new scaffold sequence
.bloom                    :: Bloom filter created by shredding the -f input into k-mers of size -k


Understanding the .scaffolds csv file
-------------------------------------

scaffold1,7484,f127Z7068k12a0.58m42_f3090z62k7a0.14m76_f1473z354

column 1: a unique scaffold identifier
column 2: the sum of all contig sizes that made it to the scaffold/supercontig
column 3: a contig chain representing the layout:

e.g.
f127Z7068k12a0.58m42_f3090z62k7a0.14m76_f1473z354

means: contig f127 (strand=f/+), size (z) 7068 (Z if contig was used as the seed sequence) has 12 links (k), link ratio of 0.58 (a) with a mean gap of 42nt (m) with reverse (r) of contig 3090 (size 62) on the right.  if m values are negative, it's just that a possible overlap was calculated using the mean distance supplied by the user and the position of the reads flanking the contig.
Negative m values imply that there's a possible overlap between the contigs.  But since the pairing distance distribution usually follows a Normal/Gaussian distribution, some distances are expected to be larger than the median size expected/observed.  In reality, if the exact size was known between each paired-reads, we wouldn't expect much negative m values unless a break occurred during the contig extension (likely due to base errors/SNPs). 


License
-------

LINKS Copyright (c) 2014-2015 Canada's Michael Smith Genome Science Centre.  All rights reserved.

LINKS is released under the GNU General Public License - GPL

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
