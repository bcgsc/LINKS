![Logo](https://github.com/warrenlr/links/blob/master/links-logo.png)

# LINKS

## Long Interval Nucleotide K-mer Scaffolder
## LINKS v1.8.7 Rene L. Warren, 2014-2019
## email: rwarren [at] bcgsc [dot] ca
## Visit www.bcgsc.ca/bioinfo/software/links for additional information


### Description
-----------

LINKS is a genomics application for scaffolding genome assemblies with long
reads, such as those produced by Oxford Nanopore Technologies Ltd.
It can be used to scaffold high-quality draft genome assemblies with any long
sequences (eg. ONT reads, PacBio reads, another draft genomes, etc).
It is also used to scaffold contig pairs linked by ARCS/ARKS.

### What's new in v1.8.7 ?
---------------------

minor bug fixes, improved documentation


### What's new in v1.8.6 ?
---------------------

When pipelined with the ARCS/ARKS scaffolder, LINKS v1.8.6 prioritizes paths with shorter gaps only when there are no ambiguous gaps distances with the neighboring sequences considered. Add linkage information (ID and size) to the output .gv graph. We now provide a utility (./tools/consolidateGraphs.pl) to highlight differences between two LINKS graphs. Implement -z min_size with an ARCS checkpoint.


### What's new in v1.8.5 ?
---------------------

When pipelined with the ARCS scaffolder, LINKS now extracts information from the tigpair_checkpoint.tsv file to generate a enhanced .gv file with additional information regarding other potential linking partners, number of supports and sequence orientation.

**WHEN USING LONG READS FOR SCAFFOLDING, WE RECOMMEND THE USE OF v1.8.5 (FOR NOW), AS IT CONSUMES LESS MEMORY. THE LATEST RELEASE CONCERNS THE ARC/KS PIPELINE.**


### What's new in v1.8.4 ?
---------------------

Changed license to GPLv3 


### What's new in v1.8.3 ?
---------------------

Fixes a bug introduced in v1.8.1 that caused sequence overuse


### What's new in v1.8.2 ?
---------------------

Implements the -z option, minimum contig size cutoff for scaffolding


### What's new in v1.8.1 ?
---------------------

Stratifies/prioritizes short-to-long distances when building the scaffold layout 


### What's new in v1.8 ?
---------------------

Native support for iterative k-mer pair extraction at distinct length
intervals


### What's new in v1.7 ?
---------------------

Support for scaffolding with MPET (jumping library) reads
Support for reading compressed long sequence [reads] and assembly files
Implemented mid-scaffolding checkpoint to:

- more quickly test certain parameters (-l min. links / -a min. links ratio)
- quickly recover from crash
- explore very large kmer spaces


### What's new in v1.6.1 ?
---------------------

Added a new output file, .assembly_correspondence.tsv 
This human-readable correspondence file lists the scaffold ID, contig ID, original assembly contig name, orientation, #linking kmer pairs, links ratio, gap or overlap


### What's new in v1.6 ?
---------------------

Incorporation of the BC Genome Sciences Centre custom Bloom filter with
the fast ntHash recursive nucleotide hash function.  This new data structure supports the creation of
Bloom filters from large genome assemblies (tested on assemblies of 3 Gbp human and 20 Gbp white spruce).

The Bloom filter data structure swap in v1.6+ offers a ~30-fold kmer insert
speed-up (~6x query speed-up) over v1.5.2, while supporting the
creation of filters from large genome assembly drafts.


### What's new in v1.5.2 ?
---------------------

LINKS outputs a scaffold graph in gv format, with highlighted merges and edge attributes  


### What's new in v1.5.1 ?
---------------------

Fixed a bug that prevented the creation of Bloom filters with a different false positive rate (FPR) than default.
Using lower FPR does not influence scaffolding itself, only run time. 
For large genomes (>1Gbp), using a higher FPR is recommended when compute memory (RAM) is limiting.


### What's new in v1.5 ?
---------------------

LINKS uses a Bloom filter to limit hashed paired k-mers to only those found in the sequence file to re-scaffold.
This feature decreases RAM usage by over 60%, while the run time is nearly unchanged.
When ran iteratively, users can re-use Bloom filters with the -r options, which results in faster run times up to half compared to v1.3 and earlier.


### What's new in v1.3 ?
---------------------

Added support for fastq files.
Added support for multiple long-reads files. With v1.3, the reads file is not supplied directly through -s, but with a file-of-filenameinstead, which is a text file listing the fullpath/FASTA or FASTQ on your system. The file-of-filenames supplied through the -s option could include a mixture of FASTA and FASTQ files.


### What's new in v1.2 ?
---------------------

Fixed bug that prevented reading traditional FASTA sequences (where a sequence is represented as a series of lines typically no longer than 120 characters)


### What's new in v1.1 ?
---------------------

Included offset option (-o option) - Enables LINKS to explore a wider k-mer space range when running iteratively
Minor fixes: IUPAC codes are now preserved


### Implementation and requirements
-------------------------------

LINKS is implemented in PERL and runs on any OS where PERL is installed.
In v1.6, there is a single dependency to the BloomFilter.pm (included) - BC Genome Sciences Centre's common Bloom filter
In v1.5, there is a single dependency to Bloom::Faster - an extension for the c library libbloom.


### Install
-------

Download the tar ball, gunzip and extract the files on your system using:
<pre>
gunzip links_v1-8-7.tar.gz
tar -xvf links_v1-8-7.tar
</pre>
In v1.6 and higher, the use of the Bloom::Faster PERL library is deprecated

### Instructions for building the BloomFilter PERL module
-------

1. DOWNLOAD the BC Genome Sciences Centre's BloomFilter: The BTL C/C++ Common
Bloom filters for bioinformatics projects, as well as any APIs created for
other programming languages.
<pre>
cd ./links_v1.8.7/lib
git clone git://github.com/bcgsc/bloomfilter.git
cd swig
</pre>

2. BUILD a PERL5 module

Make sure you have swig installed and included in your path.

http://www.swig.org/


TO BUILD a Perl5 module (run in swig/):
```
a) /home/<user>/<path to swig>/preinst-swig -Wall -c++ -perl5 BloomFilter.i
b) g++ -c BloomFilter_wrap.cxx -I/usr/lib64/perl5/CORE -fPIC -Dbool=char -O3
c) g++ -Wall -shared BloomFilter_wrap.o -o BloomFilter.so -O3
```

TO COMPILE, swig needs the following Perl5 headers:
```C++
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
```
If they are not located in /usr/lib64/perl5, you can run "perl -e 'use Config; print $Config{archlib};" to locate them.


3. VERIFY your install

[swig]$ ./test.pl


4. CHANGE the path to BloomFilter.pm in LINKS/writeBloom.pl/testBloom.pl

You only need to change if you have re-built in a relative directory different
from:
<pre>
use lib "$FindBin::Bin/./lib/bloomfilter/swig"; (for LINKS)
use lib "$FindBin::Bin/../lib/bloomfilter/swig"; (for writeBloom.pl/testBloom.pl)
</pre>


### Documentation
-------------

Refer to the LINKS-readme.txt/LINKS-readme.pdf file on how to run LINKS and the LINKS web site for information about the software and its performance 
www.bcgsc.ca/bioinfo/software/links

Questions or comments?  We would love to hear from you!
rwarren at bcgsc.ca


### Citing LINKS
------------

<pre>
Rene L. Warren, Chen Yang, Benjamin P. Vandervalk, Bahar Behsaz,
Albert Lagman, Steven J. M. Jones and Inanç Birol. 2015. LINKS: Scalable,
alignment-free scaffolding of draft genomes with long reads. GigaScience 4:35
DOI: 10.1186/s13742-015-0076-3©  Warren et al. 2015
</pre>

Thank you for using, developing and promoting this free software.


### Credits
-------

LINKS:
Rene Warren

SWIG/BloomFilter.pm:
Sarah Yeo
Justin Chu

https://github.com/bcgsc/bloomfilter:
Justin Chu
Ben Vandervalk
Hamid Mohamadi (recursive/ntHash)
Sarah Yeo
Golnaz Jahesh


### Running LINKS
-------------
<pre>
e.g. ./LINKS -f ecoliK12_abyss_illumina_contig_baseline.fa -s K12_F2D.fof -b ecoliK12-ONT_linksSingleIterationTIG

Usage: ./LINKS [v1.8.7]
-f  sequences to scaffold (Multi-FASTA format with each sequence on a single line, required)
-s  file-of-filenames, full path to long sequence reads or MPET pairs [see below] (Multi-FASTA/fastq format, required)
-m  MPET reads (default -m 1 = yes, default = no, optional
	! DO NOT SET IF NOT USING MPET. WHEN SET, LINKS WILL EXPECT A SPECIAL FORMAT UNDER -s
	! Paired MPET reads in their original outward orientation <- -> must be separated by ":"
	  >template_name
	  ACGACACTATGCATAAGCAGACGAGCAGCGACGCAGCACG:ATATATAGCGCACGACGCAGCACAGCAGCAGACGAC
-d  distance between k-mer pairs (ie. target distances to re-scaffold on. default -d 4000, optional)
	Multiple distances are separated by comma. eg. -d 500,1000,2000,3000
-k  k-mer value (default -k 15, optional)
-t  step of sliding window when extracting k-mer pairs from long reads
(default -t 2, optional)
	Multiple steps are separated by comma. eg. -t 10,5
-o  offset position for extracting k-mer pairs (default -o 0, optional)
-e  error (%) allowed on -d distance   e.g. -e 0.1  == distance +/- 10%
(default -e 0.1, optional)
-l  minimum number of links (k-mer pairs) to compute scaffold (default -l 5, optional) 
-a  maximum link ratio between two best contig pairs (default -a 0.3, optional)
	 *higher values lead to least accurate scaffolding*
-z  minimum contig length to consider for scaffolding (default -z 500, optional)
-b  base name for your output files (optional)
-r  Bloom filter input file for sequences supplied in -s (optional, if none provided will output to .bloom)
	 NOTE: BLOOM FILTER MUST BE DERIVED FROM THE SAME FILE SUPPLIED IN -f WITH SAME -k VALUE
	 IF YOU DO NOT SUPPLY A BLOOM FILTER, ONE WILL BE CREATED (.bloom)
-p  Bloom filter false positive rate (default -p 0.001, optional; increase to prevent memory allocation errors)
-x  Turn off Bloom filter functionality (-x 1 = yes, default = no, optional)
-v  Runs in verbose mode (-v 1 = yes, default = no, optional)


Notes:
-s K12_F2D.fof specifies a file-of-filenames (text file) listing: K12_full2dONT_longread.fa (see ./test folder)
   -f and -s : sequences must be on a SINGLE line with no linebreaks
   eg.  
   >LONGREAD-1
   AATACAATAGACGCACA...ATGAACGCAGACTTACAG
   >LONGREAD-2
   TGTGCTCTCTGTAATGTTC...ATACAGAACACGCAGCCAAGCGA

-x When turned on (-x 1), LINKS will run with a behaviour equivalent to v1.3 (no Bloom filters).  
This may be useful for large genome assembly drafts and when long reads are extremely high quality.
</pre>


### Tips to minimize memory usage and additional notes
--------------------------------------------------

The most important parameters for decreasing RAM usage are -t and -d.
The largest dataset used for scaffolding by our group was a draft assembly of the white spruce genome (20 Gb)* - For this, a large sliding window, -t (200) was used and was decreased as the k-mer distance -d increased. 
*refer to LINKSrecipe_pglaucaPG29-WS77111.sh in the ./test folder

Because you want want to start with a low -d for scaffolding, you have to estimate how many minimum links (-l) would fit in a -d window +/- error -e given sliding window -t.
For instance, it may not make sense to use -t 200, -d 500 at low coverages BUT if you have at least 10-fold coverage it might since, in principle, you should be able to derive sufficient k-mer pairs within same locus if there's no bias in genome sequencing.

For re-scaffolding white spruce, only 1X coverage was available (since the re-scaffolding used a draft assembly instead of long reads), but even -t 200 -d 5000 (1st iteration) did merge scaffolds even though, in theory, the -e parameter will play an important role limiting linkages outside of the target range -d (+/-) -e %. This is especially true when using raw MPET for scaffolding, to limit spurious linkages by contaminating PETs.

On the data side of things, reducing the coverage (using less long reads), and limiting to only the highest quality reads would help decrease RAM usage.

In v1.5, LINKS builds a Bloom filter that comprises all k-mer of a supplied (-f) genome draft and uses it to only hash k-mer pairs from longreads having an equivalent in the Bloom filter. When LINKS runs iteratively, the bloom filter built at the first iteration is re-used thus saving execution time.

In v1.8, Users may input multiple distances using the -d parameter, separating
each distance by a comma. 
eg. -d 500,1000,2000,4000,5000
will have for effect to extract kmer pairs at these five distance intervals.

Similarly, the window step size now accepts multiple integers, each separated
by a comma and with the order matching that in -d
However, the size of the array can be shorter, and the last valid -t will be
propagated to subsequent distances when they are not defined. 
eg. -t 20,10,5
A step size of 20, 10, 5, 5, 5 bp will be used when exploring the distances 
supplied above.

In v1.8, a single round of scaffolding is done, using the kmer pair space 
extracted at the specified distances.
Accordingly, v1.8 is not expected to yield the exact same results as
separate iterative LINKS runs. For users comfortable with the original set up,
no change is needed to make use of the previous LINKS behavior. 

Simultaneous exploration of vast kmer space is expected to yield better scaffolding results.

WARNING:
Specifying many distances will require large amount of RAM, especially with
low -t values.


### Test data
---------
<pre>
Go to ./test
</pre>
(cd test)


run:
<pre>
-------------------------------------
./runme_EcoliK12single.sh

The script will download the baseline E. coli abyss scaffold assembly and full 2D ONT reads (Quick et al 2014) and used the latter to re-scaffold the former, with default parameters (Table 1D in paper). 
NEED ~8GB RAM WITH CURRENT PARAMETERS. Increase (-t) to use less RAM.

./runme_EcoliK12singleMPET.sh will scaffold using E. coli K12 MPET reads
(~42-90GB RAM for trimmed vs raw MPET)

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
</pre>

Additional info:

The file: 
LINKSrecipe_pglaucaPG29-WS77111.sh
is provided to show the re-scaffolding recipe used to produce a re-scaffolded white spruce (P. glauca) genome assembly.

The file:
LINKSrecipe_athaliana_ectools.sh
and
LINKSrecipe_athaliana_raw.sh
are provided to show the re-scaffolding of the A. thaliana high-quality genome draft using ECTools-corrected or raw Pacific Biosciences reads.


### Testing the Bloom filters
-------------------------
<pre>
# To test insertions:
cd tools
./writeBloom.pl
Usage: ./writeBloom.pl
-f  sequences to scaffold (Multi-FASTA format, required)
-k  k-mer value (default -k 15, optional)
-p  Bloom filter false positive rate (default -p 0.0001, optional - increase to prevent memory allocation errors)

# To test queries:
cd tools
./testBloom.pl
Usage: ./testBloom.pl
-f  sequences to test (Multi-FASTA format, required)
-k  k-mer value (default -k 15, optional)
-r  Bloom filter file
</pre>

### How it works
------------

Process: nanopore/long reads are supplied as input (-s option, fasta/fastq format) and k-mer pairs are extracted using user-defined k-mer length (-k) and distance between the 5’-end of each pairs (-d) over a sliding window (-t). Unique k-mer pairs at set distance are hashed. Fasta sequences to scaffold are sup-plied as input (-f), and are shredded to k-mers on both strands, tracking the [contig] sequence of origin, k-mer positions and frequencies of observation. 

Algorithm: LINKS has two main stages, a contig pairing and a scaffold layout phase. Cycling through each k-mer pair, k-mers that are uniquely placed on contigs are identified, and putative contig pairs are formed if k-mers are placed on different contigs. Contig pairs are only considered if the calculated distances between them satisfy the mean distance provided (-d) while allowing for a deviation (-e). Contig pairs having a valid gap or over-lap are allowed to proceed to the scaffolding stage. Contigs in pairs may be ambiguous: a given contig may link to multiple contigs. To mitigate, the number of spanning k-mer pairs (links) between any given contig pair is recorded, along with a mean putative gap or overlap. Once pairing between contigs is complete, the scaffolds are built using contigs as seeds. Contigs are used in turn until all have been incorporated into a scaffold. Scaffolding is controlled by merging sequences only when a minimum number of links (-l) join two contig pairs, and when links are dominant compared to that of another possible pairing (-a). The predecessor of LINKS is the unpublished scaffolding engine in the widely-used SSAKE assembler (Warren et al. 2007), and foundation of the SSPACE-LongRead scaffolder (Boetzer and Pirovano, 2014).

Output: A summary of the scaffold layout is provided (.scaffold) as a text file and captures the linking information of successful scaffolds. A fasta file (.scaffold.fa) is generated using that information, placing sized N-pads for gaps and a single “n” in cases of overlaps between contigs. A log summary of k-mer pairing in the assembly is provided (.log) along with a text file describing possible issues in pairing (.pairing_issues) and pairing distribu-tion (.pairing_distribution.csv).


Consider the following contig pairs (AB, AC and rAD):
<pre>
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
</pre>

Two parameters control scaffolding (-l and -a).  The -l option specifies the minimum number of links (read pairs) a valid contig pair MUST have to be considered.  The -a option specifies the maximum ratio between the best two contig pairs for a given seed/contig being extended.  For example, contig A shares 4 links with B and 2 links with C, in this orientation.  contig rA (reverse) also shares 3 links with D.   When it's time to extend contig A (with the options -l and -a set to 2 and 0.7, respectively), both contig pairs AB and AC are considered.  Since C (second-best) has 2 links and B (best) has 4 (2/4) = 0.5 below the maximum ratio of 0.7, A will be linked with B in the scaffold and C will be kept for another extension. If AC had 3 links the resulting ratio (0.75), above the user-defined maximum 0.7 would have caused the extension to terminate at A, with both B and C considered for a different scaffold.  A maximum links ratio of 1 (not recommended) means that the best two candidate contig pairs have the same number of links -- LINKS will accept the first one since both have a valid gap/overlap.  When a scaffold extension is terminated on one side, the scaffold is extended on the "left", by looking for contig pairs that involve the reverse of the seed (in this example, rAD).  With AB and AC having 4 and 2 links, respectively and rAD being the only pair on the left, the final scaffolds outputted by LINKS would be:

1. rD-A-B
2. C 

LINKS outputs a .scaffolds file with linkage information between contigs (see "Understanding the .scaffolds csv file" below)
Accurate scaffolding depends on many factors.  Number and nature of repeats in your target sequence, optimum adjustments of distance (-d), deviation on the distance (-e), kmer sizes (-k), Minimum number of links (-l) and link ratio (-a) and data quality will all affect LINKS's ability to build scaffolds.

NOTE: IT IS ADVISED TO RUN LINKS WITH SMALLER DISTANCES (-d) FIRST, ESPECIALLY WHEN ASSEMBLIES ARE VERY FRAGMENTED.


### MPET INPUT
---------------------

In v1.7, a new option (-m) instructs LINKS that the long-read source (-s) is MPET. The users should prepare their input as specified in:
cd test
runme_EcoliK12singleMPET.sh

The MPET input is a custom format akin to FASTA and the sequence record must consist of read1:read2
<pre>
>template
ACGACACATCTACGCAGCGACGACGATAAATATAC:ATCAGCACAGCGACGCAGCGACAGCAGGACGACGAC
</pre>

NOTES:

- Paired MPET reads are supplied in their original outward orientation <- ->
- MPET sequences do not need to be trimmed (the Bloom filter will take care of eliminating erroneous kmers not found in the assembly)
- You CANNOT combine MPET and long reads simultaneously in the same LINKS process
- You may trim or process MPET reads if you wish (eg. with NxTrim), but remember to supply resulting MPETs in their original, outward-facing configuration (ie. <- ->). The script in ./tools/makeMPETOutput2EQUALfiles.pl does that for you.

The default behaviour is to extract kmer pairs from long-read FASTA/FASTQ files specified in -s.

Alternatively, when set to the MPET read length, the -m option will signal LINKS to
extract kmer pairs across a distance set in -d, for each MPET pair supplied in files supplied under -s

When doing so, ensure that -t is set to extract at least ~5 kmer pairs/MPET pair
As a rule of thumb, -l should be set to at least double that value (-l 10 in this case)


### Preparing the MPET input
------------------------
<pre>
For each fastq MPET file, convert in fasta:
 gunzip -c EcMG1_S7_L001_R1_001.fastq.gz | perl -ne '$ct++;if($ct>4){$ct=1;}print if($ct<3);' > mpet4k_1.fa
 gunzip -c EcMG1_S7_L001_R2_001.fastq.gz | perl -ne '$ct++;if($ct>4){$ct=1;}print if($ct<3);' > mpet4k_2.fa
 gunzip -c set1_R1.mp.fastq.gz | perl -ne '$ct++;if($ct>4){$ct=1;}print if($ct<3);' > trimmedmpet4k_1.fa
 gunzip -c set1_R2.mp.fastq.gz | perl -ne '$ct++;if($ct>4){$ct=1;}print if($ct<3);' > trimmedmpet4k_2.fa

Generate the paired input (refer to the tools folder):
Usage: ./makeMPETOutput2EQUALfiles.pl
<fasta file 1>
<fasta file 2>
<read pair orientation 0/1,  0=raw MPET (<-->)   1=PET (-><-) >
** fasta files must have the same number of records & arranged in the same order

RAW: ../tools/makeMPETOutput2EQUALfiles.pl mpet4k_1.fa mpet4k_2.fa       
TRIMMED: ../tools/makeMPETOutput2EQUALfiles.pl trimmedmpet4k_1.fa trimmedmpet4k_2.fa 1

echo mpet4k_1.fa_paired.fa > mpet.fof
echo trimmedmpet4k_1.fa_paired.fa > trimmedmpet.fof
</pre>

### OUTPUT FILES
------------------------

|Output files|                    Description|
|---|---|
|.log                         | text file; Logs execution time / errors / pairing stats|
|.pairing_distribution.csv    | comma-separated file; 1st column is the calculated distance for each pair (template) with reads that assembled logically within the same contig.  2nd column is the number of pairs at that distance|
|.pairing_issues              | text file; Lists all pairing issues encountered between contig pairs and illogical/out-of-bounds pairing |
|.scaffolds                   | comma-separated file; see below |
|.scaffolds.fa                | fasta file of the new scaffold sequence |
|.bloom                       | Bloom filter created by shredding the -f input into k-mers of size -k |
|.gv                          | scaffold graph (for visualizing merges), can be rendered in neato, graphviz, etc |
|.assembly_correspondence.tsv | correspondence file lists the scaffold ID, contig ID, original_name, #linking kmer pairs, links ratio, gap or overlap|
|.simplepair_checkpoint.tsv   | checkpoint file, contains info to rebuild datastructure for .gv graph |
|.tigpair_checkpoint.tsv      | if -b BASNAME.tigpair_checkpoint.tsv is present, LINKS will skip the kmer pair extraction and contig pairing stages. Delete this file to force LINKS to start at the beginning. This file can be used to: 1) quickly test parameters (-l min. links / -a min. links ratio, 2) quickly recover from crash 3) explore very large kmer spaces 4) scaffold with output of ARCS |


#### Interpreting .assembly_correspondence.tsv
-------------------------------------

This human-readable correspondence file lists the scaffold ID, contig ID,
original assembly contig name, contig/sequence orientation, #linking kmer pairs, links ratio,
gap or overlap(-) in this order


#### Interpreting the graph / .gv file
-------------------------------------

1. Vertices correspond to the sequences being considered for scaffolding, with the LINKS re-numbered sequences displayed in each vertex (unlinked sequences are not shown)

2. Edges are drawned between vertices when there is evidence for linking scaffolds (even if they are no ultimately scaffolded)

3. Only vertices/scaffolds highlighted in blue satisfied user-specified scaffold criteria (l and a parameters and satisfied logic/distance). These are scaffolded in the final LINKS output

4. Each edge in the graph will have 3 types of information (l=,g=,type=)

l=:number of kmer pairs linking any two vertices/sequences

g=:estimated gap or overlap (-) length between any two sequences

type=:refers to the orientation of the sequences (forward=1,reverse=0)


#### Understanding the .scaffolds csv file
-------------------------------------
<pre>
scaffold1,7484,f127Z7068k12a0.58m42_f3090z62k7a0.14m76_f1473z354
</pre>

column 1: a unique scaffold identifier
column 2: the sum of all contig sizes that made it to the scaffold/supercontig
column 3: a contig chain representing the layout:

e.g.
f127Z7068k12a0.58m42_f3090z62k7a0.14m76_f1473z354

means: contig f127 (strand=f/+), size (z) 7068 (Z if contig was used as the seed sequence) has 12 links (k), link ratio of 0.58 (a) with a mean gap of 42nt (m) with reverse (r) of contig 3090 (size 62) on the right.  if m values are negative, it's just that a possible overlap was calculated using the mean distance supplied by the user and the position of the reads flanking the contig.
Negative m values imply that there's a possible overlap between the contigs.  But since the pairing distance distribution usually follows a Normal/Gaussian distribution, some distances are expected to be larger than the median size expected/observed.  In reality, if the exact size was known between each paired-reads, we wouldn't expect much negative m values unless a break occurred during the contig extension (likely due to base errors/SNPs). 


### License
-------

LINKS Copyright (c) 2014-2019 British Columbia Cancer Agency Branch.  All rights reserved.

LINKS is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
 
For commercial licensing options, please contact
Patrick Rebstein <prebstein@bccancer.bc.ca>
