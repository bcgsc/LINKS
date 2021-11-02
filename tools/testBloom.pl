#!/usr/bin/env perl

#AUTHOR
#   Rene Warren
#   rwarren at bcgsc.ca

#NAME
#testBloom

#SYNOPSIS

#DOCUMENTATION
#   LINKS-readme.txt distributed with this software @ www.bcgsc.ca
#   http://www.bcgsc.ca/platform/bioinfo/software/links
#   We hope this code is useful to you -- Please send comments & suggestions to rwarren * bcgsc.ca
#   If you use LINKS, the LINKS code or ideas, please cite our work

#LICENSE
#   LINKS Copyright (c) 2014-2018 Canada's Michael Smith Genome Science Centre.  All rights reserved.

#LINKS is released under the BC Cancer Agency software license agreement (academic use). Details of the license can be accessed at: http://www.bcgsc.ca/platform/bioinfo/license/bcca_2010
#For commercial use, please contact rwarren@bcgsc.ca


use strict;
use POSIX;
use FindBin;
use lib "$FindBin::Bin/../lib/bloomfilter/swig";
use BloomFilter;
use Getopt::Std;
use Net::SMTP;
use vars qw($opt_f $opt_k $opt_r);
getopts('f:k:r:');
my ($k,$bf_file,)=(15,"");

#-------------------------------------------------

if(! $opt_f ){
   print "Usage: $0\n";
   print "-f  sequences to test (Multi-FASTA format, required)\n"; 
   print "-k  k-mer value (default -k $k, optional)\n";
   die "-r  Bloom filter file\n";
}

my $assemblyfile = $opt_f;
$k = $opt_k if($opt_k);
$bf_file = $opt_r if($opt_r); 

print "\nRunning:$0 -f $assemblyfile -k $k -b $bf_file\n\n";


if(! -e $assemblyfile){
   my $file_message = "\nInvalid file: $assemblyfile -- fatal\n";
   print $file_message;
   exit;
}elsif(! -e $bf_file){
   my $file_message = "\nInvalid file: $bf_file -- fatal\n";
   print $file_message;
   exit;
}else{
   my $file_message = "Checking sequence target file $assemblyfile...ok\n";
   print $file_message;
}



my $date = `date`;
chomp($date);

eval{

my $bloom = new BloomFilter::BloomFilter($bf_file);
print "Loading bloom filter from $bf_file\n";

my $date = `date`;
chomp($date);

print "$date:Shredding supplied sequence file (-f $assemblyfile) into $k-mers and testing against your Bloom filter..\n";
&contigsToBloom($assemblyfile,$k,$bloom);

my $date = `date`;
chomp($date);


};

if($@){
   my $message = $@;
   my $failure = "\nSomething went wrong running $0 $date\n$message\n";
   print $failure;
}else{
   my $success = "\n$date:$0 executed normally\n";
   print $success;
}

exit 1;


#----------------
sub contigsToBloom{
   my ($file,$k,$bloom) = @_;

   my $prevhead = "";
   my $seq = "";
   my $cttig=0;

   my ($ct_hit,$ct_total) = (0,0);
   open(IN,$file) || die "Error reading $file -- fatal.\n";

   print "Contigs processed k=$k:\n";
   ###
   while(<IN>){
       chomp;
      if(/^\>(\S+)/){
         my $head=$1;

         if ($head ne $prevhead && $seq ne '' && $prevhead ne ''){
            $cttig++;
            print "\r$cttig";
            $|++;
            ($ct_hit,$ct_total) = &kmerizeContigBloom(uc($seq),$bloom,$k,$ct_hit,$ct_total);
         }
         $seq = '';
         $prevhead = $head;
      }else{
         $seq .= $_;
      }
   }
   $cttig++;
   print "\r$cttig";
   $|++;
   ($ct_hit,$ct_total) = &kmerizeContigBloom(uc($seq),$bloom,$k,$ct_hit,$ct_total);
   ###
   close IN;

   print "\n\nFound $ct_hit out of $ct_total kmers probable in your Bloom filter\n";
}

#----------------
sub kmerizeContigBloom{
   my ($seq,$bloom,$k,$ct_hit,$ct_total) = @_;

   my $string = substr($seq, 0, 15);
   for(my $pos=0;$pos<=(length($seq)-$k);$pos++){
      $ct_total++;
      my $kmer = substr($seq,$pos,$k);
      if($bloom->contains($kmer)){
         $ct_hit++;
      }else {
          print "Missing: $ct_total - $kmer\n";
      }
   }
   return $ct_hit,$ct_total;
}
    

#-----------------------
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGCYRKMBDHV/TACGRYMKVHDB/; 
   return (reverse());
}


