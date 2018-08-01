#!/usr/bin/env perl
#RLW 2010/2016

use strict;

if($#ARGV<1){
   print "Usage: $0\n";
   print "<fasta file 1>\n";
   print "<fasta file 2>\n";
   print "<read pair orientation 0/1,  0=raw MPET (<-->)   1=PET (-><-) >\n";
   die "** fasta files must have the same number of records & arranged in the same order\n";
}

my $file1=$ARGV[0];
my $file2=$ARGV[1];

open(IN1,$file1)||die "Can't open $file1 for reading -- fatal.\n";
open(IN2,$file2)||die "Can't open $file2 for reading -- fatal.\n";
my $po = $file1 . "_paired.fa";
my $upo = $file1 . "_unpaired.fa";
open(PAIR,">$po") || die "can't open $po for writing -- fatal.\n";
open(UNP,">$upo") || die "can't open $upo for writing -- fatal.\n";

my @all2 = <IN2>;

close IN2;
my $ct=0;
my $template="";
while(<IN1>){
   chomp;
   if(/^\@(\S+)/){
     $template=$1;
   }else{
     my $rev = $ct-1;
     chomp($all2[$rev]);
     if($all2[$rev]=~/$template/){
        chomp($all2[$ct]);

        my $rd1=$_;
        my $rd2=$all2[$ct];

        ###Will re-orient PET from -><- to <--> for LINKS, if input orientation is set
        $rd1=&reverseComplement($rd1) if($ARGV[2]);
        $rd2=&reverseComplement($rd2) if($ARGV[2]);

        if($rd1 ne "" && $rd2 ne ""){
           print PAIR ">$template\n$rd1:$rd2\n";
        }else{
           if($rd1 ne ""){
              print UNP ">$template";
              print UNP "1\n$rd1\n";
           }
           if($rd2 ne ""){
              print UNP ">$template";
              print UNP "2\n$rd2\n";
           }
        }
     }
   }
   $ct++;
}
close IN1;
close PAIR;
close UNP;

exit;

#-----------------------
sub reverseComplement{
   $_ = shift;
   $_ = uc();
   tr/ATGCYRKMBDHV/TACGRYMKVHDB/;
   return (reverse());
}
