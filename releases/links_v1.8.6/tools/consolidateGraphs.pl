#!/usr/bin/env perl
#Rene Warren 2018

use strict;

if($#ARGV<1){
	die "Usage: $0 <first.gv> <second.gv>\n";
}

my $g1;
my $s1;
my $s2;
my $s3;
my @g2;
my $merged = $ARGV[0] . "-" . $ARGV[1];


open(G1,$ARGV[0]);
while(<G1>){
   chomp;
   $g1->{$_}=1;
   $s1->{$1} = $2 if(/(.*)(\[.*\])/);
}

close G1;

my $cto=0;
my $ct = 0;
open(G2,$ARGV[1]);
open(OUT,">$merged");

while(<G2>){
   chomp;
   push @g2,$_;
   $s2->{$1} = $2 if(/(.*)(\[.*\])/);
   my ($left,$right) = ($1,$2) if(/(.*)(\[.*\])/);
   my ($node1,$node2) = ($1,$2) if($left=~/(\S+\s+)\-\-\s+(\S+\s+)/);
   $node1 = "\t" . $node1;
   $node2 = "\t" . $node2;
   if($g1->{$_}==1){##common
      $g2[$ct]=~s/deepskyblue/orchid/g;### will only replace instances that have blue
      $s2->{$1} = $2 if($g2[$ct]=~/(.*)(\[.*\])/);
      #print "$_ >>> $g2[$ct]\n";
   }else{###lines differ
      if($g2[$ct]=~/blue/){##lines differ, and second graph states blue, change for red
         $g2[$ct]=~s/deepskyblue/tomato/g;
      }else{### second graph has no highlights
         $cto++;
         if(defined $s1->{$left} && $s1->{$left} ne $right){
            if(defined $s1->{$node1} && $s1->{$node1} ne $s2->{$node1}){###cases where first graph has node highlights, not second

               if(defined $s2->{$node1}){
                  print OUT $node1 . $s2->{$node1} . "\n";
               }else{
                  print OUT $node1 . $s1->{$node1} . "\n";
               }
            }
            if(defined $s1->{$node2} && $s1->{$node2} ne $s2->{$node2}){###cases where first graph has node highlights, not second
               if(defined $s2->{$node2}){
                  print OUT $node2 . $s2->{$node2} . "\n";
               }else{
                  print OUT $node2 . $s1->{$node2} . "\n";
               }
            }
            $g2[$ct] = $left . $s1->{$left};

         #}else{
         #   print "$_\n";
         }
      }
   }
   print OUT "$g2[$ct]\n";
   $ct++;
}

close G2;

close OUT;
exit;
