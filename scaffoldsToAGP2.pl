#!/usr/bin/env perl
#Rene Warren
#V0.1
#3 Feb 2017
use strict;

#format that SSAKE/LINKS/ARCS outputs
#scaffold1,2884845,f2441408z328515k9a0.22m10_r3290314z119047k11a0.18m10_r2489402z151210k16a0.12m10_f2362377z142182k10a0.2m10_f2401086z197138k12a0.16m10_f2345317z259078k19a0.10m10_f2521875Z1508258k19a0.10m10_f2462456z179417
#
#
if($#ARGV<0){
   print "Script to convert LINKS .scaffolds output into AGPv2\n";
   die "Usage: $0 <.scaffolds> <arcs -- use this if used the ARCS pipeline>\n";
}

my $type;
my $gap;

open(IN,$ARGV[0]) || die "Can't open $ARGV[0] -- fatal.\n";
while(<IN>){
   chomp;
   my @a=split(/\,/);
   my @b=split(/\_/,$a[2]);
   my $bigsize=0;
#   foreach my $component(@b){
#      if($component=~/^([fr])(\d+)z(\d+)k\d+a\d+\.\d+m(.*)/i){
#         my ($ori,$id,$sz,$gap)=($1,$2,$3,$4);

#         if($gap<0){
#            $gap=100;
#         }
#	 $bigsize += ($sz + $gap);
#      }elsif($component=~/^([fr])(\d+)z(\d+)$/i){
#         my ($ori,$id,$sz,$gap)=($1,$2,$3);
#         $bigsize += ($sz);
#      }
#   }
   my $ct=0;
   my ($begin,$end,$size)=(1,0,0);

   foreach my $component(@b){
      my ($ori,$id,$type,$sz,$gap) = ("NA","NA","N",0,"NA");
      $ct++;
      if($component=~/^([fr])(\d+)z(\d+)k\d+a\d+\.\d+m(.*)/i){
         ($ori,$id,$sz,$gap)=($1,$2,$3,$4);
         if($ori eq "f"){
            $ori="+";
         }else{
            $ori="-";
         }
         if($gap<0){
            $gap=100;
            $type="U";
         }
      }elsif($component=~/^([fr])(\d+)z(\d+)$/i){
         ($ori,$id,$sz)=($1,$2,$3);
         if($ori eq "f"){
            $ori="+";
         }else{
            $ori="-";
         }
      }
      $end = $begin + $sz;
      if($ARGV[1] eq "arcs"){
         print "$a[0]\t$begin\t$end\t$ct\tU\t$id\t100\tNA\tscaffold\tNA\tyes\t$ori\tpaired-ends\n";
         $begin = $end + 100;
      }else{
	 print "$a[0]\t$begin\t$end\t$ct\t$type\t$id\t$gap\tNA\tscaffold\tNA\tyes\t$ori\tpaired-ends\n";
         $begin = $end + $gap;
      }
   }
}
close IN;

exit;
