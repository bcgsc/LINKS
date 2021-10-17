#!/usr/bin/env perl
#Rene Warren
#https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/
#V0.1
#3 Feb 2017
#mod 24apr2019
use strict;

#format SSAKE/LINKS/ARCS outputs:
#scaffold1,14000,f34z4000k8a0.3m3_r99Z10000
#
if($#ARGV<0){
   print "Script to convert LINKS .scaffolds output into AGPv2\n";
   die "Usage: $0 <.scaffolds>\n";
}

my $type;
my $gap;

open(IN,$ARGV[0]) || die "Can't open $ARGV[0] -- fatal.\n";
while(<IN>){
   chomp;
   my @a=split(/\,/);
   my @b=split(/\_/,$a[2]);
   my $bigsize=0;

   my $ct=0;
   my ($begin,$end,$size)=(1,0,0);
   my $el=0;
   foreach my $component(@b){
      my ($ori,$id,$type,$sz,$gap) = ("NA","NA","N",0,"NA");
      $ct++;
      if($component=~/^([fr])(\d+)z(\d+)k\d+a\S+m(.*)$/i){
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

      #abyss-fatoagp
      #scaffold1	1	4000	1	W	contigtest_0	1	4000	+
      #scaffold1	4001	4003	2	N	3	scaffold	yes	paired-ends
      #scaffold1	4004	14003	3	W	contigtest_1	1	10000	+
  
      #this script:
      #scaffold1	1	4000	1	W	oriscaffold_34	1	4000	+
      #scaffold1	4001	4003	2	N	3	scaffold	yes	paired-ends
      #scaffold1	4004	14003	3	W	oriscaffold_99	1	10000	-           (this will match the orientation infered by LINKS scaffolder)


      $end = $begin + $sz - 1;

      print "$a[0]\t$begin\t$end\t$ct\tW\toriscaffold_$id\t1\t$sz\t$ori\n";###object (ie. scaffold)

      if($b[$el+1] ne ""){###component (ie. gap)
         $ct++;
         my $gapstart = $end+1;
         my $gapend = $gap + $gapstart - 1;
	 print "$a[0]\t$gapstart\t$gapend\t$ct\tN\t$gap\tscaffold\tyes\tpaired-ends\n";
         $begin = $gapend + 1;
      }

      $el++;
   }
}
close IN;

exit;
