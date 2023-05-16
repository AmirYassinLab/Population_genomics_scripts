#!/usr/bin/perl -w
#Estimate Hudson-Slatkin-Maddison (1992) Fst from sync file for population pairs
#use strict;
use List::Util qw(sum);

#Path to the folder containing the *_sort.sync file
$genomes = '';
#Name of the *_sort.sync file without the extension .sync
$sync='';
#minimal depth
$min=0;
#Name of the output .fst file
my $outfile = ">".$genomes.$sync."_min".$min.".fst";
open(O,$outfile);
#Name of the input *_sort.sync file
my $infile = $genomes.$sync.".sync";
open(M, $infile);
#Parse the input file
open(M, $infile);
while (<M>){
  chomp;
  @line=split('\s',$_);
  @size=();
  @A=();
  @T=();
  @C=();
  @G=();
  for($p=4;$p<@line;$p++){
    @pop=split(':',$line[$p]);
    @strain=@pop[0..3];
    $size=sum(@strain);
    push @size, $size;
#Create an array for population counts for each allele at a site
    unless($size < $min){
      push @A, $strain[0];
      push @T, $strain[1];
      push @C, $strain[2];
      push @G, $strain[3];
      }
    }
#Skip sites with any population with 0:0:0:0:0:0 alleles
  if(scalar(@A) == @line-4){
##  if(scalar(@A) == 3){
#Estimate between-population heterozygosity, if all population are homozygous for the same allele, $H = 0
    $H = 1 - ((sum(@A)/sum(@size))**2) - ((sum(@T)/sum(@size))**2) - ((sum(@C)/sum(@size))**2) - ((sum(@G)/sum(@size))**2);
#Retain only variable sites, i.e. $H>0
    if ($H > 0){
      print "$line[0]\t$line[1]\t$line[2]\t$line[3]";
      print O "$line[0]\t$line[1]\t$line[2]\t$line[3]";
#Compare pairs of populations
#Define allele frequencies in pop1
      for($i=0;$i<@line-4;$i++){
##      for($i=0;$i<3;$i++){
        @p1=($A[$i]/$size[$i],$T[$i]/$size[$i],$C[$i]/$size[$i],$G[$i]/$size[$i]);
        @q1=(1-($A[$i]/$size[$i]),1-($T[$i]/$size[$i]),1-($C[$i]/$size[$i]),1-($G[$i]/$size[$i]));
#Define allele frequencies in pop2
        for($j=$i+1; $j < @line-4; $j++){
##        for($j=$i+1; $j < 3; $j++){
          @Hw=();
          @Hb=();
          @p2=($A[$j]/$size[$j],$T[$j]/$size[$j],$C[$j]/$size[$j],$G[$j]/$size[$j]);
          @q2=(1-($A[$j]/$size[$j]),1-($T[$j]/$size[$j]),1-($C[$j]/$size[$j]),1-($G[$j]/$size[$j]));
          for($b=0;$b<@p1;$b++){
#For each allele, estimate average within-population heterozygosity and between-population heterozygosity in each pair of populations
            push @Hw, ($p1[$b] * $q1[$b]) + ($p2[$b] * $q2[$b]);
            push @Hb, ($p1[$b] * $q2[$b]) + ($p2[$b] * $q1[$b]);
            }
#Estimate average within-population heterozygosity for all alleles
          $Hw=sum(@Hw);
#Estimate between-population heterozygosity for all alleles
          $Hb=sum(@Hb);
#If the two populations are homozygous for the same allele (i.e. $Hb = 0), set $Fst = 0
          if($Hb==0){
            $Fst=0;
            }
#Estimate Fst
          else{
            $Fst=1-($Hw/$Hb);
            }
          print "\t$Fst";
          print O "\t$Fst";
          }
        }
        print "\n";
        print O "\n";
      }
    }
  }
exit;

