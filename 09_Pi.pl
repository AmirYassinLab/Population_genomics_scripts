#!/usr/bin/perl -w
#Estimating pi from sync files
#use strict;
use List::Util qw(sum);

#Path to the folder containing the *_sort.sync file
$genomes = '';

#Name of the *_sort.sync file without the extension .sync
$sync='';
$min=0;
#Name of the output .fst file
my $outfile = ">".$genomes.$sync."_min".$min.".pi";
open(O,$outfile);
#Name of the input *_sort.sync file
my $infile = $genomes.$sync.".sync";
open(M, $infile);
#Parse the input file
while (<M>){
#Remove new line characters from the end of each line in the sync file
  chomp;
#Divide each line of the sync file into multiple elements (columns) depending on space characters, e.g.,  , \t, \r, \n or \f. This is because the join command could substitute \r with  
  @line=split('\s',$_);
#The sum of A, T, C, G alleles at a site for a population
  @size=();
#Sum of hereterozygote frequencies for all alleles at a site for each population
  @pi=();
#For each strain/population, note the first strain read counts is $line[4] and the last strain is $line[@line-1]
  for($p=4;$p<@line;$p++){
###  for($p=4;$p<7;$p++){
#Define the strain/population as an array @pop
    @pop=split(':',$line[$p]);
#Only retain the first four elements, i.e. A, T, C and G
    @strain=@pop[0..3];
#Estimate the population size, i.e. the sum of A, T, C and G
    $size=sum(@strain);
#Hereterozygote frequencies for all alleles at a site for a population
    @H=();
#Skip populations with 0:0:0:0:0:0 alleles at a site
	unless($size < $min){
      for($n=0;$n<@strain;$n++){
#For each allele, estimate the frequency of drawing it (p) times the frequency of drawing an alternative allele (1 - p)
        $H=($strain[$n]/$size)*(1-($strain[$n]/$size));
        push @H, $H;
        }
#Heterozygosity frequency at a site for a population
      $pi=sum(@H);
      push @pi, $pi;
      }
    }
#Print the output
  if(scalar(@pi) == @line-4){
##  if(scalar(@pi) == 3){
    print "$line[0]\t$line[1]\t$line[2]\t$line[3]";
    print O "$line[0]\t$line[1]\t$line[2]\t$line[3]";
    for($p=0;$p<@pi;$p++){
      print "\t$pi[$p]";
      print O "\t$pi[$p]";
      }
    print "\n";
    print O "\n";
    }
  }
exit;

