#!/usr/bin/perl -w
#Convert read counts values in a sync file containing inbred lines into genotypes, for example 1145:1235:0:0:0:0 will become 1:1:0:0:0:0 for a heterozygous line and 21:0:0:0:0:0 will become 2:0:0:0:0:0 in a homozygous line
#use strict;
use List::Util qw(sum);

#Path to the folder containing the *.sync files
$genomes = '';
#Name of the input sync file
$sync='';
#Scaffolds to be retained, usually X, 2L, 2R, 3L and 3R
@Chr = ();
#Scaffolds length according to the reference fasta or gff files
#@Chrsize=();
#Cumulative value for each position according to the scaffolds order to be used in subsequent sort and join commands
@Chrcum=();
#Define output sync file
my $outfile = ">".$genomes.$sync."_genotype.sync";
open(O,$outfile);
#Define output sync file for triallelic sites
my $outfile1 = ">".$genomes.$sync."_genotype_t.sync";
open(O1,$outfile1);
#Enter input sync file
my $File = $genomes.$sync.'.sync';
open(M, $File);
#Define the minimum depth at a site
$min=10;
#Define the minimum ratio of an allele
$mina=0.25;
#Parse the input file
while (<M>){
#Remove new line characters from the end of each line in the sync file
  chomp;
#Divide each line of the sync file into multiple elements (columns)
#  @line=split("\t",$_);
  @line=split('\t',$_);
  for($c=0;$c<@Chr;$c++){
#Retain the desired scaffolds defined in @Chr
    if($line[0] eq $Chr[$c]){
	  @pos=();
#Push in the printable array (@pos) the cumulative value for the position
      $loc=$line[1]+$Chrcum[$c];
      push @pos, $loc;
#Then push the three first columns of the sync file, i.e. the scaffold, the position and the nucleotide at the reference genome
      for($l=0;$l<3;$l++){
        push @pos, $line[$l];
        }
#For each strain/population, note the first strain read counts is $line[3] and the last strain is $line[@line-1]
      for($p=3;$p<@line;$p++){
#Define the strain/population as an array @pop
        @pop=split(':',$line[$p]);
#Only retain the first four elements in @Nuc, i.e. A, T, C and G
        @strain=@pop[0..3];
#Do not genotype low-depth positions
        if(sum(@strain)<$min){
          push @pos, "0:0:0:0:$pop[4]:$pop[5]";
          }
#Genotype high-depth positions, alleles with read counts >= $mina ha
        else{
          $t=0;
          @genotype=();
          @homo=();
          $genotype='';
          for($b=0;$b<@strain;$b++){
            if(($strain[$b]/sum(@strain))>=$mina){
              push @genotype, 1;
              }
            else{
              push @genotype, 0;
              }
            }
#Homozygous sites
          if (sum(@genotype) == 1){
            @homo=(2*$genotype[0],2*$genotype[1],2*$genotype[2],2*$genotype[3]);
            push @homo, "$pop[4]:$pop[5]";
            $homo=join(':',@homo);
            push @pos, $homo;
			}
#Heterozygous sites
          elsif (sum(@genotype) == 2){
            push @genotype, "$pop[4]:$pop[5]";
            $genotype=join(':',@genotype);
            push @pos, $genotype;
			}
#Tri or quadri-allelic sites are omitted
          elsif (sum(@genotype) >2){
            $t++;
            push @pos, "0:0:0:0:$pop[4]:$pop[5]";
            }
          }
        }
      if($t > 0){
        print O1 "$_\n";
        }
      $pos=join("\t",@pos);
      print "$pos\n";
      print O "$pos\n";
      }
    }
  }   
exit;


