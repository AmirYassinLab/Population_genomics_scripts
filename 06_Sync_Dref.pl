#!/usr/bin/perl -w
#Detecting contamination by measuring divergence from reference genome in a sync file containing multiple genotyped strains
#use strict;
use List::Util qw(sum);

#Path to the folder containing the *.sync files
$genomes = '';
#Name of the input sync file
$sync='';
#List of the *sort.bam files that have been aligned in the sync file. Keep the same order
@FastqFile=();

#Define the nucleotides order in a sync file
@Nuc = ('A','T','C','G','N','-');
#Define output sync file
my $outfile = ">".$genomes.$sync."_genotype_sort2.dref";
open(O,$outfile);
#Enter input sync file
my $File = $genomes.$sync.'_genotype_sort.sync';
open(M, $File);
$k=0;#number of sites
$l=0;#number of 1-kb windows
$s=0;#number of evaluated sites, i.e. sites with no 0:0:0:0 nucleotides
#Parse the input file
while (<M>){
#Remove new line characters from the end of each line in the sync file
  chomp;
#Divide each line of the sync file into multiple elements (columns)
  @line=split('\t',$_);
  @pos=();
#To save computer memory (RAM), we will empty @AoA each 1000 sites
  if($k % 1000 == 0 && $k != 0){
    @diva=();
    @sitea=();
    for($f=0;$f<@FastqFile;$f++){
      $s=0;
      @divf=();
      for($i=0;$i<$k;$i++){
        unless($AoA[$i][$f] eq 'na'){
          push @divf, $AoA[$i][$f];
          $s++;
          }
        }
      push @diva, sum(@divf);
      push @sitea, $s;
      }
    push @DoD, [@diva];
    push @SoS, [@sitea];
    print "$line[1]\t$line[2]\t$DoD[$l][0]\t$SoS[$l][0]\n";
    @AoA=();
    $k=0;
    $l++;
    }
#Identify the reference nucleotide corresponding to the position $n. Ignore 'N' and '-'
  unless(uc($line[3]) eq 'N' || $line[3] eq '-'){
    for($n=0;$n<4;$n++){
      if(uc($line[3]) eq $Nuc[$n]){
        $nuc=$n;
        }
      }
    for($p=4;$p<@line;$p++){
#Define the strain/population as an array @pop
      @pop=split(':',$line[$p]);
      @strain=@pop[0..3];
      if(sum(@strain) == 0){
        push @pos, 'na';
        }
      else{
#Measure divergence and push it in the array @pos. Note if the strain has the same nucleotide as the reference divergence = 0.
#	    push @pos, 1-($strain[$nuc]/2);
	    push @pos, 1-($strain[$nuc]/sum(@strain));
	    }
	  }
#Create an array of arrays @AoA containing for all strains 'na' or divergence values
	push @AoA, [@pos];
	$k++;
	}
  }

@diva=();
@sitea=();
for($f=0;$f<@FastqFile;$f++){
  $s=0;
  @divf=();
  for($i=0;$i<$k;$i++){
    unless($AoA[$i][$f] eq 'na'){
      push @divf, $AoA[$i][$f];
      $s++;
      }
    }
  push @diva, sum(@divf);
  push @sitea, $s;
  }
push @DoD, [@diva];
push @SoS, [@sitea];
$l++;

for($f=0;$f<@FastqFile;$f++){
  @DoDa=();
  @SoSa=();
  for($i=0;$i<$l;$i++){
    push @DoDa, $DoD[$i][$f];
    push @SoSa, $SoS[$i][$f];
    }
  $sites= sum(@SoSa);
  $Dxy=(sum(@DoDa))/$sites;
  print "$FastqFile[$f]\t$Dxy\t$sites\n";
  print O "$FastqFile[$f]\t$Dxy\t$sites\n";
  }
exit;



