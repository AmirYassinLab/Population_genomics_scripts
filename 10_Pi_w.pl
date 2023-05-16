#!/usr/bin/perl -w
#Average pi estimates over windows of predefined sizes
#use strict;
use List::Util qw(sum);

#Path to the folder containing the *.sync file
$genomes = '';
#Name of the .pi file without the extension .pi
$sync='';
#Determine the window size in bp
$win=10000;
#Estimate the window size in kb
$k=$win/1000;
#Name of the output .pi file
my $outfile = ">".$genomes.$sync."_".$k."k.pi";
open(O,$outfile);
#Window increment
$w=1;
#Define the number of evaluated sites per a window, remember some sites were excluded due to 0:0:0:0:0:0 in the 09_Pi.pl script
$l=0;
$poly=0;
#Name of the input .pi file
my $infile = $genomes.$sync.".pi";
open(M, $infile);
#Parse the input file
while (<M>){
#Remove new line characters from the end of each line in the input file
  chomp;
#Divide each line of the sync file into multiple elements (columns)
  @line=split('\t',$_);
#Define a @pi array containing pi values for each population
  @pi=@line[4..(@line-1)];
#Define window limits according to position in the scaffold (i.e. $line[2]: window start = ($w * $win)-($win - 1) and window end = $w * $win, e.g. for the first ($w = 1) 10-k window: start = 1 and end = 10000 corresponding to X:1..10000 
  if($line[2] > (($w * $win)-($win-1)) && $line[2] <= ($w * $win)){
#Push all pi values for that site in an array of arrays @AoA
    push @AoA, [@pi];
    $l++;
    if(sum(@pi) != 0){
      $poly++;
      }
    $e=$line[0];
    $scaf=$line[1];
    $locus=$line[2];
    }
#We exceed the end of the window, so average pi values in @AoA and then empty it to be filled with values of the next window
  if($line[2] > ($w * $win)){
    print "$e\t$scaf\t$w\t$locus\t$l\t$poly";
    print O "$e\t$scaf\t$w\t$locus\t$l\t$poly";
#$c=columns of @AoA
    for($c=0;$c<@line-4;$c++){
#$r=rows of @AoA
      for($r=0;$r<$l;$r++){
#push pi values for all sites for a population in an array
        push @PI, $AoA[$r][$c];
        }
#Average pi values for all sites for a population in an array
      $PI=sum(@PI)/scalar(@PI);
      print "\t$PI";
      print O "\t$PI";
      @PI=();
      }
    print "\n";
    print O "\n";
    @AoA=();
#Start the new window @AoA
    $l=0;
    $poly=0;
    $n=$w;
    $w=$n+int(($line[2]-($n * $win))/$win)+1;
    push @AoA, [@pi];
    $l++;
    if(sum(@pi) != 0){
      $poly++;
      }
    $e=$line[0];
    $scaf=$line[1];
    $locus=$line[2];
    }
#When changing the scaffold, average pi values in @AoA if the window is < $win
  if($line[1] ne $scaf){
    print "$e\t$scaf\t$w\t$locus\t$l\t$poly";
    print O "$e\t$scaf\t$w\t$locus\t$l\t$poly";
    for($c=0;$c<@line-4;$c++){
      for($r=0;$r<$l;$r++){
        push @PI, $AoA[$r][$c];
        }
      $PI=sum(@PI)/scalar(@PI);
      print "\t$PI";
      print O "\t$PI";
      @PI=();
      }
    print "\n";
    print O "\n";
    @AoA=();
    $l=0;
    $poly=0;
    $w=1;
    if($line[2] > (($w * $win)-($win-1)) && $line[2] <= ($w * $win)){
      push @AoA, [@pi];
      $l++;
    if(sum(@pi) != 0){
      $poly++;
      }
      $scaf=$line[1];
      }
    if($line[2] > ($w * $win)){
      $n=$w;
      $w=$n+int(($line[2]-($n * $win))/$win);
      push @AoA, [@pi];
      $l++;
      $poly++;
      $e=$line[0];
      $scaf=$line[1];            
      $locus=$line[2];
      }
    }
  }
#When parsing is over, average pi values in @AoA for the last window for the last scaffold if it is < $win
    print "$e\t$scaf\t$w\t$locus\t$l\t$poly";
    print O "$e\t$scaf\t$w\t$locus\t$l\t$poly";
    for($c=0;$c<@line-4;$c++){
      for($r=0;$r<$l;$r++){
        push @PI, $AoA[$r][$c];
        }
      $PI=sum(@PI)/scalar(@PI);
      print "\t$PI";
      print O "\t$PI";
      @PI=();
      }
    print "\n";
    print O "\n";
exit;

