#!/usr/bin/perl -w
#Average Fst estimates over windows of predefined sizes
#use strict;
use List::Util qw(sum);

#Path to the folder containing the *.fst file
$genomes = '';
#Name of the .fst file without the extension .fst
$sync='';
#Determine the window size in bp
$win=10000;
#Estimate the window size in kb
$k=$win/1000;
#Name of the output .fst or .dxy file
my $outfile = ">".$genomes.$sync."_".$k."k.fst";

open(O,$outfile);
#Window increment
$w=1;
#Define the number of evaluated sites per a window, remember some sites were excluded due to 0:0:0:0:0:0 in the 11_Fst.pl script
$l=0;
#Name of the input .fst or dxy file
my $infile = $genomes.$sync.".fst";
open(M, $infile);
#Parse the input file
while (<M>){
  chomp;
  @line=split('\t',$_);
  @pi=@line[4..@line-1];
  if($line[2] > (($w * $win)-($win-1)) && $line[2] <= ($w * $win)){
    push @AoA, [@pi];
    $l++;
    $e=$line[0];
    $scaf=$line[1];
    $locus=$line[2];
    }
  if($line[2] > ($w * $win)){
    print "$e\t$scaf\t$w\t$locus\t$l";
    print O "$e\t$scaf\t$w\t$locus\t$l";
    for($c=0;$c<@line-4;$c++){
      for($r=0;$r<$l;$r++){
        push @PI, $AoA[$r][$c];
        }
      if(scalar(@PI) == 0){
        print "\tNA";
        print O "\tNA";
        @PI=();
        }
      else{
        $PI=sum(@PI)/scalar(@PI);
        print "\t$PI";
        print O "\t$PI";
        @PI=();
        }
      }
    print "\n";
    print O "\n";
    @AoA=();
    $l=0;
    $n=$w;
    $w=$n+int(($line[2]-($n * $win))/$win)+1;
    push @AoA, [@pi];
    $l++;
    $e=$line[0];
    $scaf=$line[1];
    $locus=$line[2];
    }
  if($line[1] ne $scaf){
    print "$e\t$scaf\t$w\t$locus\t$l";
    print O "$e\t$scaf\t$w\t$locus\t$l";
    for($c=0;$c<@line-4;$c++){
      for($r=0;$r<$l;$r++){
        push @PI, $AoA[$r][$c];
        }
      if(scalar(@PI) == 0){
        print "\tNA";
        print O "\tNA";
        @PI=();
        }
      else{
        $PI=sum(@PI)/scalar(@PI);
        print "\t$PI";
        print O "\t$PI";
        @PI=();
        }      }
    print "\n";
    print O "\n";
    @AoA=();
    $l=0;
    $w=1;
    if($line[2] > (($w * $win)-($win-1)) && $line[2] <= ($w * $win)){
      push @AoA, [@pi];
      $l++;
      $scaf=$line[1];
      }
    if($line[2] > ($w * $win)){
      $n=$w;
      $w=$n+int(($line[2]-($n * $win))/$win);
      push @AoA, [@pi];
      $l++;
      $e=$line[0];
      $scaf=$line[1];            
      $locus=$line[2];
      }
    }
  }
    print "$e\t$scaf\t$w\t$locus\t$l";
    print O "$e\t$scaf\t$w\t$locus\t$l";
    for($c=0;$c<@line-4;$c++){
      for($r=0;$r<$l;$r++){
        push @PI, $AoA[$r][$c];
        }
      if(scalar(@PI) == 0){
        print "\tNA";
        print O "\tNA";
        @PI=();
        }
      else{
        $PI=sum(@PI)/scalar(@PI);
        print "\t$PI";
        print O "\t$PI";
        @PI=();
        }      }
    print "\n";
    print O "\n";
exit;

