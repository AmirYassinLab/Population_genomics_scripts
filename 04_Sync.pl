#!/usr/bin/perl -w
#Generating a synchronized file from sorted bam files
#use strict;

$i = 0;
$cmd = "";
#Abbreviation of the reference genome to be added at the beginning of the sync file name
$ref = "";
#Abbreviation of the population name if the multiple sorted bam files are from the same population
$pop= "";

#Name of reference genome fasta without .fasta
$reference = '';
#Path to PoPoolation2
$popoolation = '';
#Path to the folder containing the *sort.bam files
$genomes = '';
#List of the *sort.bam files to be aligned. Keep the extension sort.bam
@FastqFile=();
#Create the mpileup command
$cmd = "samtools mpileup -f " . $reference . ".fasta -B ";
push @cmd, $cmd;
for ($i = 0; $i < @FastqFile; $i++){
  push @cmd, $genomes.$FastqFile[$i] . " ";
  }
push @cmd, ">";
#Create the mpileup/sync output filename to be used in the mpileup and synchronizing commands
push @mp, $genomes;
push @mp, $ref."_";
#if you want the mpileup/sync file name to include the name of all sort.bam files
#for ($i = 0; $i < @FastqFile-1; $i++){
#  push @mp, $FastqName[$i] . "_";
#  }
#push @mp, $FastqName[@FastqName-1].".mpileup";
#if you want the mpileup/sync file name to include the abbreviation of the population name
push @mp, $pop.".mpileup";
$mp=join('',@mp);
push @cmd, $mp;
$cmd=join('',@cmd);
system($cmd);

#Create the synchronizing command
@cmd=();
@mp=();
$cmd = "java -ea -Xmx7g -jar " . $popoolation . "mpileup2sync.jar --input ";
push @cmd, $cmd;
push @cmd, $mp;
$cmd = " --output ";
push @cmd, $cmd;
@mp=split('\.',$mp);
push @cmd, $mp[0];
$cmd=".sync --fastq-type sanger --min-qual 20 --threads 8";
push @cmd,$cmd;
$cmd=join('',@cmd);
system($cmd);

exit;
