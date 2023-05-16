#!/usr/bin/perl -w
#Aligning fastq files to a reference genome
#use strict;
#Assumes the samtools executables have been copied to your system path

$i = 0;
$cmd = "";
#Memory (in gigabytes) allocated to Picard
$mem = 8; 

#Name of reference genome fasta without .fasta
$reference = '';
#$reference = '/Volumes/Elements/Genomes/References/Tuberculatus_assemblage';
#Path to minimap
$minimap = '';
#Path to samtools
$samtools = '';
#Path to all Picard java modules, which should all be in a single folder
$picard = '';
#Path to the folder containing the *.fastq files
$genomes = '';
#List of the *.fastq files to be aligned. Remove the extension .fastq
@FastqFile = ();
#List extension of the fastq files
$f='';

for ($i = 0; $i < @FastqFile; $i++){
#Uncompress fastq files
#if single-end or concatenated
##  $cmd = 'gunzip ' . $genomes.$FastqFile[$i] . '.'.$f;
##  system($cmd);
#Create sam file
##  $cmd = $minimap."minimap2 -ax sr -t 16 " . $reference . ".fasta "  . $genomes.$FastqFile[$i] . ".".$f." -o " . $genomes.$FastqFile[$i] . ".sam";
##  system($cmd);

#if paired-end
##  $cmd = 'gunzip ' . $genomes.$FastqFile[$i] . '_1.'.$f;
##  system($cmd);
##  $cmd = 'gunzip ' . $genomes.$FastqFile[$i] . '_2.'.$f;
##  system($cmd);
#Create sam file
##  $cmd = $minimap."minimap2 -ax sr -t 16 " . $reference . ".fasta "  . $genomes.$FastqFile[$i] . "_1." . $f . " "  . $genomes.$FastqFile[$i] . "_2." . $f . " -o " . $genomes.$FastqFile[$i] . ".sam";
##  system($cmd);
  $cmd = 'gunzip ' . $genomes.$FastqFile[$i] . '.1.'.$f;
  system($cmd);
  $cmd = 'gunzip ' . $genomes.$FastqFile[$i] . '.2.'.$f;
  system($cmd);
#Create sam file
  $cmd = $minimap."minimap2 -ax sr -t 16 " . $reference . ".fasta "  . $genomes.$FastqFile[$i] . ".1." . $f . " "  . $genomes.$FastqFile[$i] . ".2." . $f . " -o " . $genomes.$FastqFile[$i] . ".sam";
  system($cmd);


#Create bam file
  $cmd = $samtools . "samtools view -bS " . $genomes.$FastqFile[$i] . ".sam > " . $genomes.$FastqFile[$i] . ".bam";
  system($cmd);
#Remove sam file
  $cmd = "rm " . $genomes.$FastqFile[$i] . ".sam";
  system($cmd);

#Compress fastq files
#if single-end or concatenated
##  $cmd = 'gzip ' . $genomes.$FastqFile[$i] . '.'.$f;
##  system($cmd);
#if paired-end
##  $cmd = 'gzip ' . $genomes.$FastqFile[$i] . '_1.'.$f;
##  system($cmd);
##  $cmd = 'gzip ' . $genomes.$FastqFile[$i] . '_2.'.$f;
##  system($cmd);
  $cmd = 'gzip ' . $genomes.$FastqFile[$i] . '.1.'.$f;
  system($cmd);
  $cmd = 'gzip ' . $genomes.$FastqFile[$i] . '.2.'.$f;
  system($cmd);

#Clean bam file, i.e. soft-clipping beyond-end-of-reference alignments and setting MAPQ to 0 for unmapped reads
  $cmd = "java -Xmx" . $mem . "g -jar " . $picard . "CleanSam INPUT=" . $genomes.$FastqFile[$i] . ".bam OUTPUT= " . $genomes.$FastqFile[$i] . "clean.bam";
  system($cmd);
#Remove bam file
  $cmd = "rm " . $genomes.$FastqFile[$i] . ".bam";
  system($cmd);
#Sort cleaned bam file by the reference sequence name (RNAME) field using the reference sequence dictionary (@SQ tag). Alignments within these subgroups are secondarily sorted using the left-most mapping position of the read (POS). 
  $cmd = "java -Xmx" . $mem . "g -jar " . $picard . "SortSam SORT_ORDER=coordinate INPUT=" . $genomes.$FastqFile[$i] . "clean.bam OUTPUT=" . $genomes.$FastqFile[$i] . "sort.bam";
  system($cmd);
#Remove cleaned bam file
  $cmd = "rm " . $genomes.$FastqFile[$i] . "clean.bam";
  system($cmd);
  }
exit;
