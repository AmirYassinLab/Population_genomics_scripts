#!/usr/bin/perl -w
#Aligning two reference genomes
#use strict;

#Path to minimap
$minimap = '';
#Path to the folder containing the assemblies files
$genomes = '';
#Name of the assemblies without .fasta extension
$assembly1='';
$assembly2='';
#Path to the folder to contain the paf file
$align = '';
#Name of the output files
$out="";

$cmd=$minimap."minimap2 -cx asm20 --cs ".$genomes.$assembly1.".fasta ".$genomes.$assembly2.".fasta > ".$align.$out.".paf";
system($cmd);
$cmd="sort -k6,6 -k8,8n ".$align.$out.".paf | k8 ".$minimap."paftools.js call -f ".$genomes.$assembly1.".fasta -L10000 -l1000 - > ".$align.$out.".vcf";
system($cmd);
exit;
