#!/usr/bin/perl -w
#Index the reference genome for minimap, picard and samtools
#use strict;
#Assumes the samtools executables have been copied to your system path

$cmd = "";
#Memory (in gigabytes) allocated to Picard
$mem = 8; 

#Name of reference genome fasta without .fasta
$reference = '';
#Path to minimap
$minimap = '';
#Path to samtools
$samtools = '';
#Path to all Picard java modules, which should all be in a single folder
$picard = '';

$cmd = $minimap."minimap2 -d " . $reference . ".mmi " . $reference . ".fasta";
system($cmd);
$cmd = "java -Xmx" . $mem . "g -jar " . $picard . " CreateSequenceDictionary REFERENCE=" . $reference . ".fasta  OUTPUT=" . $reference . ".dict";
system($cmd);
$cmd = "samtools faidx " . $reference . ".fasta";
system($cmd);
exit;
