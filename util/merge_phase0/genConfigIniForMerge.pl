#!/usr/bin/env perl

use strict;
use warnings;
use Carp;

my $baseDir  = $ARGV[0];
my $oPrefix  = $ARGV[1];
my $mWeight  = $ARGV[2] || 1.6; # MAPPING_WEIGHT  (1.6)
my $dWeight  = $ARGV[3] || 1.4; # DENOVO_WEIGHT   (1.4)
my $hWeight  = $ARGV[4] || 2.0; # HOMOLOGY_WEIGHT (2.0)
my $aWeight  = $ARGV[5] || 1.8; # AUGUSTUS_WEIGHT (1.8)
my $sWeight  = $ARGV[6] || 1.2; # SNAP_WEIGHT     (1.2)

my $mappingGtf     = "$baseDir/Prep/mapping/mapping/$oPrefix.gtf";
my $mappingToMerge = "$baseDir/Prep/mapping/merge/$oPrefix\_merge.gff3";
my $gmapResult     = "$baseDir/Prep/denovo/$oPrefix\_denovo.cd.fasta.gmap.filtered";
my $denovoResult   = "$baseDir/Prep/denovo/$oPrefix\_final.gff3";
my $homologyResult = "$baseDir/Prep/homology/homology_filter/all_homology_filter.gff";
my $augustusResult = "$baseDir/Prep/abinitio/augusutus/Augustus_abinitio.gtf";
my $snapResult     = "$baseDir/Prep/abinitio/snap/$oPrefix.final.gff";

my $flag = 0;
foreach my $anInputFile ($mappingGtf, $mappingToMerge, $gmapResult, $denovoResult, $homologyResult, $augustusResult, $snapResult) {
    unless (-e $anInputFile) {
        print STDERR "No input file \"$anInputFile\".\n";
        $flag = 1;
    }
}
if ($flag == 1) {
    die "Not enough input files.";
}

# file path
print "MAPPING_GTF=$mappingGtf\n";
print "MAPPING_TO_MERGE=$mappingToMerge\n";
print "GMAP_RESULT=$gmapResult\n";
print "DENOVO_RESULT=$denovoResult\n";
print "HOMOLOGY_RESULT=$homologyResult\n";
print "AUGUSTUS_RESULT=$augustusResult\n";
print "SNAP_RESULT=$snapResult\n";

# weight
printf("MAPPING_WEIGHT=%f\n",  $mWeight);
printf("DENOVO_WEIGHT=%f\n",   $dWeight);
printf("HOMOLOGY_WEIGHT=%f\n", $hWeight);
printf("AUGUSTUS_WEIGHT=%f\n", $aWeight);
printf("SNAP_WEIGHT=%f\n",     $sWeight);



