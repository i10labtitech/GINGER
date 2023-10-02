#!/usr/bin/env perl

# Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology
# 
# This file is part of GINGER.
# 
# GINGER is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# GINGER is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with GINGER; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

use strict;
use warnings;
use Carp;

my $configFile = $ARGV[0];

##################################################
# Loading weights from nextflow.config
##################################################

my $baseDir       = "";
my $oPrefix       = "";
my %w             = ();
my %additionalGff = ();
my $mappingGtf     = "";
my $mappingToMerge = "";
my $gmapResult     = "";
my $denovoResult   = "";
my $homologyResult = "";
my $augustusResult = "";
my $snapResult     = "";
my $cwd = `pwd`;
chomp($cwd);    
&loadData($configFile, $cwd, \$baseDir, \$oPrefix, \%w, \%additionalGff,
          \$mappingGtf, \$mappingToMerge, \$gmapResult, \$denovoResult,
          \$homologyResult, \$augustusResult, \$snapResult);

##################################################
# Output for config.ini
##################################################

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
foreach my $aKey (keys(%additionalGff)) {
    unless (-e $additionalGff{$aKey}) {
        print STDERR "No additional input file \"$aKey:$additionalGff{$aKey}\".\n";
    }
}

# file path
print "MAPPING_GTF=$mappingGtf\n";
print "MAPPING_TO_MERGE=$mappingToMerge\n";
print "GMAP_RESULT=$gmapResult\n";
print "DENOVO_RESULT=$denovoResult\n";
print "HOMOLOGY_RESULT=$homologyResult\n";
print "AUGUSTUS_RESULT=$augustusResult\n";
print "SNAP_RESULT=$snapResult\n";
foreach my $aKey (keys(%additionalGff)) {
    print "$aKey=$additionalGff{$aKey}\n";
}

# weight
printf("MAPPING_WEIGHT=%f\n",  $w{mappingbase});
printf("DENOVO_WEIGHT=%f\n",   $w{denovobase});
printf("HOMOLOGY_WEIGHT=%f\n", $w{homology});
printf("AUGUSTUS_WEIGHT=%f\n", $w{AUGUSTUS});
printf("SNAP_WEIGHT=%f\n",     $w{SNAP});
foreach my $aKey (keys(%additionalGff)) {
    if (($w{$aKey} =~ /^\s*\d+\s*$/) || ($w{$aKey} =~ /^\s*\d+\.\d+\s*$/)) {
        printf("$aKey\_WEIGHT=%f\n",   $w{$aKey});
    } else {
        print STDERR "No weight! --- $aKey\n"; 
    }
}

##################################################
# Sub routine
##################################################

sub loadData {
    my($configFile, $cwd, $baseDir, $oPrefix, $w, $additionalGff,
       $mappingGtf, $mappingToMerge, $gmapResult, $denovoResult,
       $homologyResult, $augustusResult, $snapResult) = @_;

    my $pdir = "";
    open(IN, $configFile) or die "can not open a file [nextflow.config].";
    while (<IN>) {
        s/\#.*$//;
        if (/^\s*PDIR\s*=\s*\"(\S+)\"/) {
            $pdir = $1;
        }
    }
    close(IN);

    open(IN, $configFile) or die "can not open a file [nextflow.config].";
    while (<IN>) {
        s/\#.*$//;
        if (/^\s*PDIR_PREP\s*=\s*\"(\S+)\"/) {
            my $fPath = &getFullpath($1, $pdir);
            $$baseDir = $fPath;
        } elsif (/^\s*OPREFIX\s*=\s*\"(\S+)\"/) {
            $$oPrefix = $1;
        } elsif (/^\s*RNASEQ_OTHER(\d+)\s*=\s*\"(\S+)\"/) {
            my $num = 1;
            my $fPath = &getFullpath($2, $pdir);
            $additionalGff->{"RNASEQ_OTHER$num"} = $fPath;
        } elsif (/^\s*HOMOLOGY_OTHER(\d+)\s*=\s*\"(\S+)\"/) {
            my $num = 1;
            my $fPath = &getFullpath($2, $pdir);
            $additionalGff->{"HOMOLOGY_OTHER$1"} = $fPath;
        } elsif (/^\s*ABINITIO_OTHER(\d+)\s*=\s*\"(\S+)\"/) {
            my $num = 1;
            my $fPath = &getFullpath($2, $pdir);
            $additionalGff->{"ABINITIO_OTHER$1"} = $fPath;
        } elsif (/^\s*MAPPING_WEIGHT\s*=\s*(\S+)/) {
            $w->{mappingbase}  = $1;
        } elsif (/^\s*DENOVO_WEIGHT\s*=\s*(\S+)/) {
            $w->{denovobase}   = $1;
        } elsif (/^\s*HOMOLOGY_WEIGHT\s*=\s*(\S+)/) {
            $w->{homology}     = $1;
        } elsif (/^\s*AUGUSTUS_WEIGHT\s*=\s*(\S+)/) {
            $w->{AUGUSTUS}     = $1;
        } elsif (/^\s*SNAP_WEIGHT\s*=\s*(\S+)/) {
            $w->{SNAP}         = $1;
        } elsif (/^\s*RNASEQ_OTHER(\d+)_WEIGHT\s*=\s*(\S+)/) {
            $w->{"RNASEQ_OTHER$1"}   = $2;
        } elsif (/^\s*HOMOLOGY_OTHER(\d+)_WEIGHT\s*=\s*(\S+)/) {
            $w->{"HOMOLOGY_OTHER$1"} = $2;
        } elsif (/^\s*ABINITIO_OTHER(\d+)_WEIGHT\s*=\s*(\S+)/) {
            $w->{"ABINITIO_OTHER$1"} = $2;
        } elsif (/^\s*STRINGTIE_RESULT\s*=\s*\"(\S+)\"/) {
            $$mappingGtf = "$cwd/$1";
        } elsif (/^\s*MAPPING_FINALRESULT\s*=\s*\"(\S+)\"/) {
            $$mappingToMerge = "$cwd/$1";
        } elsif (/^\s*GMAP_RESULT\s*=\s*\"(\S+)\"/) {
            $$gmapResult  = "$cwd/$1";
        } elsif (/^\s*DENOVO_FINALRESULT\s*=\s*\"(\S+)\"/) {
            $$denovoResult  = "$cwd/$1";
        } elsif (/^\s*HOMOLOGY_FINALRESULT\s*=\s*\"(\S+)\"/) {
            $$homologyResult  = "$cwd/$1";
        } elsif (/^\s*AUGUSTUS_RESULT\s*=\s*\"(\S+)\"/) {
            $$augustusResult = "$cwd/$1";
        } elsif (/^\s*SNAP_RESULT\s*=\s*\"(\S+)\"/) {
            $$snapResult = "$cwd/$1";
        }
    }
    close(IN);
}

sub getFullpath {
    my($path, $pdir) = @_;
    
    $path =~ s/\$\{PDIR\}/$pdir\//;
    
    return $path;
}
