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

if ($#ARGV != 1) {
    print "Usage:\n";
    print "perl weighten.pl [all.gff] [nextflow.config]\n";
    exit;
}
my $gffFile    = $ARGV[0];
my $configFile = $ARGV[1];

my $max_weight = 0;

##################################################
# Loading weights from nextflow.config
##################################################

my $flag   = 0;
my %w      = ();
open(IN, $configFile) or die "can not open a file [nextflow.config].";
while (<IN>) {
    s/\#.*$//;
    if (/MAPPING_WEIGHT\s*=\s*(\S+)/) {
        $w{mappingbase}  = $1;
    } elsif (/DENOVO_WEIGHT\s*=\s*(\S+)/) {
        $w{denovobase}   = $1;
    } elsif (/HOMOLOGY_WEIGHT\s*=\s*(\S+)/) {
        $w{homology}     = $1;
    } elsif (/AUGUSTUS_WEIGHT\s*=\s*(\S+)/) {
        $w{AUGUSTUS}     = $1;
    } elsif (/SNAP_WEIGHT\s*=\s*(\S+)/) {
        $w{SNAP}         = $1;
    } elsif (/(RNASEQ_OTHER\d+)_WEIGHT\s*=\s*(\S+)/) {
        $w{$1}           = $2;
    } elsif (/(HOMOLOGY_OTHER\d+)_WEIGHT\s*=\s*(\S+)/) {
        $w{$1}           = $2;
    } elsif (/(ABINITIO_OTHER\d+)_WEIGHT\s*=\s*(\S+)/) {
        $w{$1}           = $2;
    }
}
close(IN);

##################################################
# Setting $max_weight
##################################################

{
    foreach my $tag (keys(%w)) {
        if ($max_weight < $w{$tag}) { 
            $max_weight = $w{$tag};
        }
    }
}

##################################################
# Re-setting $w
##################################################

{
    foreach my $tag (keys(%w)) {
        $w{$tag} = $w{$tag} / $max_weight;
    }
}

##################################################
# Loading data form gff
# and setting new weights
##################################################

open(IN, $gffFile) or die "can not open a file [all.gff].";
while (<IN>) {
    my $aLine = $_;
    chomp($aLine);
    my @elem = split(/\t/, $aLine);
    if (($#elem > 8) & ($elem[0] != "#")) {
    } else {
        if ($elem[2] eq "CDS") {
            $elem[5] *= $w{$elem[1]};
            #                $elem[5] = "$elem[5]/$aVal";
            print join("\t", @elem), "\n";
        } else {
            print $_;
        }
    }
}
close(IN);

__END__
