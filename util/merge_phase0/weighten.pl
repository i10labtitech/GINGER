#!/usr/bin/perl

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

my $hwFlag = 0;

##################################################
# Loading weights from nextflow.config
##################################################

my $flag   = 0;
my %w      = ();
my %hwOrg  = ();
my $orgChk = ();
open(IN, $configFile) or die "can not open a file [nextflow.config].";
while (<IN>) {
    s/\#.*$//;
    if (/HOMOLOGY_DATA\s*=\s*\[/) {
        $flag = 1;
    } elsif (/MAPPING_WEIGHT\s*=\s*(\S+)/) {
        $w{mappingbase}  = $1;
    } elsif (/DENOVO_WEIGHT\s*=\s*(\S+)/) {
        $w{denovobase}   = $1;
    } elsif (/HOMOLOGY_WEIGHT\s*=\s*(\S+)/) {
        $w{homology}     = $1;
    } elsif (/AUGUSTUS_WEIGHT\s*=\s*(\S+)/) {
        $w{AUGUSTUS}     = $1;
    } elsif (/SNAP_WEIGHT\s*=\s*(\S+)/) {
        $w{SNAP}         = $1;
    } elsif (/HOMOLOGY_WEIGHT.(\S+)\s*=\s*(\S+)/) {
        $hwOrg{$1}       = $2;
    } else {
        if (($flag == 1) && (/\"PREFIX\"/)) {
            if (/\"PREFIX\"\s*\=\s*\"(\S+)\"/) {
                $flag = 3;
                $orgChk{$1} = 1;
            } elsif (/\"PREFIX\"\s*\:\s*\[/) {
                $flag = 2;
                if (/\[\s*\"(\S+)\"/) {
                    $orgChk{$1} = 1;
                }
            }
        } elsif ($flag == 2) {
            if (/\"(\S+)\"/) {
                $orgChk{$1} = 1;
            }
            if (/\]/) {
                $flag = 3;
            }
        }
    }
}
close(IN);

##################################################
# Checking the number or organisms,
# if HOMOLOGY_WEIGHT.* are set
##################################################

if (0 < scalar(keys(%hwOrg))) {
    $hwFlag = 1;
}

if ($hwFlag == 1) {
    my $orgChkFlag = 0;
    foreach my $anOrg (keys(%orgChk)) {
        if ($hwOrg{$anOrg} =~ /\S/)  {
        } else {
            print STDERR "No weight for $anOrg.\n";
            $orgChkFlag = 1;
        }
    }
    foreach my $anOrg (keys(%hwOrg)) {
        if ($orgChk{$anOrg} != 1) {
            print STDERR "$anOrg is not found in targets in \"Homology based\".\n";
            $orgChkFlag = 1;
        }
    }
    if ($orgChkFlag == 1) {
        die "error.";
    }
}

##################################################
# Setting $max_weight
##################################################

{
    foreach my $tag (keys(%w)) {
        if ($max_weight < $w{$tag}) { 
            $max_weight = $w{$tag};
        }
    }
    if ($hwFlag == 1) {
        foreach my $tag (keys(%hwOrg)) {
            if ($max_weight < $hwOrg{$tag}) { 
                $max_weight = $hwOrg{$tag};
            }
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
    if ($hwFlag == 1) {
        foreach my $tag (keys(%hwOrg)) {
            $hwOrg{$tag} = $hwOrg{$tag} / $max_weight;
        }
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
            if ($hwFlag == 1) {
                if ($elem[1] eq "homology") {
                    if ($elem[8] =~ /ID=cds\.mRNA\d+\.([^\;]+)/) {
                        my $aVal = $elem[5] * $hwOrg{$1};
                        #                        $elem[5] = "$elem[5]/$aVal/$1";
                        $elem[5] = $aVal;
                    } else {
                        print STDERR "Something wrong \"$aLine\"\n.";
                    }
                } else {
                    my $aVal = $elem[5] * $w{$elem[1]};
                    #                    $elem[5] = "$elem[5]/$aVal";
                    $elem[5] = $aVal;
                }
            } else {
                my $aVal = $elem[5] * $w{$elem[1]};
                #                $elem[5] = "$elem[5]/$aVal";
                $elem[5] = $aVal;
            }
            print join("\t", @elem), "\n";
        } else {
            print $_;
        }
    }
}
close(IN);

__END__
