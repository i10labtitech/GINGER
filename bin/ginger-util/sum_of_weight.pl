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

if ($#ARGV != 0) {
    print "Usage:\n";
    print "perl sum_of_weight.pl [nextflow.config]\n";
    exit;
}
my $configFile = $ARGV[0];

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
    if ($hwFlag != 1) {
        foreach my $tag (keys(%w)) {
            if ($max_weight < $w{$tag}) { 
                $max_weight = $w{$tag};
            }
        }
    } elsif ($hwFlag == 1) {
        foreach my $tag (keys(%w)) {
            next if ($tag eq "homology");
            if ($max_weight < $w{$tag}) { 
                $max_weight = $w{$tag};
            }
        }
        foreach my $tag (keys(%hwOrg)) {
            if ($max_weight < $hwOrg{$tag}) { 
                $max_weight = $hwOrg{$tag};
            }
        }
    }
}

##################################################
# 
# 
##################################################

my $sum_of_weight= 0;
{
    if ($hwFlag != 1) {
        foreach my $tag (keys(%w)) {
#            print "
#            $sum_of_weight += $w{$tag} / $max_weight;
#            ";
            $sum_of_weight += $w{$tag} / $max_weight;
        }
    } elsif ($hwFlag == 1) {
        foreach my $tag (keys(%w)) {
            next if ($tag eq "homology");
#            print "
#            $sum_of_weight += $w{$tag} / $max_weight;
#            ";
            $sum_of_weight += $w{$tag} / $max_weight;
        }
        my $homology_max_weight = 0;
        foreach my $tag (keys(%hwOrg)) {
            if ($homology_max_weight < $hwOrg{$tag}) { 
                $homology_max_weight = $hwOrg{$tag};
            }
        }
#        print "
#        $sum_of_weight += $homology_max_weight / $max_weight;
#        ";
        $sum_of_weight += $homology_max_weight / $max_weight;
    }
}

print $sum_of_weight, "\n";

__END__
