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

my %additionalGff = ();
open(IN, $configFile) or die "can not open a file [nextflow.config].";
while (<IN>) {
    s/\#.*$//;
    if (/RNASEQ_OTHER(\d+)\s*=\s*\"(\S+)\"/) {
        $additionalGff{"RNASEQ_OTHER$1"}   = $2;
    } elsif (/HOMOLOGY_OTHER(\d+)\s*=\s*\"(\S+)\"/) {
        $additionalGff{"HOMOLOGY_OTHER$1"} = $2;
    } elsif (/ABINITIO_OTHER(\d+)\s*=\s*\"(\S+)\"/) {
        $additionalGff{"ABINITIO_OTHER$1"} = $2;
    }
}
close(IN);

##################################################
# Cat
##################################################

foreach my $aKey (keys(%additionalGff)) {
    open(IN, $additionalGff{$aKey}) or die "can not open a file \"$additionalGff{$aKey}\".";
    while (<IN>) {
        print $_;
    }
    close(IN);
}
