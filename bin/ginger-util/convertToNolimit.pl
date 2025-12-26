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

my $inFile = $ARGV[0];

#less /home/tany/src/TransDecoder-master/util/cdna_alignment_orf_to_genome_orf.pl
#less /home/tany/src/TransDecoder-master/PerlLib/

my $libPath = "";
if ($inFile =~ /^(\S+)\/cdna_alignment_orf_to_genome_orf\.pl/) {
    $libPath = "$1\/..\/PerlLib\/";
} else {
    die; 
}

open(IN, $inFile) or die "can not open a file \"$inFile\".";
while (<IN>) {
    if (/^use\s+lib/) {
        print $_;
        print "use lib (\"$libPath\");\n";
    } elsif (/\#\#\s+orf\s+orient\s+is\s+\'-\'/) {
        $flag = 1;
        print $_;
    } else {
        if ($flag == 1) {
            print "#$_";
            if (/^\s*\}\s*$/) {
                $flag = 0;
            }
        } else {
            print $_;
        }
    }
}
close(IN);
