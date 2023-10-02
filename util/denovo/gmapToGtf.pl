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

my $gAcc = "";
my $tAcc = "";
my @lPos = ();
my @rPos = ();

while (<>) {
    if (/^>\s*(\S+)/) {
        my $tmpAcc = $1;
        if ($gAcc =~ /\S/) {
            &printData($gAcc, $tAcc, \@lPos, \@rPos);
        }
        $gAcc = "";
        $tAcc = $tmpAcc;
        undef(@lPos);
        undef(@rPos);
        @rPos = ();
    } elsif (/\s*[\+\-]([^\:]+):(\d+)-(\d+)\s+\(\d+-\d+\)\s+\d[^\%]+\%/) {
        $gAcc  = $1;
        push(@lPos, $2);
        push(@rPos, $3);
    }
}

if ($gAcc =~ /\S/) {
    &printData($gAcc, $tAcc, \@lPos, \@rPos);
}

sub printData {
    my($gAcc, $tAcc, $lPos, $rPos) = @_;

    my @sortedPos = sort {$a <=> $b} (@{$lPos}, @{$rPos});
    
    my $str = "";
    for (my $i = 0; $i <= $#{$lPos}; $i++) {
        if ($lPos->[$i] < $rPos->[$i]) {
            $str = "+";
            last;
        } elsif ($lPos->[$i] > $rPos->[$i]) {
            $str = "-";
            last;
        }
    }
    if ($str =~ /^\s*$/) {
        die;
    }

    print $gAcc, "\t";
    print ".\t";
    print "transcript\t";
    print $sortedPos[0], "\t";
    print $sortedPos[$#sortedPos], "\t";
    print ".\t";
    print $str, "\t";
    print ".\t";
    print "gene_id \"g\|$tAcc\"; transcript_id \"t\|$tAcc\"; name \"$tAcc\";\n";
    for (my $i = 0; $i <= $#{$lPos}; $i++) {
        print $gAcc, "\t";
        print ".\t";
        print "exon\t";
        if ($str eq "+") {
            print $lPos->[$i], "\t";
            print $rPos->[$i], "\t";
        } else {
            print $rPos->[$i], "\t";
            print $lPos->[$i], "\t";
        }
        print ".\t";
        print $str, "\t";
        print ".\t";
        print "gene_id \"g\|$tAcc\"; transcript_id \"t\|$tAcc\";\n";
    }
    print "\n";
    print "\n";
}

__END__

NC_003281.10    .       transcript      13614662        13622228        .       -       .       gene_id "g|TRINITY_DN61_c0_g1_i2"; transcript_id "t|TRINITY_DN61_c0_g1_i2"; name "TRINITY_DN61_c0_g1_i2";
NC_003281.10    .       exon    13622120        13622228        .       -       .       gene_id "g|TRINITY_DN61_c0_g1_i2"; transcript_id "t|TRINITY_DN61_c0_g1_i2";
NC_003281.10    .       exon    13621706        13622071        .       -       .       gene_id "g|TRINITY_DN61_c0_g1_i2"; transcript_id "t|TRINITY_DN61_c0_g1_i2";
NC_003281.10    .       exon    13621449        13621662        .       -       .       gene_id "g|TRINITY_DN61_c0_g1_i2"; transcript_id "t|TRINITY_DN61_c0_g1_i2";
NC_003281.10    .       exon    13615471        13615909        .       -       .       gene_id "g|TRINITY_DN61_c0_g1_i2"; transcript_id "t|TRINITY_DN61_c0_g1_i2";
NC_003281.10    .       exon    13614662        13615006        .       -       .       gene_id "g|TRINITY_DN61_c0_g1_i2"; transcript_id "t|TRINITY_DN61_c0_g1_i2";

