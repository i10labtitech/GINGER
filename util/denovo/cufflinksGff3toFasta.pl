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

my $gff3file   = $ARGV[0];
my $genomeFile = $ARGV[1];

open(IN, $genomeFile) or die "can not open a file \"$genomeFile\".";
while (<IN>) {
    if (/^\>\s*(\S+)/) {
        $id = $1;
    } else {
        my $aLine = $_;
        chomp($aLine);
        $genome{$id} .= join("", grep(/\S/, split(/\s*/, $aLine))); 
    } 
}
close(IN);

my $seq = "";
my $gId = "";
my $tId = "";

open(IN, $gff3file) or die "can not open a file \"$gff3file\".";
while (<IN>) {
    if (/^\s*$/) {
        if ($seq =~ /\S/) {
            print ">$tId $gId\n";
            print $seq, "\n";
        }
        $seq = "";
        next;
    }
    my $aLine = $_;
    chomp($aLine);
    my @elem = split(/\t/, $aLine);
    if ($elem[8] =~ /ID=GENE\^([^\,]+),TRANS\^([^\;]+)/) {
        $gId = $1;
        $tId = $2;
    } else {
        next;
    }
    my $tmpSeq .= substr($genome{$elem[0]}, $elem[3] - 1, $elem[4] - $elem[3] + 1);
    if ($elem[6] eq "+") {
        $seq .= $tmpSeq;
    } elsif ($elem[6] eq "-") {
        $tmpSeq = reverse($tmpSeq);
        $tmpSeq =~ tr/atgcATGC/tacgTACG/;
        $seq .= $tmpSeq;
    } else {
    }
}
close(IN);

if ($seq =~ /\S/) {
    print ">$tId $gId\n";
    print $seq, "\n";
}

__END__

