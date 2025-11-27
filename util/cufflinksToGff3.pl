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

my @data = ();
my $gid  = "";
my $tid  = "";

while (<>) {
    my $aLine = $_;
    chomp($aLine);
    my @elem = split(/\t/, $aLine);
    next if ($elem[2] ne "exon");
    if ($elem[8] =~ /gene_id\s+\"([^\"]+)\";\s*transcript_id\s+\"([^\"]+)\"/) {
        $gid = $1;
        $tid = $2;
    } else {
        next;
    }
    push(@{$data{$gid}}, {"tid"  => $tid,
                         "data" => $aLine,
                         "lpos" => $elem[3],
                         "str"  => $elem[6]});
}

foreach my $anGid (keys(%data)) {
    @{$data{$anGid}} = sort{$a->{lpos} <=> $b->{lpos}} @{$data{$anGid}};
    my $anTid = $data{$anGid}->[0]->{tid};
    if ($data{$anGid}->[0]->{str} eq "+") {
        my $orfLpos = 1;
        my $orfRpos;
        for (my $i = 0; $i <= $#{$data{$anGid}}; $i++) {
            my @elem = split("\t", $data{$anGid}->[$i]->{data});
            $orfRpos = $orfLpos + ($elem[4] - $elem[3] + 1) - 1;
            $elem[1] = "Cufflinks"; 
            $elem[2] = "match"; 
            $elem[5] = 100;
            $elem[8] = "ID=GENE\^$anGid,TRANS\^$anTid;Target=GENE\^$anGid,TRANS\^$anTid $orfLpos $orfRpos +";
            print join("\t", @elem), "\n";
            $orfLpos = $orfRpos + 1;
        }
    } elsif ($data{$anGid}->[0]->{str} eq "-") {
        my $orfLpos = 1;
        my $orfRpos;
        for (my $i = $#{$data{$anGid}}; 0 <= $i; $i--) {
            my @elem = split("\t", $data{$anGid}->[$i]->{data});
            $orfRpos = $orfLpos + ($elem[4] - $elem[3] + 1) - 1;
            $elem[1] = "Cufflinks"; 
            $elem[2] = "match"; 
            $elem[5] = 100;
            $elem[8] = "ID=GENE\^$anGid,TRANS\^$anTid;Target=GENE\^$anGid,TRANS\^$anTid $orfLpos $orfRpos +";
            print join("\t", @elem), "\n";
            $orfLpos = $orfRpos + 1;
        }
    }
    print "\n";
}

__END__
