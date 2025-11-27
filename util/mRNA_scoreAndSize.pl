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

my $score = "";
my $size  = "";

while (<>) {
    my $aLine = $_;
    chomp($aLine);
    my @elem = split(/\t/, $aLine);

    if ($elem[2] eq "mRNA") {
        if ($size =~ /\S/) {
            print "$score\t$size\t$ne\n";
        }
        if ($elem[8] =~ /normalized_score=([^\;]+)/) {
            $score = $1;
            $size = $elem[4] - $elem[3] - 1;
        } else {
            die;
        }
        $ne = 0;
    } elsif ($elem[2] eq "CDS") {
        $ne++;
    } 
}

print "$score\t$size\t$ne\n";

__END__

