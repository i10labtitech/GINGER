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

my @val        = ();
my @hist       = ();
my @histSmooth = ();
my $z          = 2;
my $step       = 0.05;
my $min        = undef;
my $histMax    = 0;
my $firstMinScore = undef;
    
while (<>) {
    if (/^(\S+)/) {
        my $raw = $1;
        next if ($raw == 0.0);
        my $aVal = log($raw) / log(10.0);
        push(@val, $aVal);
        if (defined($min)) {
            $min = $aVal if ($aVal < $min);
        } else {
            $min = $aVal;
        }
    }
}

for (my $i = 0; $i <= $#val; $i++) {
    my $j = int(($val[$i] - $min) / $step);
    $hist[$j]++;
}

for (my $i = 0; $i <= $#hist; $i++) {
    $histMax = $hist[$i] if ($histMax < $hist[$i]);
}

for (my $i = 0; $i <= $#hist; $i++) {
    my $x0 = $i - $z;
    my $x1 = $i + $z;
    $x0 = 0 if ($x0 < 0);
    $x1 = $#hist if ($#hist < $x1);
    my $n = $x1 - $x0 + 1;
    my $s = 0;
    for (my $j = $x0; $j <= $x1; $j++) {
        $s += $hist[$j];
    }
    $histSmooth[$i] = $s / $n;
}

my $prevHistVal = undef;
my $prevSign    = undef;
for (my $i = $#hist; 0 <= $i; $i--) {
    if (defined($prevHistVal)) {
        my $sign = undef;
        if ($histMax / 100.0 < $histSmooth[$i]) {
            if ($prevHistVal - $histSmooth[$i] < 0) {
                $sign = -1;
            } else {
                $sign = 1;
            }
        }
        if (defined($sign)) {
            if (defined($prevSign)) {
                if (($prevSign == 1) && ($sign == -1)) {
                    my $aScore = 10.0 ** (($i + 1) * $step + $step / 2.0 + $min);
                    $firstMinScore = $aScore unless (defined($firstMinScore));
                }
            }
            $prevSign = $sign;
        }
    }
    $prevHistVal = $histSmooth[$i];
}

printf("%0.2f\n", $firstMinScore);
