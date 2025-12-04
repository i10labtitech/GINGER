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

my $prevGeneId   = "";
my $prevContigId = "";
my @cds          = ();
my $modelNum     = 0;

my @elem = ();

while (<>) {
    last if (/^[\>\#]/);
    
    my $aLine = $_;
    chomp($aLine);
    my @elem = split(/\t/, $aLine);
    if ($elem[2] eq "CDS") {
        my $geneId = "";
        if ($elem[8] =~ /Parent=([^\;\s]+)/) {
            $geneId = $1;
        }
        if (($prevGeneId =~ /\S/) && ($prevGeneId ne $geneId)) {
            &outputCDS(\@cds);
            @cds = undef;
            @cds = ();
            push(@cds, \@elem);
            $prevContigId = $elem[0];
        } else {
            push(@cds, \@elem);
#            $prevContigId = $elem[0];
        }            
        $prevGeneId = $geneId;
    }
}

&outputCDS(\@cds);

sub outputCDS {
    my($c) = @_;
    
    if ($prevContigId ne $c->[0]->[0]) {
        print ">$c->[0]->[0]\n";
    }
    
    if ($c->[0]->[6] eq "+")  {
        @{$c} = sort{$a->[3] <=> $b->[3]} @{$c};  
    } elsif ($c->[0]->[6] eq "-")  {
        @{$c} = sort{$b->[3] <=> $a->[3]} @{$c};  
    } else {
        print STDERR "May be something wrong.\n";
    }
    
    if ($#{$c} == 0) {
        print "Esngl  ";
        if (($c->[0]->[6] eq "+") && ($c->[$i]->[3] <= $c->[$i]->[4])) {
            print $c->[$i]->[3], " ";
            print $c->[$i]->[4], " ";
        } elsif (($c->[0]->[6] eq "-") && ($c->[$i]->[3] <= $c->[$i]->[4])) {
            print $c->[$i]->[4], " ";
            print $c->[$i]->[3], " ";
        } else {
            print STDERR "May be something wrong.\n";
        }
        print "MODEL$modelNum\n";
    } else {
        for (my $i = 0; $i <= $#{$c}; $i++) {
            if ($i == 0) {
                print "Einit  ";
            } elsif ($i == $#{$c}) {
                print "Eterm  ";
            } else {
                print "Exon  ";
            }
            if (($c->[0]->[6] eq "+") && ($c->[$i]->[3] <= $c->[$i]->[4])) {
                print $c->[$i]->[3], " ";
                print $c->[$i]->[4], " ";
            } elsif (($c->[0]->[6] eq "-") && ($c->[$i]->[3] <= $c->[$i]->[4])) {
                print $c->[$i]->[4], " ";
                print $c->[$i]->[3], " ";
            } else {
                print STDERR "May be something wrong.\n";
            }
            print "MODEL$modelNum\n";
        }
    }
    $modelNum++;
    
=pod        
    } elsif ($c->[0]->[6] eq "-")  {
        for (my $i = $#{$c}; 0 <= $i; $i--) {
            if ($i == $#{$c}) {
                print "Einit  ";
            } elsif ($i == 0) {
                print "Eterm  ";
            } else {
                print "Exon  ";
            }
            print $c->[$i]->[4], " ";
            print $c->[$i]->[3], " ";
            print "MODEL$modelNum\n";
        }
        $modelNum++;
    } else {
        print STDERR "May be something wrong.\n";
    }
=cut

    $prevContigId = $c->[0]->[0];
} 
