package main;

=pod
Copyright (c) 2008, Brian Haas  
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=cut

## taken from CDNA::Overlap_assembler.pm
## really should be a more general purpose class as written here.

our $SEE = 0;

package Overlap_piler;


use strict;
use warnings;
use Carp;

sub new {
    my $packagename = shift;
    my $self = {
        node_list => []
        };
    bless ($self, $packagename);
    return ($self);
}


#### Static method!!!
sub simple_coordsets_collapser {
    my @coordsets = @_; ## list of coordinates [a, b], [c, d], ...
    
    my $counter = 0;
    
    my $piler = new Overlap_piler();

    my %coords_mapping;
    foreach my $coordset (@coordsets) {
        $counter++;
        $coords_mapping{$counter} = [@$coordset];
        
        my ($lend, $rend) = @$coordset;
        if ($lend !~ /\d/ || $rend !~ /\d/) {
            confess "Error, coordinates [ $lend, $rend ] include a non-number";
        }

        $piler->add_coordSet($counter, @$coordset);
    }

    my @clusters = $piler->build_clusters();
    
    my @coord_spans;
    foreach my $cluster (@clusters) {
        my @eles = @$cluster;
        my @coords;
        foreach my $ele (@eles) {
            push (@coords, @{$coords_mapping{$ele}});
        }

        @coords = sort {$a<=>$b} @coords;
        
        my $min_coord = shift @coords;
        my $max_coord = pop @coords;

        push (@coord_spans, [$min_coord, $max_coord]);
    }

    return (@coord_spans);
    
}



####
sub add_coordSet {
    my $self = shift;
    my ($acc, $end5, $end3) = @_;
    my ($lend, $rend) = sort {$a<=>$b} ($end5, $end3);
    my $node = CoordSet_node->new($acc, $lend, $rend);
    push (@{$self->{node_list}}, $node);
}




####
sub build_clusters {
    my $self = shift;
    my $node_list_aref = $self->{node_list};
    @{$node_list_aref} = sort {$a->{lend}<=>$b->{lend}} @{$node_list_aref}; #sort by lend coord.
    ## set indices
    for (my $i = 0; $i <= $#{$node_list_aref}; $i++) {
        $node_list_aref->[$i]->{myIndex} = $i;
    }
    
    my @clusters;
    my $first_node = $node_list_aref->[0];
    my $start_pos = 0;
    my ($exp_left, $exp_right) = ($first_node->{lend}, $first_node->{rend});
    print $first_node->{acc} . " ($exp_left, $exp_right)\n" if $SEE;
    for (my $i = 1; $i <= $#{$node_list_aref}; $i++) {
        my $curr_node = $node_list_aref->[$i];
        my ($lend, $rend) = ($curr_node->{lend}, $curr_node->{rend});
        print $curr_node->{acc} . " ($lend, $rend)\n" if $SEE;
        if ($exp_left <= $rend && $exp_right >= $lend) { #overlap
            $exp_left = &min($exp_left, $lend);
            $exp_right = &max($exp_right, $rend);
            print "overlap. New expanded coords: ($exp_left, $exp_right)\n" if $SEE;
        } else {
            print "No overlap; Creating cluster: " if $SEE;
            my @cluster;
            for (my $j=$start_pos; $j < $i; $j++) {
                my $acc = $node_list_aref->[$j]->{acc};
                push (@cluster, $acc);
                print "$acc, " if $SEE;
            }
            push (@clusters, [@cluster]);
            $start_pos = $i;
            ($exp_left, $exp_right) = ($lend, $rend);
            print "\nResetting expanded coords: ($lend, $rend)\n" if $SEE;
        }
    }
    
    print "# Adding final cluster.\n" if $SEE;
    if ($start_pos != $#{$node_list_aref}) {
        print "final cluster: " if $SEE;
        my @cluster;
        for (my $j = $start_pos; $j <= $#{$node_list_aref}; $j++) {
            my $acc = $node_list_aref->[$j]->{acc};
            print "$acc, " if $SEE;
            push (@cluster, $acc);
        }
        push (@clusters, [@cluster]);
        print "\n" if $SEE;
    } else {
        my $acc = $node_list_aref->[$start_pos]->{acc};
        push (@clusters, [$acc]);
        print "adding final $acc.\n" if $SEE;
    }
    return (@clusters);
}

sub min {
    my (@x) = @_;
    @x = sort {$a<=>$b} @x;
    my $min = shift @x;
    return ($min);
}

sub max {
    my @x = @_;
    @x = sort {$a<=>$b} @x;
    my $max = pop @x;
    return ($max);
}

#################################################################
package CoordSet_node;
use strict;

sub new {
    my $packagename = shift;
    my ($acc, $lend, $rend) = @_;
    my $self = { acc=>$acc,
                 lend=>$lend,
                 rend=>$rend,
                 myIndex=>undef(),
                 overlapping_indices=>[]
                 };
    bless ($self, $packagename);
    return ($self);
}

1; #EOM




