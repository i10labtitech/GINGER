#!/usr/local/bin/perl

=pod
Copyright (c) 2008, Brian Haas  
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=cut

package TiedHash;
use strict;
use warnings;
use DB_File;
use Carp;


sub new {
    my $packagename = shift;
    
    my $prefs_href = shift;
    
    if ($prefs_href && ! ref $prefs_href) {
        confess "Error, need hash reference with opts in constructor.\n";
    }
    

    my $self = { 
        index_filename => undef,
        tied_index => {},
        tie_invoked => 0,
    };
    
    bless ($self, $packagename);
    

    if (ref $prefs_href eq "HASH") {
        if (my $index_file = $prefs_href->{"create"}) {
            $self->create_index_file($index_file);
        }
        elsif ($index_file = $prefs_href->{"use"}) {
            $self->use_index_file($index_file);
        }
    }
            
    
    return ($self);
}

####
sub tie_invoked {
    my $self = shift;
    return ($self->{tie_invoked});
}


####
sub DESTROY {
    my $self = shift;
    if ($self->{index_filename}) {
        # hash must have been tied
        # so, untie it
        untie (%{$self->{tied_index}});
    }
}


####
sub create_index_file {
    my $self = shift;
    return ($self->make_index_file(@_));
}



####
sub make_index_file {
    my $self = shift;
    my $filename = shift;
    
    unless ($filename) {
        confess "need filename as parameter";
    }

    if (-e $filename) {
        unlink $filename or confess "cannot remove existing index filename $filename";
    }
    
    $self->{index_filename} = $filename;
    
    tie (%{$self->{tied_index}}, 'DB_File', $filename, O_CREAT|O_RDWR, 0666, $DB_BTREE);

    $self->{tie_invoked} = 1;
    
    return;
}


####
sub use_index_file {
    my $self = shift;
    my $filename = shift;
    
    unless ($filename) {
        confess "need filename as parameter";
    }
    
    unless (-s $filename) {
        confess "Error, cannot locate file: $filename\n";
    }
    
    $self->{index_filename} = $filename;
    
    tie (%{$self->{tied_index}}, 'DB_File', $filename, O_RDONLY, 0, $DB_BTREE);

    $self->{tie_invoked} = 1;

    my @keys = $self->get_keys();
    unless (@keys) {
        confess "Error, tried using $filename db, but couldn't perform retrievals.\n";
    }
    
    return;

}


####
sub store_key_value {
    my ($self, $identifier, $value)  = @_;
    
    #my $num_keys = scalar ($self->get_keys());
    
    unless ($self->tie_invoked()) {
        confess "Error, cannot store key/value pair since tied hash not created.\n";
    }


    my $found = 0;
    while (! $found) {
        $self->{tied_index}->{$identifier} = $value;
        
        my $val = $self->get_value($identifier);
        if (defined $val) {
            $found = 1;
        }
        else {
            warn "Berkeley DB had trouble storing ($identifier); trying again.\n";
        }
    }
    
    return;
    
}


####
sub get_value {
    my $self = shift;
    my $identifier = shift;

      
    unless ($self->tie_invoked()) {
        confess "Error, cannot retrieve value from untied hash\n";
    }

    my $value = $self->{tied_index}->{$identifier};
    
    return ($value);
}


## 
sub get_keys {
    my $self = shift;
    
    unless ($self->tie_invoked()) {
        confess "Error, cannot retrieve values from untied hash\n";
    }
    
    return (keys %{$self->{tied_index}});
}


1; #EOM
