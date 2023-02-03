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

package Gene_obj_indexer;
use strict;
use warnings;
use base qw(TiedHash);
use Gene_obj;
use Storable qw (thaw nfreeze);
use Carp;


####
sub new {
    my $packagename = shift;
    
    my $self = $packagename->SUPER::new(@_);
 
    return ($self);
   
}

####
sub store_gene {
    my ($self, $identifier, $gene_obj)  = @_;
    

    unless (ref $gene_obj) {
        confess "Error, no gene_obj as param";
    }
    
    my $blob = nfreeze ($gene_obj);
    
    my $success = 0;
    
    while (! $success) {
        $self->store_key_value($identifier, $blob);
        
        eval {
            my $gene_obj = $self->get_gene($identifier);
            
        };
        if ($@) {
            warn "error trying to store gene $identifier using berkeley db.  Trying again...\n";
        }
        else {
            # worked.
            $success = 1;
        }
    }
    
}


####
sub get_gene {
    my $self = shift;
    my $identifier = shift;

    my $blob = $self->get_value($identifier);
    
    unless ($blob) {
        confess "Error, no gene obj retrieved based on identifier $identifier";
    }

    my $gene_obj = thaw($blob);
    unless (ref $gene_obj) {
        confess "Error retrieving gene_obj based on identifier $identifier.  Data retrieved but not thawed properly.\n";
    }
    
    return ($gene_obj);
}


1; #EOM

