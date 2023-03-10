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

our $SEE;

package GTF_to_geneobjs;

use strict;
use warnings;
use Gene_obj;
use Gene_obj_indexer;
use Carp;
use Data::Dumper;

sub parse_file {
    my ($gtf_filename, $gene_obj_indexer) = @_;
    ## note can use gene_obj_indexer or provide a hash reference.
    
    my $hash_mode = 0;
    if (ref $gene_obj_indexer eq 'HASH') {
        $hash_mode = 1;
    }
    

    my %seqname_map;

    my $gene_objs = GTF_to_gene_objs($gtf_filename);
    
    for my $gene_obj (@$gene_objs) {

        my $gene_id = $gene_obj->{TU_feat_name};
        my $seqname = $gene_obj->{asmbl_id};
        
        if ($hash_mode) {
            $gene_obj_indexer->{$gene_id} = $gene_obj;
        }
        else {
            $gene_obj_indexer->store_gene($gene_id, $gene_obj);
        }
        
        # add to gene list for asmbl_id
        my $gene_list = $seqname_map{$seqname};
        unless (ref $gene_list) {
            $gene_list = $seqname_map{$seqname} = [];
        }
        push (@$gene_list, $gene_id);
    }
    return (\%seqname_map); # contig_id -> [ list of gene_ids]
}

sub GTF_to_gene_objs {
    my ($gtf_filename) = @_;

    my %gene_transcript_data;

    open (my $fh, $gtf_filename) or die "Error, cannot open $gtf_filename";
    while (<$fh>) {
        unless (/\w/) { next; }
        if (/^\#/) { next; } # comment line.
        
        chomp;
        my ($seqname, $source, $type, $lend, $rend, $score, $strand, $gtf_phase, $annot) = split (/\t/);
        
        my ($end5, $end3) = ($strand eq '+') ? ($lend, $rend) : ($rend, $lend);
        
        $annot =~ /gene_id \"([^\"]+)\"/;
        my $gene_id = $1 or confess "Error, cannot get gene_id from $annot of line\n$_";
        
        $annot =~ /transcript_id \"([^\"]+)\"/;
        my $transcript_id = $1 or confess "Error, cannot get transcript_id from $annot of line\n$_";
        
        if ($type eq 'CDS' || $type eq 'stop_codon') {
            $gene_transcript_data{$seqname}->{$gene_id}->{$transcript_id}->{CDS}->{$end5} = $end3;
            $gene_transcript_data{$seqname}->{$gene_id}->{$transcript_id}->{mRNA}->{$end5} = $end3;
        }
        if ($type =~ /UTR/) {
            $gene_transcript_data{$seqname}->{$gene_id}->{$transcript_id}->{mRNA}->{$end5} = $end3;
        }
    }
    close $fh;
    
    ## create gene objects.
 
    my @top_gene_objs;
    
    foreach my $seqname (keys %gene_transcript_data) {

        my $genes_href = $gene_transcript_data{$seqname};
        
        foreach my $gene_id (keys %$genes_href) {
            
            my $transcripts_href = $genes_href->{$gene_id};
            
            my @gene_objs;

            foreach my $transcript_id (keys %$transcripts_href) {

                my $coord_types_href = $transcripts_href->{$transcript_id};

                my $CDS_coords_href = $coord_types_href->{CDS};
                my $mRNA_coords_href = $coord_types_href->{mRNA};
                

                my $gene_obj = new Gene_obj();
                $gene_obj->populate_gene_object($CDS_coords_href, $mRNA_coords_href);
                
                $gene_obj->{TU_feat_name} = $gene_id;
                $gene_obj->{Model_feat_name} = $transcript_id;
                $gene_obj->{com_name} = $transcript_id;
                $gene_obj->{asmbl_id} = $seqname;
                
                $gene_obj->join_adjacent_exons();

                push (@gene_objs, $gene_obj);
            }
        
                    
            ## want single gene that includes all alt splice variants here
            if(scalar(@gene_objs)) {
                my $template_gene_obj = shift @gene_objs;
                foreach my $other_gene_obj (@gene_objs) {
                    $template_gene_obj->add_isoform($other_gene_obj);
                }
                push (@top_gene_objs, $template_gene_obj);       
            
                print $template_gene_obj->toString() if $SEE; 
            }
        }
    }
    return (\@top_gene_objs);
}


1; #EOM
