#!/usr/bin/env perl

# Copyright (c) 2012, The Broad Institute, Inc. All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
# 
# ·         Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
# 
# ·         Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
# 
# ·         Neither the name of the Broad Institute nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.**
# 
# THIS SOFTWARE IS PROVIDED BY THE BROAD INSTITUTE  ''AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, 
# THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE BROAD INSTITUTE 
# BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES(INCLUDING, BUT NOT LIMITED TO, 
#                                                                                              PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, 
# EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

use strict;
use warnings;

use FindBin;
use lib ("$FindBin::Bin/../../lib/perl5/mapping/");
use Gene_obj;

my $usage = "usage: $0 cufflinks.gtf\n\n";

my $cufflinks_gtf = $ARGV[0] or die $usage;


main: {
	
	my %genome_trans_to_coords;
	
	open (my $fh, $cufflinks_gtf) or die "Error, cannot open file $cufflinks_gtf";
	while (<$fh>) {
		chomp;

        if (/^\#/) { next; }
		unless (/\w/) { next; }
		
		my @x = split(/\t/);
		
		my $scaff = $x[0];
		my $type = $x[2];
		my $lend = $x[3];
		my $rend = $x[4];

		my $orient = $x[6];
		
		my $info = $x[8];
		
		unless ($type eq 'exon') { next; }

		my @parts = split(/;/, $info);
		my %atts;
		foreach my $part (@parts) {
			$part =~ s/^\s+|\s+$//;
			$part =~ s/\"//g;
			my ($att, $val) = split(/\s+/, $part);
			unless (defined $att) { next; }
            
			if (exists $atts{$att}) {
				die "Error, already defined attribute $att in $_";
			}
			
			$atts{$att} = $val;
		}

		my $gene_id = $atts{gene_id} or die "Error, no gene_id at $_";
		my $trans_id = $atts{transcript_id} or die "Error, no trans_id at $_";
		
		my ($end5, $end3) = ($orient eq '+') ? ($lend, $rend) : ($rend, $lend);

		$genome_trans_to_coords{$scaff}->{$gene_id}->{$trans_id}->{$end5} = $end3;

	}


	## Output genes in gff3 format:

	foreach my $scaff (sort keys %genome_trans_to_coords) {

		my $genes_href = $genome_trans_to_coords{$scaff};

		foreach my $gene_id (keys %$genes_href) {

			my $trans_href = $genes_href->{$gene_id};

			foreach my $trans_id (keys %$trans_href) {

				my $coords_href = $trans_href->{$trans_id};

				my $gene_obj = new Gene_obj();

				$gene_obj->{TU_feat_name} = $gene_id;
				$gene_obj->{Model_feat_name} = "$trans_id";
				$gene_obj->{com_name} = "cufflinks $gene_id $trans_id";
				
				$gene_obj->{asmbl_id} = $scaff;
				
				$gene_obj->populate_gene_object($coords_href, $coords_href);


                ## encode gene and trans id in the ID and target fields for later extraction.  (probably a better way to do this!!)
				print $gene_obj->to_alignment_GFF3_format("GENE^$gene_id,TRANS^$trans_id", "GENE^$gene_id,TRANS^$trans_id", "Cufflinks");
				
				print "\n";
			}
		}
	}
    

	exit(0);
}

