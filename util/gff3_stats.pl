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

use FindBin;
use lib ("$FindBin::Bin/perlmod");

use my_bio;
use my_gff;

(@ARGV != 2) and die "usage: $0 gene.gff3 genome.fa\n";

gff3_tree_set_relations(\%gff3_tree, $ARGV[0]);

open($in, $ARGV[1]);
while (($name, $seq) = fasta_nonwhite_get($in)) {
	$ref{$name} = $seq;
}
close $in;

open IN, $ARGV[0];
while (<IN>) {
	chomp;
	($seq_name, $feature, $st, $ed, $strand, $attr) = (split(/\t/, $_))[0, 2, 3, 4, 6, 8];
	$id = gff3_attr_get_id($attr);

	if ($feature eq 'exon') {
		@mRNA_ids = gff3_tree_get_parent_id(\%gff3_tree, $id);
		for $mRNA_id (@mRNA_ids) {
			push(@{$mRNA{$mRNA_id}{exon_pos}}, $st);
			push(@{$mRNA{$mRNA_id}{exon_pos}}, $ed);

			++$mRNA{$mRNA_id}{num_exon};
			$mRNA{$mRNA_id}{len_exon} += $ed - $st + 1;
		}
	}
	elsif ($feature eq 'CDS') {
		@mRNA_ids = gff3_tree_get_parent_id(\%gff3_tree, $id);
		for $mRNA_id (@mRNA_ids) {
			push(@{$mRNA{$mRNA_id}{CDS_pos}}, $st);
			push(@{$mRNA{$mRNA_id}{CDS_pos}}, $ed);
		}
	}
	elsif ($feature eq 'mRNA' or $feature eq 'transcript') {
		$mRNA{$id}{seq_name} = $seq_name;
		$mRNA{$id}{strand} = $strand;
	}
}
close IN;


for (values %mRNA) {
	++$num_mRNA;
	$num_exon += $_->{num_exon};
	$total_len_exon += $_->{len_exon};
	++$num_single_exon if ($_->{num_exon} == 1);

	$seq = \$ref{$_->{seq_name}};

	@pos = sort{$a <=> $b} @{$_->{exon_pos}};
	$total_len_mRNA += $pos[$#pos] - $pos[0] + 1;
	for ($i = 1; $i < @pos - 2; $i += 2) {
		++$num_intron;
		$total_len_intron += ($pos[$i + 1] - $pos[$i] - 1);

		$s = substr($$seq, $pos[$i], 2) . substr($$seq, $pos[$i + 1] - 3, 2); 
		if ($_->{strand} eq '-') {
			$s = rev_comp($s);
		}
		++$splice_site{uc($s)};
	}

	@pos = sort{$a <=> $b} @{$_->{CDS_pos}};
	$s = '';
	for ($i = 0; $i < @pos - 1; $i += 2) {

		$s .= substr($$seq, $pos[$i] - 1, ($pos[$i + 1] - $pos[$i] + 1)); 
	}
	if ($_->{strand} eq '-') {
		$s = rev_comp($s);
	}
	++$num_CDS;
	$total_len_CDS += length($s);
}

print("#mRNAs\t$num_mRNA\n");
print("total_exon_intron_length\t$total_len_mRNA\n");
print("mean_exon_intron_length\t", $total_len_mRNA / $num_mRNA, "\n");
print("#exons\t$num_exon\n");
print("mean_#exons\t", $num_exon / $num_mRNA, "\n");
print("#single_exon_mRNAs\t$num_single_exon\n");
print("total_exon_length\t$total_len_exon\n");
print("mean_exon_length\t", $total_len_exon / $num_exon, "\n");
#print("\n");

print("#CDSs\t$num_CDS\n");
print("total_CDS_length\t", $total_len_CDS, "\n");
print("mean_CDS_length\t", $total_len_CDS / $num_CDS, "\n");
#print("\n");

print("#introns\t$num_intron\n");
print("total_intron_length\t$total_len_intron\n");
print("mean_intron_length\t", $total_len_intron / $num_intron, "\n");
print("%GTAG_splice_sites\t", $splice_site{GTAG} / $num_intron * 100, "\n");
