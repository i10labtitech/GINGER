#!/usr/bin/perl

use FindBin;
use lib ("$FindBin::Bin/perlmod");

use my_bio;
use my_gff;

(@ARGV != 4) and die "usage: $0 genes.gff3 genome.fa gene_name_prefix source_name";

$gene_prefix = $ARGV[2];
$source_name = $ARGV[3];

open($in, $ARGV[1]);
while (($name, $seq) = fasta_nonwhite_get($in)) {
	push(@seq_name_order, $name);
	$seq_name2id{$name} = $n_seq;
	++$n_seq;
}
close $in;

gff3_tree_set_relations(\%gff3_tree, $ARGV[0]);

open($in, $ARGV[0]);
while (chomp($l = <$in>)) {
	($seq_name, $feature, $st, $ed, $strand, $attr) = (split(/\t/, $l))[0, 2, 3, 4, 6, 8];
	$id = gff3_attr_get_id($attr);
	if ($feature eq 'gene') {
		$genes{$id}{seq_id} = $seq_name2id{$seq_name};
		$genes{$id}{start} = $st;
		$genes{$id}{end} = $ed;
		$genes{$id}{strand} = $strand;
	}
	elsif ($feature eq 'mRNA') {
		$gene_id = gff3_tree_get_root_id(\%gff3_tree, $id);
		push(@{$genes{$gene_id}{children}}, $id);
		$mRNAs{$id}{start} = $st;
		$mRNAs{$id}{end} = $ed;
	}
	elsif ($feature eq 'exon' or $feature eq 'CDS') {
		@mRNA_ids = gff3_tree_get_parent_id(\%gff3_tree, $id);
		for $mRNA_id (@mRNA_ids) {
			push(@{$mRNAs{$mRNA_id}{$feature}{children}}, $id . $st . $ed);
			$elements{$feature}{$id . $st . $ed}{start} = $st;
			$elements{$feature}{$id . $st . $ed}{end} = $ed;
		}
	}
}
close $in;

print "##gff-version 3\n";

$gene_id = 0;
for $gene (sort{($a->{seq_id} <=> $b->{seq_id}) or ($a->{start} <=> $b->{start})} values %genes) {
	++$gene_id;
	$gene_name = sprintf("%s%07d", $gene_prefix, $gene_id);
	print(
		join("\t", ($seq_name_order[$gene->{seq_id}], $source_name, 'gene', $gene->{start}, $gene->{end}, '.', $gene->{strand}, '.', sprintf("ID=%s", $gene_name)))
		, "\n")
	;

	$mRNA_id = 0;
	for $mRNA (map{$mRNAs{$_}} @{$gene->{children}}) {
		++$mRNA_id;
		$mRNA_name = sprintf("%s.mrna%d", $gene_name, $mRNA_id);
		print(
			join("\t", ($seq_name_order[$gene->{seq_id}], $source_name, 'mRNA', $mRNA->{start}, $mRNA->{end}, '.', $gene->{strand}, '.', sprintf("ID=%s;Parent=%s", $mRNA_name, $gene_name)))
			, "\n")
		;

		for $element_sp ('exon', 'CDS') {
			@tmp_elements = map{$elements{$element_sp}{$_}} @{$mRNA->{$element_sp}{children}};
			if ($gene->{strand} eq '+') {
				@tmp_elements = sort{$a->{start} <=> $b->{start}} @tmp_elements;
			}
			else {
				@tmp_elements = sort{$b->{end} <=> $a->{end}} @tmp_elements;
			}
				
			$element_id = 1;
			for $element (@tmp_elements) {
				$element_name = sprintf("%s.%s%d", $mRNA_name, lc($element_sp), $element_id);
				print(
					join("\t", ($seq_name_order[$gene->{seq_id}], $source_name, $element_sp, $element->{start}, $element->{end}, '.', $gene->{strand}, '.', sprintf("ID=%s;Parent=%s", $element_name, $mRNA_name)))
					, "\n")
				;

				if ($element_sp ne 'CDS') {
					++$element_id;
				}
			}
		}
	}
}
