package my_bio;

use Exporter;

@ISA = qw(Exporter);
@EXPORT = qw/rev_comp print_seq fasta_get fasta_nonwhite_get nt2aa nt2aa_yeast_mt/;

sub rev_comp
{
	my $seq = shift;
	my $rev;

	$rev = reverse $seq;
	$rev =~ tr/ATGCatgc/TACGtacg/;

	return $rev;
}

sub print_seq
{
	my $seq = shift;
	my $line_len = shift;
	my $seq_len = length $seq;

	unless ($line_len) {
		$line_len = 80;
	}

	for (my $i = 0; $i < $seq_len; $i += $line_len) {
		print(substr($seq, $i, $line_len) . "\n");
	}
}

sub print_seq_file_handle
{
	my $out = shift;
	my $seq = shift;
	my $line_len = shift;
	my $seq_len = length $seq;

	unless ($line_len) {
		$line_len = 80;
	}

	for (my $i = 0; $i < $seq_len; $i += $line_len) {
		print($out substr($seq, $i, $line_len) . "\n");
	}
}

sub fasta_get
{
	my $in = shift;
	my($l, $name, $seq);

	while ($l = <$in>) {
		if (substr($l, 0, 1) eq '>') {
			chomp $l;
			last;
		}
	}
	if ($l eq '') {
		return ();
	}

	$name = substr($l, 1);

	while ($l = <$in>) {
		if (substr($l, 0, 1) eq '>') {
			seek($in, -length($l), 1);
			return ($name, $seq);
		}
		chomp $l;
		$seq .= $l 
	}
	return ($name, $seq);
}

sub fasta_nonwhite_get
{
	my $in = shift;
	my($l, $name, $seq);

	while ($l = <$in>) {
		if ($l =~ /^>(\S*)/) {
			$name = $1;
			last;
		}
	}
	if ($l eq '') {
		return ();
	}

	while ($l = <$in>) {
		if (substr($l, 0, 1) eq '>') {
			seek($in, -length($l), 1);
			return ($name, $seq);
		}
		chomp $l;
		$seq .= $l 
	}
	return ($name, $seq);
}

{
	%codon_table = (
	'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TTA' => 'L',
	'TTG' => 'L', 'CTT' => 'L', 'CTC' => 'L', 'CTA' => 'L', 'CTG' => 'L',
	'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R',
	'AGG' => 'R', 'AAA' => 'K', 'AAG' => 'K', 'AAT' => 'N', 'AAC' => 'N',
	'ATG' => 'M', 'GAT' => 'D', 'GAC' => 'D', 'TTT' => 'F', 'TTC' => 'F',
	'TGT' => 'C', 'TGC' => 'C', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P',
	'CCG' => 'P', 'CAA' => 'Q', 'CAG' => 'Q', 'TCT' => 'S', 'TCC' => 'S',
	'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S', 'GAA' => 'E',
	'GAG' => 'E', 'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
	'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G', 'TGG' => 'W',
	'CAT' => 'H', 'CAC' => 'H', 'TAT' => 'Y', 'TAC' => 'Y', 'ATT' => 'I',
	'ATC' => 'I', 'ATA' => 'I', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V',
	'GTG' => 'V', 'TAA' => '*', 'TAG' => '*', 'TGA' => '*',
	);

	%codon_table_yeast_mt = (
	'GCT' => 'A', 'GCC' => 'A', 'GCA' => 'A', 'GCG' => 'A', 'TTA' => 'L',
	'TTG' => 'L', 'CTT' => 'T', 'CTC' => 'T', 'CTA' => 'T', 'CTG' => 'T',
	'CGT' => 'R', 'CGC' => 'R', 'CGA' => 'R', 'CGG' => 'R', 'AGA' => 'R',
	'AGG' => 'R', 'AAA' => 'K', 'AAG' => 'K', 'AAT' => 'N', 'AAC' => 'N',
	'ATG' => 'M', 'GAT' => 'D', 'GAC' => 'D', 'TTT' => 'F', 'TTC' => 'F',
	'TGT' => 'C', 'TGC' => 'C', 'CCT' => 'P', 'CCC' => 'P', 'CCA' => 'P',
	'CCG' => 'P', 'CAA' => 'Q', 'CAG' => 'Q', 'TCT' => 'S', 'TCC' => 'S',
	'TCA' => 'S', 'TCG' => 'S', 'AGT' => 'S', 'AGC' => 'S', 'GAA' => 'E',
	'GAG' => 'E', 'ACT' => 'T', 'ACC' => 'T', 'ACA' => 'T', 'ACG' => 'T',
	'GGT' => 'G', 'GGC' => 'G', 'GGA' => 'G', 'GGG' => 'G', 'TGG' => 'W',
	'CAT' => 'H', 'CAC' => 'H', 'TAT' => 'Y', 'TAC' => 'Y', 'ATT' => 'I',
	'ATC' => 'I', 'ATA' => 'M', 'GTT' => 'V', 'GTC' => 'V', 'GTA' => 'V',
	'GTG' => 'V', 'TAA' => '*', 'TAG' => '*', 'TGA' => 'W',
	);

	sub nt2aa
	{
		my $nt = shift;
		my($i, $aa);

		for ($i = 0; $i < length($nt); $i += 3) {
			my $c = $codon_table{uc(substr($nt, $i, 3))};
			if ($c) {
				$aa .= $c;
			}
			else {
				$aa .= 'X';
			}
		}

		return $aa;
	}

	sub nt2aa_yeast_mt
	{
		my $nt = shift;
		my($i, $aa);

		for ($i = 0; $i < length($nt); $i += 3) {
			my $c = $codon_table_yeast_mt{uc(substr($nt, $i, 3))};
			if ($c) {
				$aa .= $c;
			}
			else {
				$aa .= 'X';
			}
		}

		return $aa;
	}
}

1;
