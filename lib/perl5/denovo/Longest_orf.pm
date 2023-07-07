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

package main;
our $SEE; 

package Longest_orf;

use strict;
use warnings;
use Nuc_translator;
use Carp;

## if allow_partials is set, partial orfs are included in the analysis.  

# below used to be static, now instance vars.
#my $ALLOW_5PRIME_PARTIALS = 0; #allow for lacking start codon in logest orf.
#my $ALLOW_3PRIME_PARTIALS = 0; #allow for lacking stop codon in longest orf.
#my $FORWARD_STRAND = 1; #default set to true (analyze forward strand)
#my $REVERSE_STRAND = 1; #default set to true.
#my $ALLOW_NON_MET_STARTS = 0; #allow for non-methionine start codons.


sub new {
    shift;
	
    ## This object stores the longest ORF identified.
    my @stop_codons = &Nuc_translator::get_stop_codons(); # live call, depends on current genetic code.
    print "Stop codons in use: @stop_codons, set dynamically via current Nuc_translator settings.\n" if $SEE;
    unless (@stop_codons) {
	  confess "Fatal, no stop codons set";
	}

	my $obj = { pep_seq => undef,
				nt_seq => undef,
				length => undef, #length of nt_seq
				end5 => undef,
				end3 => undef,
				all_ORFS=>[], #container holds all ORFs found in order of decreasing length. Use orfs() method to retrieve them.
				stop_codons => [@stop_codons],
				
				## ORF settings
				ALLOW_5PRIME_PARTIALS => 0,
				ALLOW_3PRIME_PARTIALS => 0,
				FORWARD_STRAND => 1,
				REVERSE_STRAND => 1,
				ALLOW_NON_MET_STARTS => 0
					
				};
    bless ($obj);
    return ($obj);
}

## can include partial orfs at end of sequence.
sub allow_partials {
    my $self = shift;
    die unless (ref $self);
    $self->{ALLOW_5PRIME_PARTIALS} = 1;
    $self->{ALLOW_3PRIME_PARTIALS} = 1;
	
    if ($SEE) {
		print "Longest_orf: allowing both 5' and 3' partials.\n";
    }
}

sub allow_5prime_partials {
    my $self = shift;
    die unless (ref $self);
    $self->{ALLOW_5PRIME_PARTIALS} = 1;
    if ($SEE) {
		print "Longest_orf: allowing 5prime partials.\n";
    }
}

sub allow_3prime_partials {
    my $self = shift;
    die unless (ref $self);
    $self->{ALLOW_3PRIME_PARTIALS} = 1;
    if ($SEE) {
		print "Longest_orf: allowing 3prime partials\n";
    }
}

sub forward_strand_only {
    my $self = shift;
    die unless (ref $self);
    $self->{REVERSE_STRAND} = 0;
    if ($SEE) {
		print "Longest_orf: forward strand only.\n";
    }
	
}

sub reverse_strand_only {
    my $self = shift;
    die unless (ref $self);
    $self->{FORWARD_STRAND} = 0;
    if ($SEE) {
		print "Longest_orf: reverse strand only.\n";
    }
    
}

sub allow_non_met_starts {
    my $self = shift;
    $self->{ALLOW_NON_MET_STARTS} = 1;
    if ($SEE) {
		print "Longest_orf: allowing non Met start codons.\n";
    }
}


sub get_longest_orf {
    my $self = shift;
    my $input_sequence = shift;
    
    unless ($input_sequence) {
		print STDERR "I require a cDNA nucleotide sequence as my only parameter\n";
		return;
    }
    unless (length ($input_sequence) >= 3) {
		print STDERR "Sequence must code for at least a codon. Your seq_length is too short\n";
		return;
    }
    my @orfList = $self->capture_all_ORFs($input_sequence);
	# print "Found " . scalar @orfList . " orfs.\n";
	if (@orfList) {
		return ($orfList[0]); # longest ORF found is first in the sorted list.
	}
	else {
		## no ORFs found
		return (undef);
	}
}


sub capture_all_ORFs {
    
    my $self = shift;
    my $input_sequence = shift;

    unless ($input_sequence) {
		print STDERR "I require a cDNA nucleotide sequence as my only parameter\n";
		return;
    }
    unless (length ($input_sequence) >= 3) {
		print STDERR "Sequence must code for at least a codon. Your seq_length is too short\n";
		return;
    }
    
    $input_sequence = lc ($input_sequence);
    
    my (@starts, @stops, @orfs);

    if ($self->{FORWARD_STRAND}) {
		## analyse forward position
		@stops = $self->identify_putative_stops($input_sequence);
		@starts = $self->identify_putative_starts($input_sequence,\@stops);
		@orfs = $self->get_orfs (\@starts, \@stops, $input_sequence, '+');
    }
    
    if ($self->{REVERSE_STRAND}) {
		## reverse complement sequence and do again
		$input_sequence = &revcomp ($input_sequence);
		@stops = $self->identify_putative_stops($input_sequence);
		@starts = $self->identify_putative_starts($input_sequence, \@stops);
		push (@orfs,  $self->get_orfs (\@starts, \@stops, $input_sequence, '-'));
    }
	
    if (@orfs) {
		## set in order of decreasing length
		@orfs = reverse sort {$a->{length} <=> $b->{length}} @orfs;
		my $longest_orf = $orfs[0];
		my $seq = $longest_orf->{sequence};
    my $length = length($seq);
    my $protein = &translate_sequence($seq, 1);
		
		my $start = $longest_orf->{start};
		my $stop = $longest_orf->{stop};
		$self->{end5} = $start;  ## now coord is seq_based instead of array based.
		$self->{end3} = $stop;
		$self->{length} = $length;
		$self->{nt_seq} = $seq;
		$self->{pep_seq} = $protein;
		$self->{all_ORFS} = \@orfs;
	}

	return (@orfs);
}

sub check_for_fp_5prime_partials {
  my($self,$min_len,$orfs_ref,$pct_cutoff) = @_;
  my @orfs = @$orfs_ref;
  print "IN fp $#orfs\n" if($SEE);
  $pct_cutoff /= 100;
  for(my $i = 0; $i <= $#orfs; $i++) {
    use Data::Dumper;

    next if($orfs[$i]->{type} ne '5prime_partial');
    my $fl_start = index($orfs[$i]->{protein}, 'M');
    next if($fl_start == -1);

 
    my $fl_protein = substr($orfs[$i]->{protein}, $fl_start);
    print "partial with M\nnew prot: $fl_protein\n" if($SEE);
    print Dumper($orfs[$i]) if($SEE);

    next if(length($fl_protein) < $min_len);
    next if(length($fl_protein) < $pct_cutoff*length($orfs[$i]->{protein}));
    next if($i < $#orfs && length($orfs[$i+1]->{protein}) > length($fl_protein));
 
    print "seq start: $orfs[$i]->{start} new start: ".($orfs[$i]->{start} + $fl_start*3)." length: $orfs[$i]->{length}\n" if($SEE);
    $orfs[$i]->{sequence} = substr($orfs[$i]->{sequence},$fl_start*3);
    $orfs[$i]->{length} = length($orfs[$i]->{sequence});
    if($orfs[$i]->{orient} eq '+') {
      $orfs[$i]->{start} = $orfs[$i]->{start} +$fl_start*3;
    } else {
      $orfs[$i]->{start} = $orfs[$i]->{start} -$fl_start*3;
    }

    $orfs[$i]->{protein} = $fl_protein;
    $orfs[$i]->{type} = 'complete';
    
    print Dumper($orfs[$i]) if($SEE);
  }

}

sub orfs {
    my $self = shift;
    return (@{$self->{all_ORFS}});
}

#####################
# supporting methods
#####################

sub get_end5_end3 {
    my $self = shift;
    return ($self->{end5}, $self->{end3});
}

sub get_peptide_sequence {
    my $self = shift;
    return ($self->{pep_seq});
}

sub get_nucleotide_sequence {
    my $self = shift;
    return ($self->{nt_seq});
}


sub toString {
    my $self = shift;
    my ($end5, $end3) = $self->get_end5_end3();
    my $protein = $self->get_peptide_sequence();
    my $nt_seq = $self->get_nucleotide_sequence();
    my $ret_string = "Coords: $end5, $end3\n" 
	. "Protein: $protein\n"
	    . "Nucleotides: $nt_seq\n";
    return ($ret_string);
}


#################################

#Private methods:


sub get_orfs {
    my ($self, $starts_ref, $stops_ref, $seq, $direction) = @_;
    
	unless ($starts_ref && $stops_ref && $seq && $direction) {
		confess "Error, params not appropriate";
	}
	
	my %last_delete_pos = ( 0=>-1,
							1=>-1,
							2=>-1); #store position of last chosen stop codon in spec reading frame.
    my @orfs;
    my $seq_length = length ($seq);
	
    if ($SEE) {
		print "Potential Start codons: " . join (", ", @$starts_ref) . "\n";
		print "Potential Stop codons: " . join (", ", @$stops_ref) . "\n";
    }
    
    
    foreach my $start_pos (@{$starts_ref}) {
		my $start_pos_frame = $start_pos % 3;
		foreach my $stop_pos (@{$stops_ref}) {
		  # print "Comparing start: $start_pos to stop: $stop_pos, $direction\n";
		  if ( ($stop_pos > $start_pos)   && #end3 > end5
				 ( ($stop_pos - $start_pos) % 3 == 0) #must be in-frame
				 && ($start_pos > $last_delete_pos{$start_pos_frame})) #only count each stop once.
			{
				
				$last_delete_pos{$start_pos_frame} = $stop_pos;
				my ($start_pos_adj, $stop_pos_adj) = ( ($start_pos+1), ($stop_pos+1+2));
				#print "Startposadj: $start_pos_adj\tStopPosadj: $stop_pos_adj\n";
				# sequence based position rather than array-based
				
				my ($start, $stop) = ($direction eq '+') ? ($start_pos_adj, $stop_pos_adj) 
					: (&revcomp_coord($start_pos_adj, $seq_length), &revcomp_coord($stop_pos_adj, $seq_length));
				
				print "Retrieving ORF, Start: $start\tStop: $stop\n" if $SEE;
				my $orfSeq =  substr ($seq, $start_pos, ($stop_pos - $start_pos + 3)); #include the stop codon too.
				my $protein = &translate_sequence($orfSeq, 1);
				if ($protein =~ /\*.*\*/) {
				  confess "Fatal Error: Longest_orf: ORF returned which contains intervening stop(s): ($start-$stop, $direction\nProtein:\n$protein\nOf Nucleotide Seq:\n$seq\n";
				}

        my $got_start = 0;
        my $got_stop = 0;
        if ($protein =~ /^M/) {
            $got_start = 1;
        } 
        if ($protein =~ /\*$/) {
            $got_stop = 1;
        }
        
        my $prot_type = "";
        if ($got_start && $got_stop) {
            $prot_type = "complete";
        } elsif ($got_start) {
            $prot_type = "3prime_partial";
        } elsif ($got_stop) {
            $prot_type = "5prime_partial";
        } else {
            $prot_type = "internal";
        }
        print "Setting type $prot_type\n" if($SEE);

				my $orf = { sequence => $orfSeq,
							protein => $protein,
              type => $prot_type,
							start=>$start,
							stop=>$stop,
							length=>length($orfSeq),
							orient=>$direction
							};
				push (@orfs, $orf);
				last;
			}
		}
    }
    return (@orfs);
}


sub identify_putative_starts {
    my ($self, $seq, $stops_aref) = @_;
    my %starts;
    my %stops;
    foreach my $stop (@$stops_aref) {
		$stops{$stop} = 1;
    }
	
    if ($self->{ALLOW_5PRIME_PARTIALS} || $self->{ALLOW_NON_MET_STARTS}) {
		$starts{0} = 1 unless $stops{0};
		$starts{1} = 1 unless $stops{1};
		$starts{2} = 1 unless $stops{2};
    }
    
    if (! $self->{ALLOW_NON_MET_STARTS}) { #Look for ATG start codons.
		my $start_pos = index ($seq, "atg");
		while ($start_pos != -1) {
			$starts{$start_pos} = 1;
			#print "Start: $start_pos\n";
			$start_pos = index ($seq, "atg", ($start_pos + 1));
		}
    } else {
		# find all residues just subsequent to a stop codon, in-frame:
		foreach my $stop (@$stops_aref) {
			my $candidate_non_met_start = $stop +3;
			unless ($stops{$candidate_non_met_start}) {
				$starts{$candidate_non_met_start} = 1;
			}
		}
    }
    my @starts = sort {$a<=>$b} keys %starts;
    return (@starts);
}


sub identify_putative_stops {
  my ($self, $seq) = @_;
  my %stops;
  if ($self->{ALLOW_3PRIME_PARTIALS}) {
	## count terminal 3 nts as possible ORF terminators.
	my $seq_length = length ($seq);
	$stops{$seq_length} = 1;
	$seq_length--;
	$stops{$seq_length} = 1;
	$seq_length--;
	$stops{$seq_length} = 1;
  }
  my @stop_codons = @{$self->{stop_codons}};
  foreach my $stop_codon (@stop_codons) {
	$stop_codon = lc $stop_codon;
	print "Searching for stop codon: ($stop_codon).\n" if $SEE;
	my $stop_pos = index ($seq, $stop_codon);
	while ($stop_pos != -1) {
	  $stops{$stop_pos} = 1;
	  $stop_pos = index ($seq, $stop_codon, ($stop_pos + 1)); #include the stop codon too.
	}
  }
  my @stops = sort {$a<=>$b} keys %stops;
  return (@stops);
}


sub revcomp {
    my ($seq) = @_;
    my $reversed_seq = reverse ($seq);
    $reversed_seq =~ tr/ACGTacgtyrkm/TGCAtgcarymk/;
    return ($reversed_seq);
}


sub revcomp_coord {
    my ($coord, $seq_length) = @_;
    return ($seq_length - $coord + 1);
}




   
1;


