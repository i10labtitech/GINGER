#!/usr/local/bin/perl -w

=pod
Copyright (c) 2008, Brian Haas  
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
=cut

# lightweight fasta reader capabilities:
package Fasta_reader;

use strict;

sub new {
    my ($packagename, $fastaFile) = @_;

	## note: fastaFile can be a filename or an IO::Handle
	

    my $self = { fastaFile => undef,,
				 fileHandle => undef };

    bless ($self, $packagename);
    
    ## create filehandle
    my $filehandle = undef;
    
	if (ref $fastaFile eq 'IO::Handle') {
		$filehandle = $fastaFile;
	}
	else {
		
		open ($filehandle, $fastaFile) or die "Error: Couldn't open $fastaFile\n";
		$self->{fastaFile} = $fastaFile;
	}
	
	$self->{fileHandle} = $filehandle;

    return ($self);
}



#### next() fetches next Sequence object.
sub next {
    my $self = shift;
    my $orig_record_sep = $/;
    $/="\n>";
    my $filehandle = $self->{fileHandle};
    my $next_text_input = <$filehandle>;
    
	if (defined($next_text_input) && $next_text_input !~ /\w/) {
		## must have been some whitespace at start of fasta file, before first entry.
		## try again:
		$next_text_input = <$filehandle>;
	}
	
	my $seqobj = undef;
    
	if ($next_text_input) {
		$next_text_input =~ s/^>|>$//g; #remove trailing > char.
		$next_text_input =~ tr/\t\n\000-\037\177-\377/\t\n/d; #remove cntrl chars
		my ($header, @seqlines) = split (/\n/, $next_text_input);
		my $sequence = join ("", @seqlines);
		$sequence =~ s/\s//g;
		
		$seqobj = Sequence->new($header, $sequence);
    }
    
    $/ = $orig_record_sep; #reset the record separator to original setting.
    
    return ($seqobj); #returns null if not instantiated.
}


#### finish() closes the open filehandle to the query database.
sub finish {
    my $self = shift;
    my $filehandle = $self->{fileHandle};
    close $filehandle;
    $self->{fileHandle} = undef;
}

####
sub retrieve_all_seqs_hash {
	my $self = shift;

	my %acc_to_seq;
	
	while (my $seq_obj = $self->next()) {
		my $acc = $seq_obj->get_accession();
		my $sequence = $seq_obj->get_sequence();

		$acc_to_seq{$acc} = $sequence;
	}

	return(%acc_to_seq);
}



##############################################
package Sequence;
use strict;

sub new {
    my ($packagename, $header, $sequence) = @_;
    
    ## extract an accession from the header:
    my ($acc, $rest) = split (/\s+/, $header, 2);
        
    my $self = { accession => $acc,
		 header => $header,
		 sequence => $sequence,
		 filename => undef };
    bless ($self, $packagename);
    return ($self);
}

####
sub get_accession {
    my $self = shift;
    return ($self->{accession});
}

####
sub get_header {
    my $self = shift;
    return ($self->{header});
}

####
sub get_sequence {
    my $self = shift;
    return ($self->{sequence});
}

#### 
sub get_FASTA_format {
    my $self = shift;
    my $header = $self->get_header();
    my $sequence = $self->get_sequence();
    $sequence =~ s/(\S{60})/$1\n/g;
    my $fasta_entry = ">$header\n$sequence\n";
    return ($fasta_entry);
}


####
sub write_fasta_file {
    my $self = shift;
    my $filename = shift;

    my ($accession, $header, $sequence) = ($self->{accession}, $self->{header}, $self->{sequence});
    
	my $fasta_entry = $self->get_FASTA_format();
	
    my $tempfile;
    if ($filename) {
		$tempfile = $filename;
    } else {
		my $acc = $accession;
		$acc =~ s/\W/_/g;
		$tempfile = "$acc.fasta";
    }
    
    open (TMP, ">$tempfile") or die "ERROR! Couldn't write a temporary file in current directory.\n";
    print TMP $fasta_entry;
    close TMP;
    return ($tempfile);
}

1; #EOM


