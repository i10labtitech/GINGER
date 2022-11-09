#!/usr/bin/perl

my $prevGeneId   = "";
my $prevContigId = "";
my @cds          = ();
my $modelNum     = 0;

my @elem = ();

while (<>) {
    last if (/^[\>\#]/);
    
    my $aLine = $_;
    chomp($aLine);
    my @elem = split(/\t/, $aLine);
    if ($elem[2] eq "CDS") {
        my $geneId = "";
        if ($elem[8] =~ /Parent=([^\;\s]+)/) {
            $geneId = $1;
        }
        if (($prevGeneId =~ /\S/) && ($prevGeneId ne $geneId)) {
            &outputCDS(\@cds);
            @cds = undef;
            @cds = ();
            push(@cds, \@elem);
            $prevContigId = $elem[0];
        } else {
            push(@cds, \@elem);
#            $prevContigId = $elem[0];
        }            
        $prevGeneId = $geneId;
    }
}

&outputCDS(\@cds);

sub outputCDS {
    my($c) = @_;
    
    if ($prevContigId ne $c->[0]->[0]) {
        print ">$c->[0]->[0]\n";
    }
    
    if ($c->[0]->[6] eq "+")  {
        @{$c} = sort{$a->[3] <=> $b->[3]} @{$c};  
    } elsif ($c->[0]->[6] eq "-")  {
        @{$c} = sort{$b->[3] <=> $a->[3]} @{$c};  
    } else {
        print STDERR "May be something wrong.\n";
    }
    
    if ($#{$c} == 0) {
        print "Esngl  ";
        if (($c->[0]->[6] eq "+") && ($c->[$i]->[3] <= $c->[$i]->[4])) {
            print $c->[$i]->[3], " ";
            print $c->[$i]->[4], " ";
        } elsif (($c->[0]->[6] eq "-") && ($c->[$i]->[3] <= $c->[$i]->[4])) {
            print $c->[$i]->[4], " ";
            print $c->[$i]->[3], " ";
        } else {
            print STDERR "May be something wrong.\n";
        }
        print "MODEL$modelNum\n";
    } else {
        for (my $i = 0; $i <= $#{$c}; $i++) {
            if ($i == 0) {
                print "Einit  ";
            } elsif ($i == $#{$c}) {
                print "Eterm  ";
            } else {
                print "Exon  ";
            }
            if (($c->[0]->[6] eq "+") && ($c->[$i]->[3] <= $c->[$i]->[4])) {
                print $c->[$i]->[3], " ";
                print $c->[$i]->[4], " ";
            } elsif (($c->[0]->[6] eq "-") && ($c->[$i]->[3] <= $c->[$i]->[4])) {
                print $c->[$i]->[4], " ";
                print $c->[$i]->[3], " ";
            } else {
                print STDERR "May be something wrong.\n";
            }
            print "MODEL$modelNum\n";
        }
    }
    $modelNum++;
    
=pod        
    } elsif ($c->[0]->[6] eq "-")  {
        for (my $i = $#{$c}; 0 <= $i; $i--) {
            if ($i == $#{$c}) {
                print "Einit  ";
            } elsif ($i == 0) {
                print "Eterm  ";
            } else {
                print "Exon  ";
            }
            print $c->[$i]->[4], " ";
            print $c->[$i]->[3], " ";
            print "MODEL$modelNum\n";
        }
        $modelNum++;
    } else {
        print STDERR "May be something wrong.\n";
    }
=cut

    $prevContigId = $c->[0]->[0];
} 
