#!/usr/bin/perl

my $score = "";
my $size  = "";

while (<>) {
    my $aLine = $_;
    chomp($aLine);
    my @elem = split(/\t/, $aLine);

    if ($elem[2] eq "mRNA") {
        if ($size =~ /\S/) {
            print "$score\t$size\t$ne\n";
        }
        if ($elem[8] =~ /normalized_score=([^\;]+)/) {
            $score = $1;
            $size = $elem[4] - $elem[3] - 1;
        } else {
            die;
        }
        $ne = 0;
    } elsif ($elem[2] eq "CDS") {
        $ne++;
    } 
}

print "$score\t$size\t$ne\n";

__END__
NC_000001.11    candidate       mRNA    11774   13467   268     +       .       ID=group_num_1_1.mrna1;Parent=group_num_1_1;normalized_score=0.16;sup
