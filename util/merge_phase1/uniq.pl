#!/usr/bin/perl

my $key = "";
my $data = "";
my %exist = ();

while (<>) {
    my $aLine = $_;
    chomp($aLine);
    my @elem = split(/\t/, $aLine);
    
    if ($elem[2] eq "gene") {
        if ($exist{$key} == 1) {
        } else {
            print $data;
            $exist{$key} = 1;
        }
        $key   = "$elem[0];$elem[2];$elem[3];$elem[4];$elem[5];$elem[6];$elem[7]";
        $data  = join("\t", @elem) . "\n";
    } else {
        $key  .= "$elem[0];$elem[2];$elem[3];$elem[4];$elem[5];$elem[6];$elem[7]";
        $data .= join("\t", @elem) . "\n";
    }
}

if ($exist{$key} == 1) {
} else {
    print $data;
}
        
__END__
