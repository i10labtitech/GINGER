#!/usr/bin/perl

my @val        = ();
my @hist       = ();
my @histSmooth = ();
my $z          = 2;
my $step       = 0.05;
my $min        = undef;
my $histMax    = 0;
my $firstMinScore = undef;
    
while (<>) {
    if (/^(\S+)/) {
        my $aVal = log($1) / log(10.0);
        push(@val, $aVal);
        if (defined($min)) {
            $min = $aVal if ($aVal < $min);
        } else {
            $min = $aVal;
        }
    }
}

for (my $i = 0; $i <= $#val; $i++) {
    my $j = int(($val[$i] - $min) / $step);
    $hist[$j]++;
}

for (my $i = 0; $i <= $#hist; $i++) {
    $histMax = $hist[$i] if ($histMax < $hist[$i]);
}

for (my $i = 0; $i <= $#hist; $i++) {
    my $x0 = $i - $z;
    my $x1 = $i + $z;
    $x0 = 0 if ($x0 < 0);
    $x1 = $#hist if ($#hist < $x1);
    my $n = $x1 - $x0 + 1;
    my $s = 0;
    for (my $j = $x0; $j <= $x1; $j++) {
        $s += $hist[$j];
    }
    $histSmooth[$i] = $s / $n;
#    print "$i\t$histSomooth[$i]\n";
}

my $prevHistVal = undef;
my $prevSign    = undef;
for (my $i = $#hist; 0 <= $i; $i--) {
    if (defined($prevHistVal)) {
        my $sign = undef;
        if ($histMax / 100.0 < $histSmooth[$i]) {
            if ($prevHistVal - $histSmooth[$i] < 0) {
#                print "-";
                $sign = -1;
            } else {
#                print "+";
                $sign = 1;
            }
#            print "\n";
        }
        if (defined($sign)) {
            if (defined($prevSign)) {
                if (($prevSign == 1) && ($sign == -1)) {
                    my $aScore = 10.0 ** (($i + 1) * $step + $step / 2.0 + $min);
#                    print "*\t";
#                    printf("%0.4f", 10.0 ** (($i + 1) * $step + $min));
#                    print "\t";
#                    printf("%0.4f", $aScore);
#                    print "\n";
                    $firstMinScore = $aScore unless (defined($firstMinScore));
                }
            }
            $prevSign = $sign;
        }
    }
#    printf("%0.4f", $i * $step + $min);
#    print "\t";
#    printf("%0.4f", 10.0 ** ($i * $step + $min));
#    print "\t";
#    print $hist[$i];
#    print "\t";
#    print $histSmooth[$i];
#    print "\n";
    $prevHistVal = $histSmooth[$i];
}

printf("%0.2f\n", $firstMinScore);
