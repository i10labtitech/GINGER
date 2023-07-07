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

######################################################################

use Getopt::Long;

my $opt_mapping = "";  # Preparation phase Mapping-based method only
my $opt_denovo = "";   # Preparation phase de novo-based method only
my $opt_homology = ""; # Preparation phase homology-based method only
my $opt_abinitio = ""; # Preparation phase ab initio-based method only
my $opt_phase0 = "";   # Merge phase 0 only
my $opt_phase1 = "";   # Merge phase 1 only
my $opt_phase1manual = ""; # Merge phase 1 only, need the threshold of gene 
my $opt_phase2 = "";   # Merge phase 2 only
my $opt_totalcds = ""; # total CDS minimum length in Merge phase 2 (i.e. threshold)
my $opt_summary = "";  # Merge phase summary only
my $opt_help = "";     # help

Getopt::Long::Configure("no_auto_abbrev");
GetOptions("mapping" => \$opt_mapping,   # Preparation phase mapping-based method only
           "denovo" => \$opt_denovo,     # Preparation phase de novo-based method only
           "homology" => \$opt_homology, # Preparation phase homology-based method only
           "abinitio" => \$opt_abinitio, # Preparation phase ab initio-based method only
           "phase0" => \$opt_phase0,     # Merge phase 0 only
           "phase1" => \$opt_phase1,     # Merge phase 1 only
           "phase1m=f" => \$opt_phase1manual, # Merge phase 1 only, need the threshold of gene 
           "phase2" => \$opt_phase2,     # Merge phase 2 only
           "totalcds=i" => \$opt_totalcds, # Total CDS minimum length in Merge phase 2 (i.e. threshold)
           "summary" => \$opt_summary,   # Merge summary only
           "help" => \$opt_help          # Help
    );

if ($opt_help) {
    print <<HELP;

Usage: .\/runGINGER.pl [Netflow configuration file] 

  --mapping        Preparation phase Mapping-based method only
  --denovo         Preparation phase de novo-based method only
  --homology       Preparation phase homology-based method only
  --abinitio       Preparation phase ab initio-based method only
  --phase0         Merge phase 0 only
  --phase1         Merge phase 1 only
  --phase1manual   Merge phase 1 only, need the threshold of gene 
  --phase2         Merge phase 2 only
  --totalcds       Total CDS minimum length in Merge phase 2 (i.e. threshold)
  --summary        Merge phase summary only
  --help           This help message

HELP
}

my $tc_threshold = 50;
if ($opt_totalcds =~ /\S/) {
    $tc_threshold = $opt_totalcds;
}

######################################################################

my $gingerDirInDockerImage   = "/GINGER_v1.0.1/pipeline";

######################################################################

my $inFile = $ARGV[0];
die "./runGINGER.pl [configuration file for user specific settings]"
    if ($#ARGV < 0);

my $nextflowConfigFileByUser = `realpath $inFile`;
chomp($nextflowConfigFileByUser);
die "No configuration file for user specific settings.\n($nextflowConfigFileByUser)"
    unless (-e $nextflowConfigFileByUser);

######################################################################

my $oprefix = "";
my $pdir    = "";
my $vOpt    = "";
{
    if ($nextflowConfigFileByUser =~ /^(\S+)\/[^\/]+$/) {
        $vOpt = "-v $1:$1";
    } else {
        die "No configuration file for user specific settings.\n($nextflowConfigFileByUser)";
    }
}
&loadConfigInfo($nextflowConfigFileByUser, \$oprefix, \$pdir, \$vOpt);

######################################################################

my $cmd = "";

if ($opt_mapping =~ /\S/) {
    $cmd = <<"GINGER";
mkdir -p $pdir
cd $pdir
cp $nextflowConfigFileByUser $pdir

docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; /nextflow -c $gingerDirInDockerImage/nextflow.config.base -c $nextflowConfigFileByUser -log $pdir/nextflow.log run $gingerDirInDockerImage/mapping.nf -work-dir $pdir/work" > $pdir/mapping.stdout 2> $pdir/mapping.stderr
GINGER

}elsif ($opt_denovo =~ /\S/) {
    $cmd = <<"GINGER";
mkdir -p $pdir
cd $pdir
cp $nextflowConfigFileByUser $pdir

docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; /nextflow -c $gingerDirInDockerImage/nextflow.config.base -c $nextflowConfigFileByUser -log $pdir/nextflow.log run $gingerDirInDockerImage/denovo.nf -work-dir $pdir/work" > $pdir/denovo.stdout 2> $pdir/denovo.stderr
GINGER

}elsif ($opt_homology =~ /\S/) {
    $cmd = <<"GINGER";
mkdir -p $pdir
cd $pdir
cp $nextflowConfigFileByUser $pdir

docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; /nextflow -c $gingerDirInDockerImage/nextflow.config.base -c $nextflowConfigFileByUser -log $pdir/nextflow.log run $gingerDirInDockerImage/homology.nf -work-dir $pdir/work" > $pdir/homology.stdout 2> $pdir/homology.stderr
GINGER

}elsif ($opt_abinitio =~ /\S/) {
    $cmd = <<"GINGER";
mkdir -p $pdir
cd $pdir
cp $nextflowConfigFileByUser $pdir

docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; /nextflow -c $gingerDirInDockerImage/nextflow.config.base -c $nextflowConfigFileByUser -log $pdir/nextflow.log run $gingerDirInDockerImage/abinitio.nf -work-dir $pdir/work" > $pdir/abinitio.stdout 2> $pdir/abinitio.stderr
GINGER

}elsif ($opt_phase0 =~ /\S/) {
    $cmd = <<"GINGER";
cd $pdir
cp $nextflowConfigFileByUser $pdir

docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; cd $pdir; bash $gingerDirInDockerImage/phase0.sh $nextflowConfigFileByUser" > $pdir/phase0.stdout 2> $pdir/phase0.stderr 
GINGER

}elsif ($opt_phase1 =~ /\S/) {
    die "No output directory ($pdir)." unless (-d $pdir);
    $cmd = <<"GINGER";
cd $pdir
cp $nextflowConfigFileByUser $pdir

docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; cd $pdir; bash $gingerDirInDockerImage/phase1.sh $nextflowConfigFileByUser" > $pdir/phase1.stdout 2> $pdir/phase1.stderr 
GINGER

}elsif ($opt_phase1manual =~ /\S/) {
    die "No output directory ($pdir)." unless (-d $pdir);
    $cmd = <<"GINGER";
cd $pdir
cp $nextflowConfigFileByUser $pdir

docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; cd $pdir; bash $gingerDirInDockerImage/phase1_manual.sh $nextflowConfigFileByUser $opt_phase1manual" > $pdir/phase1_manual.stdout 2> $pdir/phase1_manual.stderr 
GINGER

} elsif ($opt_phase2 =~ /\S/) {
    die "No output directory ($pdir)." unless (-d $pdir);
    $cmd = <<"GINGER";
cd $pdir
cp $nextflowConfigFileByUser $pdir

docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; cd $pdir; bash $gingerDirInDockerImage/phase2.sh $tc_threshold" > $pdir/phase2.stdout 2> $pdir/phase2.stderr
GINGER

} elsif ($opt_summary =~ /\S/) {
    die "No output directory ($pdir)." unless (-d $pdir);
    $cmd = <<"GINGER";
cd $pdir
cp $nextflowConfigFileByUser $pdir

docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; cd $pdir; bash $gingerDirInDockerImage/summary.sh $nextflowConfigFileByUser" > $pdir/summary.stdout 2> $pdir/summary.stderr 
GINGER

} else {
    $cmd = <<"GINGER";
mkdir -p $pdir
cd $pdir
cp $nextflowConfigFileByUser $pdir

docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; /nextflow -c $gingerDirInDockerImage/nextflow.config.base -c $nextflowConfigFileByUser -log $pdir/nextflow.log run $gingerDirInDockerImage/ginger.nf -work-dir $pdir/work" > $pdir/ginger.stdout 2> $pdir/ginger.stderr
#docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; /nextflow -c $gingerDirInDockerImage/nextflow.config.base -c $nextflowConfigFileByUser -log $pdir/nextflow.log run $gingerDirInDockerImage/ginger.nf -work-dir $pdir/work" > $pdir/ginger.stdout 2> $pdir/ginger.stderr

docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; cd $pdir; bash $gingerDirInDockerImage/phase0.sh $nextflowConfigFileByUser"  > $pdir/phase0.stdout 2> $pdir/phase0.stderr
docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; cd $pdir; bash $gingerDirInDockerImage/phase1.sh $nextflowConfigFileByUser"  > $pdir/phase1.stdout 2> $pdir/phase1.stderr
docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; cd $pdir; bash $gingerDirInDockerImage/phase2.sh $tc_threshold"              > $pdir/phase2.stdout 2> $pdir/phase2.stderr
docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; cd $pdir; bash $gingerDirInDockerImage/summary.sh $nextflowConfigFileByUser" > $pdir/summary.stdout 2> $pdir/summary.stderr
GINGER
}

######################################################################

system($cmd);

######################################################################

sub loadConfigInfo {
    my($nextflowConfigFileByUser,
       $oprefix, $pdir, $vOpt) = @_;

    my $scratch = "";
    my %vTarget = ();

    open(IN, $nextflowConfigFileByUser) 
        or die "can not open a file \"$nextflowConfigFileByUser\".";
    while (<IN>) {
        s/\#.*$//;
        if (/^\s*PDIR\s*\=\s*\"(\S+)\"/) {
            $$pdir = "$1";
            $vTarget->{$1} = 1;
        }
    }
    close(IN);
  
    my $flag = 0;
    open(IN, $nextflowConfigFileByUser) 
        or die "can not open a file \"$nextflowConfigFileByUser\".";
    while (<IN>) {
#        if (/^\s*PDIR_PREP\s*=\s*\"(\S+)\"/) {
#            my $fPath = &getFullpath($1, $$pdir);
#            $vTarget->{$fPath} = 1;
        if (/^\s*OPREFIX\s*\=\s*\"(\S+)\"/) {
            $$oprefix = $1;
        } elsif (/^\s*SCRATCH\s*\=\s*\"(\S+)\"/) {
            my $fPath = &getFullpath($1, $$pdir);
            $vTarget->{$fPath} = 1;
        } elsif (/^\s*INPUT_GENOME\s*\=\s*\"(\S+)\/[^\/]+\"/) {
            my $fPath = &getFullpath($1, $$pdir);
            $vTarget->{$fPath} = 1;
        } elsif (/^\s*INPUT_MASKEDGENOME\s*\=\s*\"(\S+)\/[^\/]+\"/) {
            my $fPath = &getFullpath($1, $$pdir);
            $vTarget->{$fPath} = 1;
        } elsif (/^\s*INPUT_REPOUT\s*\=\s*\"(\S+)\/[^\/]+\"/) {
            my $fPath = &getFullpath($1, $$pdir);
            $vTarget->{$fPath} = 1;
        } elsif (/^\s*INPUT_RNASEQR1\s*\=\s*\"(\S+)\/[^\/]+\"/) {
            my $fPath = &getFullpath($1, $$pdir);
            $vTarget->{$fPath} = 1;
        } elsif (/^\s*INPUT_RNASEQR2\s*\=\s*\"(\S+)\/[^\/]+\"/) {
            my $fPath = &getFullpath($1, $$pdir);
            $vTarget->{$fPath} = 1;
        } elsif (/^\s*RNASEQ_OTHER\d+\s*=\s*\"(\S+)\/[^\/]+\"/) {
            my $fPath = &getFullpath($1, $pdir);
            $vTarget->{$fPath} = 1;
        } elsif (/^\s*HOMOLOGY_OTHER\d+\s*=\s*\"(\S+)\/[^\/]+\"/) {
            my $fPath = &getFullpath($1, $pdir);
            $vTarget->{$fPath} = 1;
        } elsif (/^\s*ABINITIO_OTHER\d+\s*=\s*\"(\S+)\/[^\/]+\"/) {
            my $fPath = &getFullpath($1, $pdir);
            $vTarget->{$fPath} = 1;
        } elsif (/\"PROTEIN\"/) {
            $flag = 1;
        }
        if ($flag == 1) {
            if (/\"(\S+)\/[^\/]+\"/) {
                my $fPath = &getFullpath($1, $$pdir);
                $vTarget->{$fPath} = 1;
            }
        }
        if ((/\]/) && ($flag == 1)) {
            $flag = 0;
        }
    }
    close(IN);

    $$vOpt .= " -v /tmp:/tmp";
    foreach my $aPath (keys(%{$vTarget})) {
        $$vOpt .= " -v $aPath:$aPath";
    }
}

sub getFullpath {
    my($path, $pdir) = @_;
    
    $path =~ s/\$\{PDIR\}/$pdir\//;
    
    return $path;
}

__END__


