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

use File::Spec;

######################################################################

my $gingerDirInDockerImage   = "/GINGER_v1.0.1/pipeline";

######################################################################

my $validationDataFile = $ARGV[0];
my $predictionDataFile = $ARGV[1];
die "./runGINGER.pl [validation data file] [prediction data File]"
    if ($#ARGV < 1);

my $fullPathA = File::Spec->rel2abs($validationDataFile);
my $fullPathB = File::Spec->rel2abs($predictionDataFile);
my $vOpt = "";
if ($fullPathA =~ /^(\S+)\/[^\/]+\s*$/) {
    $vOpt .= " -v $1:$1";
} else {
    die "Somethig wrong maybe.";
} 
if ($fullPathB =~ /^(\S+)\/[^\/]+\s*$/) {
    $vOpt .= " -v $1:$1";
} else {
    die "Somethig wrong maybe.";
} 

my $cmd = <<"GINGER";
docker run $vOpt --rm i10labtitech/tools:GINGER_v1.0.1 /bin/bash -c ". /root/bashrc4exec ; cd $pdir; bash $gingerDirInDockerImage/evaluatePred.sh $fullPathA $fullPathB mRNA CDS mRNA CDS"
GINGER

system($cmd);

__END__


