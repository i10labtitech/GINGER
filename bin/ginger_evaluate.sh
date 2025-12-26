#!/usr/bin/env bash

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

function abort
{
	echo "$@" 1>&2
	exit 1
}

ePath=`dirname $(readlink -f $0)`/ginger-util

# Checking input files
[ -f $1 ] || abort "$1 file not exist." # Validation data
[ -f $2 ] || abort "$2 file not exist." # Prediction data

# Checking arguments
if [ $# -ne 6 ]; then
    echo "" 1>&2
	echo "Error! evaluatioPred.sh needs 6 arguments!" 1>&2
    echo "" 1>&2
	echo "  Arg 1: GFF file containing Validation data" 1>&2
	echo "  Arg 2: GFF file containing Prediction data" 1>&2
    echo "" 1>&2
    echo "  Followings are type in GFF (column 2)" 1>&2
    echo "  <Related args to the file specified by arg 1>" 1>&2
	echo "  Arg 3: Delimiters for gene descriptions in first file (may be \"mRNA\" or \"gene\") " 1>&2
	echo "  Arg 4: Evaluation target in first file (may be \"CDS\") " 1>&2
    echo "" 1>&2
    echo "  <Related args to the file specified by arg 2>" 1>&2
	echo "  Arg 5: Delimiters for gene descriptions in first file (may be \"mRNA\" or \"gene\") " 1>&2
	echo "  Arg 6: Evaluation target in first file (may be \"CDS\") " 1>&2
    echo "" 1>&2
	echo " ex) ./evaluatePred.sh annotation.gff annotation2.gff mRNA CDS mRNA CDS" 1>&2
	echo ""
	echo " You should delete line(type==\"exon\")" 1>&2
	exit 1
fi

tmpFile1=`mktemp`
tmpFile2=`mktemp`

$ePath/preevaluate $1 $tmpFile1 $3 $4
$ePath/preevaluate $2 $tmpFile2 $5 $6
$ePath/evaluate $tmpFile1 $tmpFile2

rm $tmpFile1
rm $tmpFile2
