#!/bin/sh

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

TMPFILE=$(mktemp)

function rm_tmpfile {
    rm -f ${TMPFILE}.*
    rm -f ${TMPFILE}
}

LEMON_PATH=`realpath "$0"`
SCRIPT=`dirname $LEMON_PATH`/../util/summary

if test $# -ne 3 ; then
    echo "1. phase2.gff"
    echo "2. genome.fa"
    echo "3. output prefix"
else
    ${SCRIPT}/final_reform $1 ${TMPFILE}.gff
    #calcurate stats
    ${SCRIPT}/gff3_reformat.pl ${TMPFILE}.gff $2 tmp LEMON > ${TMPFILE}.reformat.gff
    ${SCRIPT}/gff3_stats.pl ${TMPFILE}.gff $2 > ${3}_stats.tsv
    #make CDS seq. and protein seq.
    python ${SCRIPT}/make_cds_from_gff.py $1 $2 $3
fi

trap "rm_tmpfile" EXIT INT PIPE TERM
