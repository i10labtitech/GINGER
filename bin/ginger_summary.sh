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

TMPFILE=$(mktemp)

NEXTFLOWCONFIG=$1

function rm_tmpfile {
    rm -f ${TMPFILE}.*
    rm -f ${TMPFILE}
}

GINGER_PATH=`dirname $(readlink -f $0)`
SCRIPT=${GINGER_PATH}/ginger-util

GenomePre=`perl -ne 'if (/INPUT_GENOME\s*\=\s*\"(\S+)\"/) {print "$1\n";}' ${NEXTFLOWCONFIG}`
Genome=`readlink -f ${GenomePre}`
prefix="ginger"
phase2output="${prefix}_phase2.gff" 

${SCRIPT}/final_reform ${phase2output} ${TMPFILE}.gff
#calcurate stats
${SCRIPT}/gff3_reformat.pl ${TMPFILE}.gff ${Genome} tmp GINGER > ${TMPFILE}.reformat.gff
${SCRIPT}/gff3_stats.pl ${TMPFILE}.gff ${Genome} > ${prefix}_stats.tsv
#make CDS seq. and protein seq.
python ${SCRIPT}/make_cds_from_gff.py ${phase2output} ${Genome} ${prefix}

trap "rm_tmpfile" EXIT INT PIPE TERM
