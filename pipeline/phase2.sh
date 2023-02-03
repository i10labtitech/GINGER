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

LEMON_PATH=`realpath "$0"`
SCRIPT=`dirname $LEMON_PATH`/../util/merge_phase2

if test $# -ne 4 ; then
    echo "--- LEMON pipeline; phase2 ---"
    echo "1: phase1.gff"
    echo "2: Group.gff"
    echo "3: total CDS minimum length"
    echo "4: output prefix"
else
    python ${SCRIPT}/multi_exon_trim_from_group.py $2 > ${4}_tmp_sin.gff
    ${SCRIPT}/Grouping -f ${4}_tmp_sin.gff -o ${4}_tmp_gro.gff > /dev/null
    python ${SCRIPT}/single_group_filtering.py ${4}_tmp_gro.gff 50 > ${4}_tmp_fin.gff
    ${SCRIPT}/geneadd $1 ${4}_tmp_fin.gff > ${4}_tmp_phase2_prefilter.gff
    python ${SCRIPT}/cdslen_filter.py ${4}_tmp_phase2_prefilter.gff $3 > ${4}_phase2.gff
    rm ${4}_tmp_*.gff
fi
echo "finished"
