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

GINGER_PATH=`dirname $(readlink -f $0)`
SCRIPT=${GINGER_PATH}/ginger-util
prefix="ginger"
phase1output="${prefix}_phase1.gff" 
phase1group="${prefix}_phase1_result/Group.gff" 

if test $# -ne 1 ; then
    echo "--- GINGER pipeline; phase2 ---"
    echo "1: total CDS minimum length"
else
    echo "Total CDS minimum length:" $1

    python ${SCRIPT}/multi_exon_trim_from_group.py ${phase1group} > ${prefix}_tmp_sin.gff
    ${SCRIPT}/grouping_v1 -f ${prefix}_tmp_sin.gff -o ${prefix}_tmp_gro.gff > /dev/null
    python ${SCRIPT}/single_group_filtering.py ${prefix}_tmp_gro.gff 50 > ${prefix}_tmp_fin.gff
    ${SCRIPT}/geneadd ${phase1output} ${prefix}_tmp_fin.gff > ${prefix}_tmp_phase2_prefilter.gff
    python ${SCRIPT}/cdslen_filter.py ${prefix}_tmp_phase2_prefilter.gff $1 > ${prefix}_phase2.gff
    rm ${prefix}_tmp_*.gff
fi
echo "finished"
