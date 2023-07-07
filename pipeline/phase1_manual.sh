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

NEXTFLOWCONFIG=$1
THRESHOLD=$2

GINGER_PATH=`dirname $(readlink -f $0)`
SCRIPT=${GINGER_PATH}/../util/merge_phase1

Gtool=${SCRIPT}/Grouping
GSub=${SCRIPT}/subgroup
NSub=${SCRIPT}/new_subgroup
Tool=${SCRIPT}/Searchalgo
Edit=${SCRIPT}/gff_editor
Polish=${SCRIPT}/initial_exon_polish
filter=${SCRIPT}/score_filtering.py
info=${SCRIPT}/info_annotate.py
ext_score=${SCRIPT}/extract_score.py
cutoff=${SCRIPT}/cutoff.pl
sum=${SCRIPT}/sum_of_weight.pl
intronLen=${SCRIPT}/intron_length.py
intronDist=${SCRIPT}/intron_distribution.py
refine=${SCRIPT}/refine_mrna_pos.py

if test $# -ne 2 ; then
    echo "---GINGER pipeline; phase1 (manual filtering mode) ---"
    echo "1: Nextflow configuration file"
    echo "2: score threshold"
else
    
    all=`readlink -f ginger_all.gff`
    GenomePre=`perl -ne 'if (/INPUT_GENOME\s*\=\s*\"(\S+)\"/) {print "$1\n";}' ${NEXTFLOWCONFIG}`
    Genome=`readlink -f ${GenomePre}`
    w=`perl ${sum} ${NEXTFLOWCONFIG}`
    prefix="ginger"
    n=${THRESHOLD}
    
    mkdir ${prefix}_phase1_result
    cd ${prefix}_phase1_result
    
    #------------Grouping
    echo "[1/3] Grouping"
    echo "$Gtool -f $all"
    $Gtool -f $all > /dev/null
    python ${intronLen} $all mRNA CDS > ${prefix}_intron.txt
    intron_limit=`python ${intronDist} ${prefix}_intron.txt | grep "predicted" | awk '{printf "%d\n", $3}'`
    echo "intron_limit"
    echo $intron_limit
    $GSub Group.gff 0.25 $w > Group_sub_pre.gff
    python ${refine} Group_sub_pre.gff > Group_sub_pre_refine.gff
    $NSub Group_sub_pre_refine.gff $intron_limit > Group_sub.gff
    
    #------------Searching
    echo "[2/3] DP search"
    $Tool -fa $Genome -gff Group_sub.gff -o ${prefix}_raw.gff
    
    #------------filtering low score mRNAs
    echo "[3/3] Filtering"
    #python $info ${prefix}_raw.gff Group_sub.gff ${prefix}_score.gff
    #python $ext_score ${prefix}_score.gff > ${prefix}_score.txt
    #n=`python $cutoff ${prefix}_score.txt`
    echo "score threshold (specified by the user) = ${n}"
    python $filter ${prefix}_raw.gff Group_sub.gff $n ${prefix}_filter.gff
    
    #------------filtering single exon genes
    
    $Edit -f ${prefix}_filter.gff -num 1 -v
    mv filtered.gff ${prefix}_filter_nosig.gff
    
    #------------polishing initial exons
    
    $Polish -i ${prefix}_filter_nosig.gff -c $all -f $Genome -o ${prefix}_filter_nosig_polish.gff
    
    #------------Combining
    
    cd ../
    ln -s `pwd`/${prefix}_phase1_result/${prefix}_filter_nosig_polish.gff ./${prefix}_phase1.gff
    
fi
echo "finished"
