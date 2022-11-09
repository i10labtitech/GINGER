#!/bin/sh

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
