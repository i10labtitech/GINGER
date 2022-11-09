#!/bin/sh

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
