#!/bin/sh

#homology.shを動かすのに必要なscriptの位置
SCRIPT=/data2/fumiya/hetero_annotate/Lemon_20191029/Prep/homology/homologyscript
SEQKIT=seqkit

if test $# -ne 3 ; then
    echo "./homology_filter.sh all_spaln_result.gff genome.fa [output file]"
else
    outputspaln=$3
    $SEQKIT fx2tab -nl $2 | awk '{print $1"\t"$2}' > tmp_genome_stats
    python ${SCRIPT}/outer_trim.py $1 tmp_genome_stats > ${1}.trim
    sed -e 's/cds/CDS/g' ${1}.trim > ${1}.rename
    ${SCRIPT}/gff_2_proteinfasta -i ${1}.rename -g $2 -o prefilter01.fa
    g++ ${SCRIPT}/flameshiftgrep.cpp -std=c++0x -O3 -o flameshiftfilter
    ./flameshiftfilter prefilter01.fa flameshiftlist.txt
    sed -e 's/mRNA/gene/g' flameshiftlist.txt > flameshiftlist.rename
    python ${SCRIPT}/gff_trim.py ${1}.rename flameshiftlist.rename > ${outputspaln}.rename
    sed -e 's/CDS/cds/g' ${outputspaln}.rename > $outputspaln
    rm flameshiftfilter prefilter01.fa flameshiftlist.txt
    rm flameshiftlist.rename
    rm ${1}.trim
    rm ${1}.rename
    rm ${outputspaln}.rename
    rm tmp_genome_stats
fi

