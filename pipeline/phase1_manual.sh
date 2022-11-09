#!/bin/sh

LEMON_PATH=`realpath "$0"`
SCRIPT=`dirname $LEMON_PATH`/../util/merge_phase1

Gtool=${SCRIPT}/Grouping
GSub=${SCRIPT}/subgroup
NSub=${SCRIPT}/new_subgroup
Tool=${SCRIPT}/Searchalgo
Edit=${SCRIPT}/gff_editor
Polish=${SCRIPT}/initial_exon_polish
filter=${SCRIPT}/score_filtering.py
info=${SCRIPT}/info_annotate.py
ext_score=${SCRIPT}/extract_score.py
cutoff=${SCRIPT}/cutoff.py
sum=${SCRIPT}/sum_of_weight.py
intronLen=${SCRIPT}/intron_length.py
intronDist=${SCRIPT}/intron_distribution.py
refine=${SCRIPT}/refine_mrna_pos.py

if test $# -ne 5 ; then
    echo "---LEMON pipeline; phase1 (manual filtering mode) ---"
    echo "1: all.gff"
    echo "2: genome.fa"
    echo "3: config.ini"
    echo "4: output prefix"
    echo "5: score threshold"
else
    
    all=`readlink -f $1`
    Genome=`readlink -f $2`
    w=`python ${sum} $3`
    prefix=$4
    n=$5
    
    mkdir ${prefix}_phase1_result
    cd ${prefix}_phase1_result
    
    #------------Grouping
    echo "[1/3] Grouping"
    $Gtool -f $all
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
