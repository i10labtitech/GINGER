#!/bin/sh

TMPFILE=$(mktemp)
#TMPFILE="tmp"

function rm_tmpfile {
    rm -f ${TMPFILE}.*
    rm -f ${TMPFILE}
}

LEMON_PATH=`realpath "$0"`
SCRIPT=`dirname $LEMON_PATH`/../util/merge_phase0

if test $# -ne 2 ; then
    echo "--- LEMON pipeline; phase0 ---"
    echo "1.config.ini"
    echo "2.prefix"
else
    #check file existence
    [ -f $MAPPING_GTF ] || abort "$MAPPING_GTF file not exist"
    [ -f $MAPPING_TO_MERGE ] || abort "$MAPPING_TO_MERGE file not exist"
    [ -f $GMAP_RESULT ] || abort "$GMAP_RESULT file not exist"
    [ -f $DENOVO_RESULT ] || abort "$DENOVO_RESULT file not exist"
    [ -f $HOMOLOGY_RESULT ] || abort "$HOMOLOGY_RESULT file not exist"
    [ -f $AUGUSTUS_RESULT ] || abort "$AUGUSTUS_RESULT file not exist"
    [ -f $SNAP_RESULT ] || abort "$SNAP_RESULT file not exist"
    
    #import config
    source ./$1
    python ${SCRIPT}/make_weight.py $MAPPING_WEIGHT $DENOVO_WEIGHT $HOMOLOGY_WEIGHT $AUGUSTUS_WEIGHT $SNAP_WEIGHT > ${TMPFILE}.weight
    echo "making temporary dir."
    echo ${TMPFILE}

    #mapping
    echo "RNA-seq mapping based"
    python ${SCRIPT}/cov_from_gtf.py $MAPPING_GTF > ${TMPFILE}.cov_stats
    python ${SCRIPT}/binning_mapping.py $MAPPING_TO_MERGE ${TMPFILE}.cov_stats > ${TMPFILE}.merge_cov
    python ${SCRIPT}/scoring_mapping.py ${TMPFILE}.merge_cov ${TMPFILE}.cov_stats > ${TMPFILE}.mapping
    ${SCRIPT}/RNA-seq_reform ${TMPFILE}.mapping ${TMPFILE}.mapping.reform
    ${SCRIPT}/Row2_rename ${TMPFILE}.mapping.reform ${TMPFILE}.mapping.rename mappingbase
    awk '$3 == "mRNA" || $3 == "CDS"' ${TMPFILE}.mapping.rename > ${TMPFILE}.mapping.correction
    ${SCRIPT}/gff_editor -f ${TMPFILE}.mapping.correction -CDSfix > /dev/null
    #${SCRIPT}/ORFframefix -i Gene_region_repair.gff -o ${2}_mappingbase.gff
    mv Gene_region_repair.gff ${2}_mappingbase.gff
    
    #denovo
    echo "RNA-seq de novo based"
    python ${SCRIPT}/gmap_alignment_rate.py $GMAP_RESULT > ${TMPFILE}.exon_id
    python ${SCRIPT}/binning_denovo.py $DENOVO_RESULT ${TMPFILE}.exon_id > ${TMPFILE}.merge_id
    python ${SCRIPT}/scoring_denovo.py ${TMPFILE}.merge_id > ${TMPFILE}.denovo
    ${SCRIPT}/RNA-seq_reform ${TMPFILE}.denovo ${TMPFILE}.denovo.reform
    ${SCRIPT}/Row2_rename ${TMPFILE}.denovo.reform ${TMPFILE}.denovo.rename denovobase
    awk '$3 == "mRNA" || $3 == "CDS"' ${TMPFILE}.denovo.rename > ${TMPFILE}.denovo.correction
    ${SCRIPT}/gff_editor -f ${TMPFILE}.denovo.correction -CDSfix > /dev/null
    #${SCRIPT}/ORFframefix -i Gene_region_repair.gff -o ${2}_denovobase.gff
    mv Gene_region_repair.gff ${2}_denovobase.gff

    #homology
    echo "Homology based"
    python ${SCRIPT}/scoring_homology.py $HOMOLOGY_RESULT > ${TMPFILE}.homology
    ${SCRIPT}/Spaln_reform ${TMPFILE}.homology ${TMPFILE}.homology.reform
    ${SCRIPT}/Row2_rename ${TMPFILE}.homology.reform ${TMPFILE}.homology.rename homology
    awk '$3 == "mRNA" || $3 == "CDS"' ${TMPFILE}.homology.rename > ${TMPFILE}.homology.correction
    ${SCRIPT}/gff_editor -f ${TMPFILE}.homology.correction -CDSfix > /dev/null
    #${SCRIPT}/ORFframefix -i Gene_region_repair.gff -o ${2}_homology.gff
    mv Gene_region_repair.gff ${2}_homology.gff

    #Augustus
    echo "Ab initio based; Augustus"
    ${SCRIPT}/Augustus_reform $AUGUSTUS_RESULT ${TMPFILE}.augustus.reform
    ${SCRIPT}/Row2_rename ${TMPFILE}.augustus.reform ${TMPFILE}.augustus.rename AUGUSTUS
    awk '$3 == "mRNA" || $3 == "CDS"' ${TMPFILE}.augustus.rename > ${TMPFILE}.augustus.correction
    ${SCRIPT}/gff_editor -f ${TMPFILE}.augustus.correction -CDSfix > /dev/null
    python ${SCRIPT}/scoring_augustus.py Gene_region_repair.gff > ${2}_augustus.gff
    rm Gene_region_repair.gff

    #Snap
    echo "Ab initio based; SNAP"
    python ${SCRIPT}/score_separator_v2.py $SNAP_RESULT > ${TMPFILE}.sep
    python ${SCRIPT}/scoring_snap.py ${TMPFILE}.sep > ${TMPFILE}.snap
    ${SCRIPT}/Row2_rename ${TMPFILE}.snap ${TMPFILE}.snap.rename SNAP
    awk '$3 == "mRNA" || $3 == "CDS"' ${TMPFILE}.snap.rename > ${TMPFILE}.snap.correction
    ${SCRIPT}/gff_editor -f ${TMPFILE}.snap.correction -CDSfix > /dev/null
    mv Gene_region_repair.gff ${2}_snap.gff

    #merge
    echo "merging all gff"
    cat ${2}_mappingbase.gff ${2}_denovobase.gff ${2}_homology.gff ${2}_augustus.gff ${2}_snap.gff > ${TMPFILE}.all
    python ${SCRIPT}/weighten.py ${TMPFILE}.all ${TMPFILE}.weight > ${2}_all.gff
fi

trap "rm_tmpfile" EXIT INT PIPE TERM
echo "finished"
