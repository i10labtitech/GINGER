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
trap rm_tmpfile EXIT INT PIPE TERM

GINGER_PATH=`dirname $(readlink -f $0)`
SCRIPT=${GINGER_PATH}/ginger-util
CDIR=`pwd`

PREFIX="ginger"

#import config
echo "perl ${SCRIPT}/genConfigIniForMerge.pl ${NEXTFLOWCONFIG} > config.ini"
perl ${SCRIPT}/genConfigIniForMerge.pl ${NEXTFLOWCONFIG} > config.ini
source ${CDIR}/config.ini

#mapping
echo "RNA-seq mapping based"
python ${SCRIPT}/cov_from_gtf.py $MAPPING_GTF > ${TMPFILE}.cov_stats
python ${SCRIPT}/binning_mapping.py $MAPPING_TO_MERGE ${TMPFILE}.cov_stats > ${TMPFILE}.merge_cov
python ${SCRIPT}/scoring_mapping.py ${TMPFILE}.merge_cov ${TMPFILE}.cov_stats > ${TMPFILE}.mapping
${SCRIPT}/RNA-seq_reform ${TMPFILE}.mapping ${TMPFILE}.mapping.reform
${SCRIPT}/Row2_rename ${TMPFILE}.mapping.reform ${TMPFILE}.mapping.rename mappingbase
awk '$3 == "mRNA" || $3 == "CDS"' ${TMPFILE}.mapping.rename > ${TMPFILE}.mapping.correction
${SCRIPT}/phase0_gff_editor -f ${TMPFILE}.mapping.correction -CDSfix > /dev/null
mv Gene_region_repair.gff ${PREFIX}_mappingbase.gff

#denovo
echo "RNA-seq de novo based"
python ${SCRIPT}/gmap_alignment_rate.py $GMAP_RESULT > ${TMPFILE}.exon_id
python ${SCRIPT}/binning_denovo.py $DENOVO_RESULT ${TMPFILE}.exon_id > ${TMPFILE}.merge_id
python ${SCRIPT}/scoring_denovo.py ${TMPFILE}.merge_id > ${TMPFILE}.denovo
${SCRIPT}/RNA-seq_reform ${TMPFILE}.denovo ${TMPFILE}.denovo.reform
${SCRIPT}/Row2_rename ${TMPFILE}.denovo.reform ${TMPFILE}.denovo.rename denovobase
awk '$3 == "mRNA" || $3 == "CDS"' ${TMPFILE}.denovo.rename > ${TMPFILE}.denovo.correction
${SCRIPT}/phase0_gff_editor -f ${TMPFILE}.denovo.correction -CDSfix > /dev/null
mv Gene_region_repair.gff ${PREFIX}_denovobase.gff

#homology
echo "Homology based"
python ${SCRIPT}/scoring_homology.py $HOMOLOGY_RESULT > ${TMPFILE}.homology
${SCRIPT}/Spaln_reform ${TMPFILE}.homology ${TMPFILE}.homology.reform
${SCRIPT}/Row2_rename ${TMPFILE}.homology.reform ${TMPFILE}.homology.rename homology
awk '$3 == "mRNA" || $3 == "CDS"' ${TMPFILE}.homology.rename > ${TMPFILE}.homology.correction
${SCRIPT}/phase0_gff_editor -f ${TMPFILE}.homology.correction -CDSfix > /dev/null
mv Gene_region_repair.gff ${PREFIX}_homology.gff

#Augustus
echo "Ab initio based; Augustus"
${SCRIPT}/Augustus_reform $AUGUSTUS_RESULT ${TMPFILE}.augustus.reform
${SCRIPT}/Row2_rename ${TMPFILE}.augustus.reform ${TMPFILE}.augustus.rename AUGUSTUS
awk '$3 == "mRNA" || $3 == "CDS"' ${TMPFILE}.augustus.rename > ${TMPFILE}.augustus.correction
${SCRIPT}/phase0_gff_editor -f ${TMPFILE}.augustus.correction -CDSfix > /dev/null
python ${SCRIPT}/scoring_augustus.py Gene_region_repair.gff > ${PREFIX}_augustus.gff
rm Gene_region_repair.gff

#Snap
echo "Ab initio based; SNAP"
python ${SCRIPT}/score_separator_v2.py $SNAP_RESULT > ${TMPFILE}.sep
python ${SCRIPT}/scoring_snap.py ${TMPFILE}.sep > ${TMPFILE}.snap
${SCRIPT}/Row2_rename ${TMPFILE}.snap ${TMPFILE}.snap.rename SNAP
awk '$3 == "mRNA" || $3 == "CDS"' ${TMPFILE}.snap.rename > ${TMPFILE}.snap.correction
${SCRIPT}/phase0_gff_editor -f ${TMPFILE}.snap.correction -CDSfix > /dev/null
mv Gene_region_repair.gff ${PREFIX}_snap.gff

#merge
echo "merging all gff"
cat ${PREFIX}_mappingbase.gff ${PREFIX}_denovobase.gff ${PREFIX}_homology.gff ${PREFIX}_augustus.gff ${PREFIX}_snap.gff > ${TMPFILE}.all
perl ${SCRIPT}/addAdditionalGff.pl ${NEXTFLOWCONFIG} >> ${TMPFILE}.all
perl ${SCRIPT}/weighten.pl ${TMPFILE}.all ${NEXTFLOWCONFIG} > ${PREFIX}_all.gff

echo "finished"
