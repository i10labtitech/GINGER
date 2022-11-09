ROOT=$(dirname $( readlink -f $0 ))
Spaln=${ROOT}/Spaln_reform
Augus=${ROOT}/Augustus_reform
Trans=${ROOT}/RNA-seq_reform
Edit=${ROOT}/gff_editor
Row2=${ROOT}/Row2_rename
ORF=${ROOT}/ORFframefix
Frame=${ROOT}/frame_reform

if [ $# -lt 2 ]; then
    echo "1 ...file  2 ...mode  3 ...rowname" 
    echo "[Mode Option]"
    echo "1 homology 2 mapping/denovo 3 Augustus 4 SNAP 5 highppv 6 homo highppv"
    exit 1
fi

file="${3}.gff"
ln -s $1 ${file}

if test $2 -eq 1; then
    $Spaln $file ${file/gff/fix1.gff}
    $Row2 ${file/gff/fix1.gff} ${file/gff/fix2.gff} $3
    rm ${file/gff/fix1.gff}
elif test $2 -eq 2; then
    $Trans $file ${file/gff/fix1.gff}
    $Row2 ${file/gff/fix1.gff} ${file/gff/fix2.gff} $3
    rm ${file/gff/fix1.gff}
elif test $2 -eq 3; then
    $Augus $file ${file/gff/fix1.gff}
    $Row2 ${file/gff/fix1.gff} ${file/gff/fix2.gff} $3
    rm ${file/gff/fix1.gff}
elif test $2 -eq 4; then
    cp $file ${file/gff/fix1.gff}
    $Row2 ${file/gff/fix1.gff} ${file/gff/fix2.gff} $3
    rm ${file/gff/fix1.gff}
elif test $2 -eq 5; then
    cp $file ${file/gff/fix1.gff}    
    $Row2 ${file/gff/fix1.gff} ${file/gff/fix2.gff} $3
    rm ${file/gff/fix1.gff}
elif test $2 -eq 6; then
    $Spaln $file ${file/gff/fix1.gff}
    $Row2 ${file/gff/fix1.gff} ${file/gff/fix2.gff} $3
    rm ${file/gff/fix1.gff}
fi

awk '{if($3 == "mRNA" || $3 == "CDS"){print $0}}' ${file/gff/fix2.gff} > ${file/gff/fix3.gff}
$Edit -f ${file/gff/fix3.gff} -CDSfix
mv Gene_region_repair.gff ${file/gff/fix4.gff}

if test $2 -eq 3 -o $2 -eq 4 -o $2 -eq 5; then
    $Frame -i ${file/gff/fix4.gff} -o ${file/gff/fix5.gff}
    rm ${file/gff/fix2.gff} ${file/gff/fix3.gff} ${file/gff/fix4.gff}
else
    $ORF -i ${file/gff/fix4.gff} -o ${file/gff/fix5.gff}
    rm ${file/gff/fix2.gff} ${file/gff/fix3.gff} ${file/gff/fix4.gff}
fi

$Edit -f ${file/gff/fix5.gff} -num 1 -v
mv filtered.gff `basename ${file/gff/fix5.gff}.rmsin`

$Edit -f ${file/gff/fix5.gff} -num 1
mv filtered.gff `basename ${file/gff/fix5.gff}.sin`
