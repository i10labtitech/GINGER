#!/bin/sh

util="/data/fumiya/annotation/mapping_tool_ver2/util"
if test $# -ne 3 -a $# -ne 4; then
    echo '1: target.gff <STRING>'
    echo '2: output file name'
    echo '3: the number of sample size'
    echo '4: seed number (default = 14)'
else
    #gffからmRNAのリストを作成
    awk '$3 == "mRNA"{print $9}' ${1} > ${2}.tag
    ${util}/tag_trimmer ${2}.tag ID > ${2}.id

    #ランダマイズ
    if test $# -eq 3; then
        ${util}/random ${2}.id ${3} > ${2}.rand
    else
        ${util}/random ${2}.id ${3} ${4} > ${2}.rand
    fi

    #gffを復元
    ${util}/gff_trimmer ${1} ${2}.rand ${2}

    #中間ファイルの削除
    rm -f ${2}.tag
    rm -f ${2}.id
    rm -f ${2}.rand
fi
