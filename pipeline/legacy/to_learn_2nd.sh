#!/bin/sh

# --- 以下の変数は適宜変更する ---
ROOT=$(dirname `which $0`)
util=${ROOT}/util
#util="/data2/fumiya/annotation/mapping_tool_ver2/util"
# --- ----

if test $# -ne 4; then
    echo '1: learn_1st.gff3 <STRING>'
    echo '2: genome.fasta <STRING>'
    echo '3: genome.fa.out <STRING>'
    echo '4: prefix in output file name <STRING>'
else
    #repeat databaseの構築
    awk '$11 != "Simple_repeat" && $11 != "Low_complexity" {print $5"\t"$6"\t"$7"\t"$11}' ${3} | sed -e '1,3d' > repeat_DB.tsv

    #genome databaseの構築
    ${util}/seqkit fx2tab -inl ${2} | awk '{print $1"\t"$2}' | sed -e '/^$/d' > genome_DB.tsv

    #annotation databaseの構築
    awk '$3 == "exon" {print $1"\t"$4"\t"$5"\t"$9}' ${1} | sed -e '/^$/d' > annotation_DB.tsv

    #repeat重複配列の検索
    ${util}/repeat_checker genome_DB.tsv repeat_DB.tsv annotation_DB.tsv repeat_output.tsv

    #repeatを含むexonのリストを作成
    awk '$2 != 0 {print $1}' repeat_output.tsv > repeat_output.list

    #ID抽出
    ${util}/tag_trimmer repeat_output.list Parent | sort | uniq > repeat_output.id

    #全transcriptのID抽出
    awk '$3 == "mRNA" {print $9}' ${1} > ${4}_all.tag
    ${util}/tag_trimmer ${4}_all.tag ID > ${4}_all.id

    #差集合をとってnon-repeatのmRNA配列を抽出
    ${util}/set_difference ${4}_all.id repeat_output.id norep_output.id

    #gff作成
    ${util}/gff_trimmer ${1} norep_output.id ${4}_learn_2nd.gff3

    #中間ファイルの削除
    rm -f repeat_DB.tsv
    rm -f genome_DB.tsv
    rm -f annotation_DB.tsv
    rm -f repeat_output.tsv
    rm -f repeat_output.list
    rm -f repeat_output.id
    rm -f norep_output.id
    rm -f ${4}_all.tag
    rm -f ${4}_all.id
fi
