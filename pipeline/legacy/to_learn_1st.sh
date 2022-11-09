#!/bin/sh

# --- 以下の変数は適宜変更する ---
ROOT=$(dirname `which $0`)
util=${ROOT}/util
GFFREAD="gffread"
CDHIT="cd-hit"

#MINはtranscriptの最小base長
MIN=300
# --- ---

if test $# -ne 3; then
    echo '1: stringtie.gtf <STRING> (output from mapping.sh)'
    echo '2: genome.fa <STRING>'
    echo '3: prefix in output file name <STRING>'
else
    #genome配列をupper caseで統一
    ${util}/seqkit seq -u ${2} > ${3}_tmpgenome.fa

    #multi exon transcriptの取得
    ${util}/exon_num_filter ${1} ${3}_multi.gtf 2 10000

    #gtfからmRNA配列fastaを作成
    ${util}/gtf_genome_to_cdna_fasta.pl ${3}_multi.gtf ${3}_tmpgenome.fa > ${3}_multi.fa

    #gtf -> gff変換
    ${util}/gtf_to_alignment_gff3.pl ${3}_multi.gtf > ${3}_multi.aln.gff3

    #最長ORF検出
    ${util}/ORF_finder ${3}_multi.fa ${3}_multi_orf ${MIN} true
    
    #1 gene 1 transcriptにするためのstats取得
    ${util}/seqkit fx2tab -nl ${3}_multi_orf.cds | awk '{print $1"\t"$2"\t"$5}' > ${3}_multi_orf.list
    
    #1 gene 1 transcriptのリスト作成
    ${util}/longest_transcript ${3}_multi_orf.list > ${3}_multi_orf_longest.list
    
    #1 gene 1 transcriptのgff作成
    ${util}/gff_trimmer ${3}_multi_orf.gff3 ${3}_multi_orf_longest.list ${3}_multi_orf_longest.gff3
    
    #1 gene 1 transcriptのgff3からfastaを作成
    ${GFFREAD} -x ${3}_multi_orf_longest.cds -g ${3}_multi.fa ${3}_multi_orf_longest.gff3

    #重複配列の削除
    ${CDHIT} -c 1.0 -i ${3}_multi_orf_longest.cds -o ${3}_multi_orf_longest_nr.cds
    
    #配列名のリストを作成
    ${util}/seqkit seq -n ${3}_multi_orf_longest_nr.cds | awk '{print $1}' > ${3}_multi_orf_longest_nr.list
    
    #リストからgffを作成
    ${util}/gff_trimmer ${3}_multi_orf_longest.gff3 ${3}_multi_orf_longest_nr.list ${3}_multi_orf_longest_nr.gff3
    
    #genome baseのgff3を作成
    ${util}/cdna_alignment_orf_to_genome_orf.pl ${3}_multi_orf_longest_nr.gff3 ${3}_multi.aln.gff3 ${3}_multi.fa > ${3}_learn_1st.gff3

    #中間ファイルの削除
    rm -f ${3}_multi.aln.gff3
    rm -f ${3}_multi.fa
    rm -f ${3}_multi.fa.fai
    rm -f ${3}_multi.gtf
    rm -f ${3}_multi_orf.cds
    rm -f ${3}_multi_orf.gff3
    rm -f ${3}_multi_orf.pep
    rm -f ${3}_multi_orf.list
    rm -f ${3}_multi_orf_longest.list
    rm -f ${3}_multi_orf_longest.gff3
    rm -f ${3}_multi_orf_longest.cds
    rm -f ${3}_multi_orf_longest_nr.cds
    rm -f ${3}_multi_orf_longest_nr.cds.clstr
    rm -f ${3}_multi_orf_longest_nr.list
    rm -f ${3}_multi_orf_longest_nr.gff3
    rm -f ${3}_tmpgenome.fa
fi
