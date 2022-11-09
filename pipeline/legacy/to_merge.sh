#!/bin/sh

# --- 以下の変数は適宜変更する ---
util="/data1/LEMON/newLemon_v2.1/Prep/mapping/util"
TD_LONGORFS="TransDecoder.LongOrfs"
TD_PREDICT="TransDecoder.Predict"
#MINはtranscriptの最小base長
MIN=90
# --- ---

if test $# -ne 3; then
    echo '1: stringtie.gtf <STRING> (output from mapping.sh)'
    echo '2: genome.fasta <STRING>'
    echo '3: prefix in output file name <STRING>'
else
    #genome配列をupper caseで統一
    ${util}/seqkit seq -u ${2} > ${3}_tmpgenome.fa

    #strandの付いていないsingle exonを分離
    awk '$7 == "."' ${1} > ${3}_single.gtf

    #gtfからmRNA配列fastaを作成
    ${util}/gtf_genome_to_cdna_fasta.pl ${3}_single.gtf ${3}_tmpgenome.fa > ${3}_single.fa
    
    #single exon最長ORF検出
    ${util}/ORF_finder ${3}_single.fa ${3}_single_orf ${MIN} false
    
    #single exonのstats作成
    ${util}/seqkit fx2tab -nl ${3}_single_orf.cds | awk '{print $1"\t"$2"\t"$5}' > ${3}_single_orf.list

    #1 gene 1 transcriptのリスト作成
    ${util}/longest_transcript ${3}_single_orf.list > ${3}_single_orf_longest.list
    
    #1 gene 1 transcriptになるようにgffをトリム
    ${util}/gff_trimmer ${3}_single_orf.gff3 ${3}_single_orf_longest.list ${3}_single_orf_longest.gff3

    #strandの追記
    ${util}/strand_replace ${3}_single_orf_longest.gff3 ${1} ${3}_replace.gtf

    #orfが取れないsingle exonを削除
    awk '$7 != "."' ${3}_replace.gtf > ${3}_new.gtf
    
    #gtfからmRNA配列fastaを作成
    ${util}/gtf_genome_to_cdna_fasta.pl ${3}_new.gtf ${3}_tmpgenome.fa > ${3}.fa

    #gtf -> gff変換
    ${util}/gtf_to_alignment_gff3.pl ${3}_new.gtf > ${3}.gff3

    #最長ORF検出
    #${util}/ORF_finder ${3}.fa ${3}_orf ${MIN} true
    $TD_LONGORFS -t ${3}.fa -m $MIN -S
    $TD_PREDICT -t ${3}.fa 

    #genome baseのgff3に変換
    ${util}/cdna_alignment_orf_to_genome_orf.pl ${3}.fa.transdecoder.gff3 ${3}.gff3 ${3}.fa > ${3}_merge.gff3
    
    #中間ファイルの削除
    rm -f ${3}.fa
    rm -f ${3}.gff3
    rm -f ${3}_new.gtf
    rm -f ${3}_orf.cds
    rm -f ${3}_orf.gff3
    rm -f ${3}_replace.gtf
    rm -f ${3}_single.fa
    rm -f ${3}_single.gtf
    rm -f ${3}_single_orf.cds
    rm -f ${3}_single_orf.gff3
    rm -f ${3}_single_orf.list
    rm -f ${3}_single_orf_longest.gff3
    rm -f ${3}_single_orf_longest.list
    rm -f ${3}_tmpgenome.fa
fi
