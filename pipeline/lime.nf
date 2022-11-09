#!/usr/bin/env nextflow

if (file(params.PDIR_PREP_MAPPING).isDirectory()) {
    error "The publish directory already exists: \"${params.PDIR_PREP_MAPPING}\""
} else {
    print "No publish directory : \"${params.PDIR_PREP_MAPPING}\""
    if (file(params.PDIR_PREP).isEmpty()) {
        print "Making dpublish directory : \"${params.PDIR_PREP}\""
        file(params.PDIR_PREP).mkdir()
    }
    if (file(params.PDIR_PREP_MAPPING).isEmpty()) {
        print "Making dpublish directory : \"${params.PDIR_PREP_MAPPING}\""
        file(params.PDIR_PREP_MAPPING).mkdir()
    }
    if (file(params.PDIR_PREP_DENOVO).isEmpty()) {
        print "Making dpublish directory : \"${params.PDIR_PREP_DENOVO}\""
        file(params.PDIR_PREP_DENOVO).mkdir()
    }
    if (file(params.PDIR_PREP_HOMOLOGY).isEmpty()) {
        print "Making dpublish directory : \"${params.PDIR_PREP_HOMOLOGY}\""
        file(params.PDIR_PREP_HOMOLOGY).mkdir()
    }
    if (file(params.PDIR_PREP_ABINITIO).isEmpty()) {
        print "Making dpublish directory : \"${params.PDIR_PREP_ABINITIO}\""
        file(params.PDIR_PREP_ABINITIO).mkdir()
    }
}

if (file(params.INPUT_GENOME).isFile())    {} else { error "No path : \"${params.INPUT_GENOME}\"" }
if (file(params.INPUT_REPOUT).isFile())    {} else { error "No path : \"${params.INPUT_REPOUT}\"" }
if (file(params.INPUT_RNASEQR1).isFile())  {} else { error "No path : \"${params.INPUT_RNASEQR1}\"" }
if (file(params.INPUT_RNASEQR2).isFile())  {} else { error "No path : \"${params.INPUT_RNASEQR2}\"" }
if (file("${params.AUGUSTUS_SPEC_DIR}/${params.AUGUSTUS_SPEC}").isDirectory()) {error "\"${params.AUGUSTUS_SPEC}\" is in ${params.AUGUSTUS_SPEC_DIR}!\n Rename or remove \"${params.AUGUSTUS_SPEC_DIR}/${params.AUGUSTUS_SPEC}\"."} 

hDataChPrefix  = Channel.from(params.HOMOLOGY_DATA['PREFIX'])
hDataChProtein = Channel.from(params.HOMOLOGY_DATA['PROTEIN'])
hDataChSpalndb = Channel.from(params.HOMOLOGY_DATA['SPALNDB'])

//if (file(params.AUGUSTUS_TRAINING_DATA).isFile()) {} else { error "No path : \"${params.AUGUSTUS_TRAINING_DATA}\"" }
size = params.AUGUSTUS_TRAINING_SIZE
dict = 'Abinitio_predict' //Augustus_prediction以下にできる作業ディレクトリの名前

SCRATCH = params.SCRATCH

// --- Mapping of RNA-seq reads ---
process mapping {
    
    cpus params.N_THREAD
    scratch "${SCRATCH}"

    publishDir "${params.PDIR_PREP_MAPPING_MAPPING}", mode: 'copy', overwrite: false
    
    input:
    path genomeFasta from params.INPUT_GENOME
    path rnaSeqRead1 from params.INPUT_RNASEQR1
    path rnaSeqRead2 from params.INPUT_RNASEQR2
    
    output:
    file "${params.OPREFIX}.gtf" into gtf
    
    shell:
    '''
    mkdir !{params.OPREFIX}_tmp_index
    mkdir tmp_mapping_!{params.OPREFIX}
    
    # mapping
    !{params.HISAT2BUILD} -p !{params.N_THREAD} !{genomeFasta} !{params.OPREFIX}_tmp_index/!{params.OPREFIX}_index
    !{params.HISAT2} -p !{params.N_THREAD} -x !{params.OPREFIX}_tmp_index/!{params.OPREFIX}_index -1 !{rnaSeqRead1} -2 !{rnaSeqRead2} -S tmp_mapping_!{params.OPREFIX}/!{params.OPREFIX}.sam --dta --no-discordant --no-mixed
    !{params.SAMTOOLS} view -bS -@ !{params.N_THREAD} tmp_mapping_!{params.OPREFIX}/!{params.OPREFIX}.sam | !{params.SAMTOOLS} sort -@ !{params.N_THREAD} -o tmp_mapping_!{params.OPREFIX}/!{params.OPREFIX}.sorted.bam
    
    # gene prediction
    !{params.STRINGTIE} -p !{params.N_THREAD} -o !{params.OPREFIX}.gtf tmp_mapping_!{params.OPREFIX}/!{params.OPREFIX}.sorted.bam
    '''
}

// --- Merge ---
process to_merge {

    publishDir "${params.PDIR_PREP_MAPPING_MERGE}", mode: 'copy', overwrite: false

    input:
    path genomeFasta from params.INPUT_GENOME
    path gtf
        
    output:
    file "${params.OPREFIX}_merge.gff3" into gff3toMerge
    
    shell:
    '''
    #genome配列をupper caseで統一
    !{params.UTILPATH_MAPPING}/seqkit seq -u !{genomeFasta} > !{params.OPREFIX}_tmpgenome.fa
    
    #strandの付いていないsingle exonを分離
    awk '$7 == "."' !{gtf} > !{params.OPREFIX}_single.gtf
    
    #gtfからmRNA配列fastaを作成
    !{params.UTILPATH_MAPPING}/gtf_genome_to_cdna_fasta.pl !{params.OPREFIX}_single.gtf !{params.OPREFIX}_tmpgenome.fa > !{params.OPREFIX}_single.fa
    
    #single exon最長ORF検出
    !{params.UTILPATH_MAPPING}/ORF_finder !{params.OPREFIX}_single.fa !{params.OPREFIX}_single_orf !{params.MIN0} false
    
    #single exonのstats作成
    !{params.UTILPATH_MAPPING}/seqkit fx2tab -nl !{params.OPREFIX}_single_orf.cds | awk '{print $1"\t"$2"\t"$5}' > !{params.OPREFIX}_single_orf.list
    
    #1 gene 1 transcriptのリスト作成
    !{params.UTILPATH_MAPPING}/longest_transcript !{params.OPREFIX}_single_orf.list > !{params.OPREFIX}_single_orf_longest.list
    
    #1 gene 1 transcriptになるようにgffをトリム
    !{params.UTILPATH_MAPPING}/gff_trimmer !{params.OPREFIX}_single_orf.gff3 !{params.OPREFIX}_single_orf_longest.list !{params.OPREFIX}_single_orf_longest.gff3
    
    #strandの追記
    !{params.UTILPATH_MAPPING}/strand_replace !{params.OPREFIX}_single_orf_longest.gff3 !{gtf} !{params.OPREFIX}_replace.gtf
    
    #orfが取れないsingle exonを削除
    awk '$7 != "."' !{params.OPREFIX}_replace.gtf > !{params.OPREFIX}_new.gtf
    
    #gtfからmRNA配列fastaを作成
    !{params.UTILPATH_MAPPING}/gtf_genome_to_cdna_fasta.pl !{params.OPREFIX}_new.gtf !{params.OPREFIX}_tmpgenome.fa > !{params.OPREFIX}.fa
    
    #gtf -> gff変換
    !{params.UTILPATH_MAPPING}/gtf_to_alignment_gff3.pl !{params.OPREFIX}_new.gtf > !{params.OPREFIX}.gff3
    
    #最長ORF検出
    # !{params.UTILPATH_MAPPING}/ORF_finder !{params.OPREFIX}.fa !{params.OPREFIX}_orf !{params.MIN0} true
    !{params.TD_LONGORFS} -t !{params.OPREFIX}.fa -m !{params.MIN0} -S
    !{params.TD_PREDICT} -t !{params.OPREFIX}.fa 
    
    #genome baseのgff3に変換
    !{params.UTILPATH_MAPPING}/cdna_alignment_orf_to_genome_orf.pl !{params.OPREFIX}.fa.transdecoder.gff3 !{params.OPREFIX}.gff3 !{params.OPREFIX}.fa > !{params.OPREFIX}_merge.gff3
    '''
}

// --- Learn(1st) ---
process to_learn_1st {

    publishDir "${params.PDIR_PREP_MAPPING_TOLEARN1ST}", mode: 'copy', overwrite: false
    
    input:
    path genomeFasta from params.INPUT_GENOME
    path gtf


    output:
    file "${params.OPREFIX}_learn_1st.gff3" into gff3_1st

    shell:
    '''
    #genome配列をupper caseで統一
    !{params.UTILPATH_MAPPING}/seqkit seq -u !{genomeFasta} > !{params.OPREFIX}_tmpgenome.fa
    
    #multi exon transcriptの取得
    !{params.UTILPATH_MAPPING}/exon_num_filter !{gtf} !{params.OPREFIX}_multi.gtf 2 10000
    
    #gtfからmRNA配列fastaを作成
    !{params.UTILPATH_MAPPING}/gtf_genome_to_cdna_fasta.pl !{params.OPREFIX}_multi.gtf !{params.OPREFIX}_tmpgenome.fa > !{params.OPREFIX}_multi.fa
    
    #gtf -> gff変換
    !{params.UTILPATH_MAPPING}/gtf_to_alignment_gff3.pl !{params.OPREFIX}_multi.gtf > !{params.OPREFIX}_multi.aln.gff3
    
    #最長ORF検出
    !{params.UTILPATH_MAPPING}/ORF_finder !{params.OPREFIX}_multi.fa !{params.OPREFIX}_multi_orf !{params.MIN1} true
    
    #1 gene 1 transcriptにするためのstats取得
    !{params.UTILPATH_MAPPING}/seqkit fx2tab -nl !{params.OPREFIX}_multi_orf.cds | awk '{print $1"\t"$2"\t"$5}' > !{params.OPREFIX}_multi_orf.list
    
    #1 gene 1 transcriptのリスト作成
    !{params.UTILPATH_MAPPING}/longest_transcript !{params.OPREFIX}_multi_orf.list > !{params.OPREFIX}_multi_orf_longest.list
    
    #1 gene 1 transcriptのgff作成
    !{params.UTILPATH_MAPPING}/gff_trimmer !{params.OPREFIX}_multi_orf.gff3 !{params.OPREFIX}_multi_orf_longest.list !{params.OPREFIX}_multi_orf_longest.gff3
    
    #1 gene 1 transcriptのgff3からfastaを作成
    !{params.GFFREAD} -x !{params.OPREFIX}_multi_orf_longest.cds -g !{params.OPREFIX}_multi.fa !{params.OPREFIX}_multi_orf_longest.gff3
    
    #重複配列の削除
    !{params.CD_HIT} -c 1.0 -i !{params.OPREFIX}_multi_orf_longest.cds -o !{params.OPREFIX}_multi_orf_longest_nr.cds
    
    #配列名のリストを作成
    !{params.SEQKIT} seq -n !{params.OPREFIX}_multi_orf_longest_nr.cds | awk '{print $1}' > !{params.OPREFIX}_multi_orf_longest_nr.list
    
    #リストからgffを作成
    !{params.UTILPATH_MAPPING}/gff_trimmer !{params.OPREFIX}_multi_orf_longest.gff3 !{params.OPREFIX}_multi_orf_longest_nr.list !{params.OPREFIX}_multi_orf_longest_nr.gff3
    
    #genome baseのgff3を作成
    !{params.UTILPATH_MAPPING}/cdna_alignment_orf_to_genome_orf.pl !{params.OPREFIX}_multi_orf_longest_nr.gff3 !{params.OPREFIX}_multi.aln.gff3 !{params.OPREFIX}_multi.fa > !{params.OPREFIX}_learn_1st.gff3
    '''
}

// --- Learn(2nd) ---
process to_learn_2nd {

    publishDir "${params.PDIR_PREP_MAPPING_TOLEARN2ND}", mode: 'copy', overwrite: false

    input:
    path genomeFasta  from params.INPUT_GENOME
    path genomeRepOut from params.INPUT_REPOUT
    path gff3_1st

    output:
    file "${params.OPREFIX}_learn_2nd.gff3" into gff3_2nd

    shell:
    '''
    #repeat databaseの構築
    awk '$11 != "Simple_repeat" && $11 != "Low_complexity" {print $5"\t"$6"\t"$7"\t"$11}' !{genomeRepOut} | sed -e '1,3d' > repeat_DB.tsv
    
    #genome databaseの構築
    !{params.SEQKIT} fx2tab -inl !{genomeFasta} | awk '{print $1"\t"$2}' | sed -e '/^$/d' > genome_DB.tsv
    
    #annotation databaseの構築
    awk '$3 == "exon" {print $1"\t"$4"\t"$5"\t"$9}' !{gff3_1st} | sed -e '/^$/d' > annotation_DB.tsv
    
    #repeat重複配列の検索
    !{params.UTILPATH_MAPPING}/repeat_checker genome_DB.tsv repeat_DB.tsv annotation_DB.tsv repeat_output.tsv
    
    #repeatを含むexonのリストを作成
    awk '$2 != 0 {print $1}' repeat_output.tsv > repeat_output.list
    
    #ID抽出
    !{params.UTILPATH_MAPPING}/tag_trimmer repeat_output.list Parent | sort | uniq > repeat_output.id
    
    #全transcriptのID抽出
    awk '$3 == "mRNA" {print $9}' !{gff3_1st} > !{params.OPREFIX}_all.tag
    !{params.UTILPATH_MAPPING}/tag_trimmer !{params.OPREFIX}_all.tag ID > !{params.OPREFIX}_all.id
    
    #差集合をとってnon-repeatのmRNA配列を抽出
    !{params.UTILPATH_MAPPING}/set_difference !{params.OPREFIX}_all.id repeat_output.id norep_output.id
    
    #gff作成
    !{params.UTILPATH_MAPPING}/gff_trimmer !{gff3_1st} norep_output.id !{params.OPREFIX}_learn_2nd.gff3
    '''
}

// --- Trinity ---
process trinity {

    cpus params.N_THREAD
    scratch "${SCRATCH}"

    publishDir "${params.PDIR_PREP_DENOVO_TRINITY}", mode: 'copy', overwrite: false
    
    input:
    path rnaSeqRead1 from params.INPUT_RNASEQR1
    path rnaSeqRead2 from params.INPUT_RNASEQR2
    
    output:
//    file "${params.OPREFIX_TRINITY}"
    file "${params.OPREFIX_TRINITY}.Trinity.fasta" into trinityFasta
    file "${params.OPREFIX_TRINITY}.Trinity.fasta.gene_trans_map"

    script:
    preprocessing = ""
    if (params.SRA_FLAG == 1) {
        preprocessing = """
        awk '{{print (NR%4 == 1) ? "@1_" ++i "/1": \$0}}' ${rnaSeqRead1} > rnaseq_1_renamed.fastq
        awk '{{print (NR%4 == 1) ? "@1_" ++i "/2": \$0}}' ${rnaSeqRead2} > rnaseq_2_renamed.fastq
        """
    } else {
        preprocessing = """
        ln -s ${rnaSeqRead1} rnaseq_1_renamed.fastq
        ln -s ${rnaSeqRead2} rnaseq_2_renamed.fastq
        """
    }

    shell:
    '''
    !{preprocessing}

    # trinityでassemble
    mkdir !{params.OPREFIX_TRINITY}
    # ---- 1. Trinityでde novo transcriptome
    time !{params.TRINITY} --seqType fq --left rnaseq_1_renamed.fastq --right rnaseq_2_renamed.fastq --output !{params.OPREFIX_TRINITY} --CPU !{params.N_THREAD} --max_memory !{params.MAX_MEMORY} --full_cleanup
    '''
}

// --- Oases ---
process oases {

    publishDir "${params.PDIR_PREP_DENOVO_OASES}", mode: 'copy', overwrite: false

    input:
    path rnaSeqRead1 from params.INPUT_RNASEQR1
    path rnaSeqRead2 from params.INPUT_RNASEQR2
        
    output:
    file "${params.OPREFIX}/transcripts.fa" into oasesFasta
    
    shell:
    '''
    !{params.VELVETH} !{params.OPREFIX} 31 -fastq -short -separate !{rnaSeqRead1} !{rnaSeqRead2}
    !{params.VELVETG} !{params.OPREFIX} -read_trkg yes 
    !{params.OASES} !{params.OPREFIX}
    '''
}

// --- De novo ---
process denovo {

    cpus params.N_THREAD

    publishDir "${params.PDIR_PREP_DENOVO}", mode: 'copy', overwrite: false
    
    input:
    path genomeFasta from params.INPUT_GENOME
    path trinityFasta
    path oasesFasta

    output:
    file "${params.OPREFIX}_denovo.cd.fasta.gmap.filtered"
    file "${params.OPREFIX}_final.gff3"

    shell:
    '''
    ### sequences ID simplification ###
    touch !{params.OPREFIX}_denovo.cds
    !{params.SEQKIT} seq -i -m 300 !{trinityFasta} >> !{params.OPREFIX}_denovo.cds
    !{params.SEQKIT} seq -i -m 300 !{oasesFasta} >> !{params.OPREFIX}_denovo.cds

    ### de novo sequences compression ###
    !{params.CD_HIT_EST} -i !{params.OPREFIX}_denovo.cds -c 1.0 -o !{params.OPREFIX}_denovo.cd.fasta -T !{params.N_THREAD} -M 100000

    ### GMAP alignment ###
    !{params.GMAP_BUILD} -D . -d !{params.OPREFIX}_genome_gmap_DB !{genomeFasta} > !{params.OPREFIX}_gmap_build.log 2>&1
    !{params.GMAP} -S -t !{params.N_THREAD} -n 1 -D . -d !{params.OPREFIX}_genome_gmap_DB !{params.OPREFIX}_denovo.cd.fasta > !{params.OPREFIX}_denovo.cd.fasta.gmap 2> !{params.OPREFIX}_gmap.stderr
    
    ### filtering alignment result ###
    !{params.DENOVO_PYTHON} !{params.UTILPATH_DENOVO}/filtering.py !{params.OPREFIX}_denovo.cd.fasta.gmap 95 > !{params.OPREFIX}_denovo.cd.fasta.gmap.filtered

    ### format conversion for ORF prediction
    !{params.UTILPATH_DENOVO}/gmap_native_to_format_converter.pl !{params.OPREFIX}_denovo.cd.fasta.gmap.filtered GTF > !{params.OPREFIX}_denovo.cd.fasta.gtf
    !{params.UTILPATH_DENOVO}/cufflinks_gtf_to_alignment_gff3.pl !{params.OPREFIX}_denovo.cd.fasta.gtf > !{params.OPREFIX}_denovo.cd.fasta.gff3
    !{params.UTILPATH_DENOVO}/cufflinks_gtf_genome_to_cdna_fasta.pl !{params.OPREFIX}_denovo.cd.fasta.gtf !{genomeFasta} > !{params.OPREFIX}_denovo_genome.fa

    ### ORF prediction ###
    !{params.TD_LONGORFS} -t !{params.OPREFIX}_denovo_genome.fa -m 90
    !{params.TD_PREDICT} --no_refine_starts --single_best_only -t !{params.OPREFIX}_denovo_genome.fa
    
    ### format conversion ###
    !{params.UTILPATH_DENOVO}/cdna_alignment_orf_to_genome_orf_nolimit.pl !{params.OPREFIX}_denovo_genome.fa.transdecoder.gff3 !{params.OPREFIX}_denovo.cd.fasta.gff3 !{params.OPREFIX}_denovo_genome.fa > !{params.OPREFIX}_final.gff3
    '''
}

// --- Homology ---
process homology {

    cpus params.N_THREAD
    scratch "${SCRATCH}"

    input:
    path genomeFasta from params.INPUT_GENOME
//    val hd from hDataCh
    val prefix  from hDataChPrefix
    path protein from hDataChProtein
    val spalndb from hDataChSpalndb

    publishDir "${params.PDIR_PREP_HOMOLOGY_HOMOLOGY}/${prefix}", mode: 'copy', overwrite: false

    script:
    println "${prefix},${protein},${spalndb}"
    
    output:
    file "${prefix}_spalnresult.gff" into spalnResults

    shell:
    '''
#    echo "start homology"
#    echo "!{genomeFasta} !{protein} !{spalndb} !{params.N_THREAD} !{params.OPREFIX}"

    #spaln動かすためスクリプトの用意
#    cp !{params.UTILPATH_HOMOLOGY}/Makefile .
#    cp !{params.UTILPATH_HOMOLOGY}/catchr.pl .
    cp !{params.UTILPATH_HOMOLOGY}/makblk.pl . # used in makeidx.pl
#    cp !{params.UTILPATH_HOMOLOGY}/makeidx.pl .
#    cp !{params.UTILPATH_HOMOLOGY}/function.hpp .

#    echo "start spaln"

    #spalnを動かす前処理
    !{params.UTILPATH_HOMOLOGY}/fastarepair !{genomeFasta} refer.mfa
    !{params.UTILPATH_HOMOLOGY}/fastarepair2 !{protein} relate.faa

    #spalnを動かす
    export PATH=./:\$PATH
    !{params.MAKEIDX} -ip refer.mfa

    #exon identity計測用splan result
    !{params.SPALN} -Q7 -O4 -ospalnresult_o4 -M -yS# -T!{spalndb} -yB# -yZ -t!{params.N_THREAD} -drefer relate.faa 
    mv spalnresult_o4 !{params.OPREFIX}_spalnresult_alignment.tsv
    
    #gff作成
    python !{params.UTILPATH_HOMOLOGY}/o4_to_gff.py !{params.OPREFIX}_spalnresult_alignment.tsv !{prefix} > !{prefix}_spalnresult.gff

#    echo "finish_homology"
    '''
}

// --- Merge results maed by homology ---
allHomology = spalnResults
    .collectFile(name: 'all_homology.gff')

// --- Homology filter---
process homology_filter {

    publishDir "${params.PDIR_PREP_HOMOLOGY_HOMOLOGYFILTER}", mode: 'copy', overwrite: false

    input:
    path genomeFasta from params.INPUT_GENOME
    path allHomology

    output:
    file "all_homology_filter.gff"
    
    shell:
    '''
    !{params.SEQKIT} fx2tab -nl !{genomeFasta} | awk '{print $1"\t"$2}' > tmp_genome_stats
    python !{params.UTILPATH_HOMOLOGY}/outer_trim.py !{allHomology} tmp_genome_stats > !{allHomology}.trim
    sed -e 's/cds/CDS/g' !{allHomology}.trim > !{allHomology}.rename
    !{params.UTILPATH_HOMOLOGY}/gff_2_proteinfasta -i !{allHomology}.rename -g !{genomeFasta} -o prefilter01.fa
    !{params.UTILPATH_HOMOLOGY}/flameshiftfilter prefilter01.fa flameshiftlist.txt
    sed -e 's/mRNA/gene/g' flameshiftlist.txt > flameshiftlist.rename
    python !{params.UTILPATH_HOMOLOGY}/gff_trim.py !{allHomology}.rename flameshiftlist.rename > all_homology_filter.gff.rename
    sed -e 's/CDS/cds/g' all_homology_filter.gff.rename > all_homology_filter.gff
    '''
}

// --- Augustus ---
process augustus {

    cpus params.N_THREAD
    scratch "${SCRATCH}"
    stageOutMode 'copy'
    
    publishDir "${params.PDIR_PREP_ABINITIO_AUGUSTUS}", mode: 'copy', overwrite: false
//    publishDir "${params.PDIR_PREP_ABINITIO}", mode: 'copy', overwrite: false
    
    input:
    path trainingData from gff3_2nd
    path genomeFASTA  from params.INPUT_GENOME
    path genomeRepOUT from params.INPUT_REPOUT
    
    output:
//    file "${CWD0}/PASS_information_for_snap_t${size}.txt" into augusutusOutput // To be ...
//    file "PASS_information_for_snap_t${size}.txt" into augusutusOutput // To be ...
//    path "Augustus_prediction/${dict}" into augustusDictDir
    path "Augustus_abinitio.gtf"
    path "second.gb.train.train" into augustusTrainTrain
    
    shell:
    '''
    export AUGUSTUS_CONFIG_PATH="!{params.AUGUSTUS_CONFIG_DIR}"
    export PATH="$PATH:!{params.AUGUSTUS_SCRIPT_DIR}"

    ###----事前環境-------      //
    script=!{params.UTILPATH_ABINITIO} #用いるスクリプトのパス #/data/yuasa/script/abinitio/script
    ###-------------------
    #マスクしたゲノム配列作成
    ${script}/simple_low_norepeatmask !{genomeFASTA} !{genomeRepOUT} !{genomeFASTA}.masked_without_SR_LC.fa
    genome=!{genomeFASTA}.masked_without_SR_LC.fa
    ########################
    
    CWD0=`pwd`

    #Augustusの出力ディレクトリの作成・移動
    mkdir -p Augustus_prediction
    cd Augustus_prediction
    
    #作業ディレクトリの作成・移動
    mkdir !{dict}
    cd !{dict}
    
    #SNAPに渡すファイルのディレクトリーを示したファイルの作成
    CWD1=`pwd`
    
    #gff3をgtfに変換
    python ${script}/gff2gtf_better.py ${CWD0}/!{trainingData} > ./input.gtf
    #学習に用いるCDSの情報のみにする
    awk '$3=="CDS"' ./input.gtf > ./input2.gtf
    #gtfをgbに変換(遺伝子周辺領域をmax1000bpとして切り出す)
    perl !{params.AUGUSTUS_SCRIPT_DIR}/gff2gbSmallDNA.pl input2.gtf ${CWD0}/${genome} 1000 first.gb
    
    #trainingしたパラメータを入れるディレクトリ作成
    perl !{params.AUGUSTUS_SCRIPT_DIR}/new_species.pl --species=!{params.AUGUSTUS_SPEC}
    
    #Augustusの入力にふさわしくない遺伝子情報のフィルタリング
    #スタートコドンがなかったり、終始コドンが無いようなエラーを取り除く(genericはAugustusのマーカーセット)参考：Scipio
    #--stopCodonExcludedFromCDS=trueはCDSにストップコドンを含めているか否かの設定がetrainingでできる。train.errファイルを確認しエラーが多く検出されていれば確認した方が良い。
    !{params.ETRAINING} --species=generic first.gb 2> train.err
    fgrep "gene" train.err | cut -f 2 -d " " > bad.etraining-test.lst
    perl !{params.AUGUSTUS_SCRIPT_DIR}/filterGenesOut_mRNAname.pl bad.etraining-test.lst first.gb > second.gb
    
    #テスト用のファイル切り出し(テストファイルサイズ100)
    perl !{params.AUGUSTUS_SCRIPT_DIR}/randomSplit.pl second.gb 100
    
    line_num=`grep "LOCUS" second.gb.train | wc -l`
    num=`echo $((${line_num} - !{size}))`
    
    perl !{params.AUGUSTUS_SCRIPT_DIR}/randomSplit.pl second.gb.train ${num}
	
    #training1
    !{params.ETRAINING} --species=!{params.AUGUSTUS_SPEC} second.gb.train.train
	
    #テスト
    !{params.AUGUSTUS} --species=!{params.AUGUSTUS_SPEC} second.gb.test | tee firsttest.out
    grep -A 22 Evaluation firsttest.out > test_result1.txt &
	    
    #最適化
    perl !{params.AUGUSTUS_SCRIPT_DIR}/optimize_augustus.pl --species=!{params.AUGUSTUS_SPEC} second.gb.train.train --cpus=!{params.N_THREAD} --kfold=!{params.N_THREAD}
	
    #training2
    !{params.ETRAINING} --species=!{params.AUGUSTUS_SPEC} second.gb.train.train
    
    #テスト2
    !{params.AUGUSTUS} --species=!{params.AUGUSTUS_SPEC} second.gb.test | tee second.out
    grep -A 22 Evaluation second.out > test_result2.txt &
	    
    #fastaファイルの分割
    python ${script}/pre_assembly_first_step.py ${CWD0}/${genome} !{params.N_THREAD}
    sort -nrk 1 -t "@" list.tMp > list_sort.tMp
    wait
    python ${script}/distribute_second_step.py list_sort.tMp !{params.N_THREAD}
    
    for y in $(seq 1 !{params.N_THREAD})
    do
	    python ${script}/tab_separate_third_step.py sPlitfile${y}.tMp > sPlitfile${y}.fasta &
    done
    wait

    #分割したfastaファイルで予測
    for k in $(seq 1 !{params.N_THREAD})
    do
        !{params.AUGUSTUS} --noInFrameStop=true --species=!{params.AUGUSTUS_SPEC} sPlitfile${k}.fasta  > abinitio${k}.gff &
    done
    wait
    
    #分割して予測したファイルを統合
    for l in $(seq 1 !{params.N_THREAD})
    do
        grep -v "#" abinitio${l}.gff > abinitio${l}_shaped.gff &
    done
    wait

    #予測結果のファイルを統合(IDの付け方はアイソフォームを考慮していない)
    cat abinitio*_shaped.gff > abinitio_shaped_cat.gff
    python ${script}/re_name_for_AugustusGFF.py abinitio_shaped_cat.gff > Augustus_abinitio_tmp.gtf

    #in-frame stopcodon geneを取り除く
    ${script}/inframe_stopcodon_exclude -i Augustus_abinitio_tmp.gtf -g ${CWD0}/!{genomeFASTA} -d1 gene -d2 CDS -o Augustus_abinitio.gtf

    # Finish
    cp Augustus_abinitio.gtf ../../Augustus_abinitio.gtf
    cp second.gb.train.train ../../second.gb.train.train
    '''
}

// --- SNAP ---
process SNAP {

    scratch "${SCRATCH}"

    publishDir "${params.PDIR_PREP_ABINITIO_SNAP}", mode: 'copy', overwrite: false
//    publishDir "${params.PDIR_PREP_ABINITIO}", mode: 'copy', overwrite: false
    
    input:
    //    path trainingData       from params.SNAP_TRAINING_DATA // mapping.nfのアウトプット　mapping/to_learn_2nd/lemon_learn_2nd.gff3
    path trainingData       from gff3_2nd
    path genomeFASTA        from params.INPUT_GENOME
    path genomeMaskedFASTA  from params.INPUT_MASKEDGENOME
    path augustusTrainTrain
    
    output:
    file "SNAP/snap_work/final_result/${params.SNAP_SPECIES}.original.gff"
    file "SNAP/snap_work/final_result/${params.SNAP_SPECIES}.hmm"
    file "SNAP/snap_work/final_result/${params.SNAP_SPECIES}.faa"
    file "SNAP/snap_work/final_result/${params.SNAP_SPECIES}.final.gff"
    file "${params.SNAP_SPECIES}.final.gff"

    shell:
    '''
    export PATH="!{params.SNAP_DIR}:$PATH" #SNAPが使えるようにパスを通す。

    ####事前環境設定####
    script=!{params.UTILPATH_ABINITIO}

    ####

    CWD0=`pwd`

    mkdir -p SNAP
    cd SNAP

    CWD1=`pwd`

    python ${script}/down_size.py ${CWD0}/!{augustusTrainTrain} ${CWD0}/!{trainingData} 
    num_locus=`grep "LOCUS" ${CWD0}/!{augustusTrainTrain} | wc -l`
    gff=down_size${num_locus}.gff

    #出力先のディレクトリ作成
    mkdir snap_work
    cd snap_work
    
    #入力データ作成のためのディレクトリ作成
    mkdir -p make_data_set
    cd make_data_set
    
    #ゲノムに含まれる配列名でリスト作成
    awk '{print $1}' ${CWD1}/${gff} > genome.pre.txt
    python ${script}/same_line_remover.py genome.pre.txt > genome.txt

    #改行を取り除いたfastaファイルの作成
    python ${script}/rmk_fasta_maker_2018NOV11.py ${CWD0}/!{genomeMaskedFASTA}

    echo "##FASTA" > tmp.txt

    #配列名ごとにMaker独自のgffフォーマットを作成
    #gffを配列ごとに作成、sortして、末尾に配列をくっ付ける。
    while read line
    do
    	#gffから対象の配列だけ抜き出す。
	    python ${script}/first_gff_picker.py ${CWD1}/${gff} ${line} > ${line}.st.tmp.gff
		#mRNAの位置でソートするために、同じ遺伝子データを一行に加工
		python ${script}/second_gff_liner.py ${line}.st.tmp.gff > ${line}.one_line.tmp.gff
		#ソート
		sort -nk 4 ${line}.one_line.tmp.gff > ${line}.sort.tmp.gff
		#一行にしていたデータを元に戻す
		python ${script}/third_gff_reset.py ${line}.sort.tmp.gff > ${line}.en.tmp.gff
		#genome配列から対象の配列のみを抽出
		python ${script}/fasta_pickup_2018NOV11.py ${CWD0}/!{genomeMaskedFASTA} ${line} > ${line}.fasta
#		#Maker独自のgffフォーマットを作り出す
#		cat ${line}.en.tmp.gff tmp.txt ${line}.fasta > ${line}.gff
#		#中間ファイルの掃除
#		rm *.tmp.gff &
#		#Maker形式のgffからzff形式に変換
#		!{params.MAKER2ZFF} ${line}.gff
        perl ${script}/gff2ann.pl ${line}.en.tmp.gff > ${line}.ann
        mv ${line}.fasta ${line}.dna
#		#中間gffの削除
#		rm ${line}.gff &
#		#名前が混同しないようにファイル名genomeを配列名に変更
#		mv genome.ann ${line}.ann
#		mv genome.dna ${line}.dna
#		wait
	done < ./genome.txt

    cd ../
    
    #SNAPの遺伝子予測の結果（HMMとgffなど）を出力するディレクトリの作成
    mkdir snap_prediction
    
    #配列ごとにばらけていたファイルを統合し、最終ファイルを作り出す
    cat make_data_set/*.ann > snap_prediction/genome.ann
    cat make_data_set/*.dna > snap_prediction/genome.dna
    wait

    cd snap_prediction/
        
    #学習セットの遺伝子の特徴を記述[スキップ可]
    !{params.FATHOM} genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
    #学習セットとして不適当な可能性のある遺伝子のエラーログをつける。[スキップ可]
    !{params.FATHOM} genome.ann genome.dna -validate > validate.log 2>&1
    #一領域一遺伝子となるようにカットする。
    !{params.FATHOM} genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
    #方向をplusに統一する。
    !{params.FATHOM} uni.ann uni.dna -export 1000 -plus> uni-plus.log 2>&1
    #遺伝子予測用のファイル構築
    !{params.FORGE} export.ann export.dna > forge.log 2>&1
    #学習
    perl !{params.SNAP_DIR}/hmm-assembler.pl !{params.SNAP_SPECIES} . > !{params.SNAP_SPECIES}.hmm
    #遺伝子予測(gff作成;出力ファイルはCDS )
    !{params.SNAP} !{params.SNAP_SPECIES}.hmm ${CWD0}/!{genomeMaskedFASTA} -gff -aa !{params.SNAP_SPECIES}.faa > !{params.SNAP_SPECIES}.gff
    
    #SNAPのgffをEVMのgffに変換
    python ${script}/SNAP_to_GFF3.py !{params.SNAP_SPECIES}.gff > !{params.SNAP_SPECIES}.lemonfmt.gff
    
    #gffにフレーム情報を追加&inframe遺伝子を取り除く(2020/10/6変更箇所：Inframe stop codonのチェックを抜く。そのような遺伝子をSNAPが予測しないため。)
    python ${script}/frame_checker_and_inframe_eliminator.snap.lemonfmt.py ${CWD0}/!{genomeMaskedFASTA} !{params.SNAP_SPECIES}.faa !{params.SNAP_SPECIES}.lemonfmt.gff > !{params.SNAP_SPECIES}.lemonfmt.added_frame_infomation.gff
    sed '/^$/d' !{params.SNAP_SPECIES}.lemonfmt.added_frame_infomation.gff > !{params.SNAP_SPECIES}.pre_final.gff
    ${script}/makefasta -f -i !{params.SNAP_SPECIES}.pre_final.gff -g ${CWD0}/!{genomeFASTA} -o !{params.SNAP_SPECIES}.pre_final.fna
    python ${script}/terminal_exon_to_cds_trouble_fix_atSNAP.py !{params.SNAP_SPECIES}.pre_final.fna !{params.SNAP_SPECIES}.faa !{params.SNAP_SPECIES}.pre_final.gff > !{params.SNAP_SPECIES}.final.gff
    ##########################
    cd ../
    #最終結果をまとめる
    mkdir -p final_result
    cd final_result
    cp ../snap_prediction/!{params.SNAP_SPECIES}.gff ./!{params.SNAP_SPECIES}.original.gff
    cp ../snap_prediction/!{params.SNAP_SPECIES}.hmm ./!{params.SNAP_SPECIES}.hmm
    cp ../snap_prediction/!{params.SNAP_SPECIES}.faa ./!{params.SNAP_SPECIES}.faa
    cp ../snap_prediction/!{params.SNAP_SPECIES}.final.gff ./!{params.SNAP_SPECIES}.final.gff
    cp  ./!{params.SNAP_SPECIES}.final.gff ../../../!{params.SNAP_SPECIES}.final.gff
    '''
}
