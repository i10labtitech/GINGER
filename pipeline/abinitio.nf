#!/usr/bin/env nextflow

if (file(params.PDIR_PREP_ABINITIO).isDirectory()) {
    error "The publish directory already exists: \"${params.PDIR_PREP_ABINITIO}\""
} else {
    print "No publish directory : \"${params.PDIR_PREP_ABINITIO}\""
    if (file(params.PDIR_PREP).isEmpty()) {
        print "Making dpublish directory : \"${params.PDIR_PREP}\""
        file(params.PDIR_PREP).mkdir()
    }
    if (file(params.PDIR_PREP_ABINITIO).isEmpty()) {
        print "Making dpublish directory : \"${params.PDIR_PREP_ABINITIO}\""
        file(params.PDIR_PREP_ABINITIO).mkdir()
    }
}

if (file(params.AUGUSTUS_TRAINING_DATA).isFile()) {} else { error "No path : \"${params.AUGUSTUS_TRAINING_DATA}\"" }
if (file(params.INPUT_GENOME).isFile())           {} else { error "No path : \"${params.INPUT_GENOME}\"" }
if (file(params.INPUT_REPOUT).isFile())           {} else { error "No path : \"${params.INPUT_REPOUT}\"" }
//if (file(params.INPUT_TOLEARN1ST).isFile())       {} else { error "No path : \"${params.TOLEARN1ST}\"" }
if (file("${params.AUGUSTUS_SPEC_DIR}/${params.AUGUSTUS_SPEC}").isDirectory()) {error "\"${params.AUGUSTUS_SPEC}\" is in ${params.AUGUSTUS_SPEC_DIR}!"} 

SCRATCH = params.SCRATCH
size = params.AUGUSTUS_TRAINING_SIZE
dict = 'Abinitio_predict' //Augustus_prediction以下にできる作業ディレクトリの名前

// --- Augustus ---
process augustus {

    cpus params.N_THREAD
    scratch "${SCRATCH}"
    stageOutMode 'copy'
    
    publishDir "${params.PDIR_PREP_ABINITIO_AUGUSTUS}", mode: 'copy', overwrite: false
//    publishDir "${params.PDIR_PREP_ABINITIO}", mode: 'copy', overwrite: false
    
    input:
    path trainingData from params.AUGUSTUS_TRAINING_DATA
    path genomeFASTA  from params.INPUT_GENOME
    path genomeRepOUT from params.INPUT_REPOUT
    
    output:
//    file "${CWD0}/PASS_information_for_snap_t${size}.txt" into augusutusOutput // To be ...
//    file "PASS_information_for_snap_t${size}.txt" into augusutusOutput // To be ...
    path "Augustus_prediction/${dict}" into augustusDictDir
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
    path trainingData       from params.SNAP_TRAINING_DATA // mapping.nfのアウトプット　mapping/to_learn_2nd/lemon_learn_2nd.gff3
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
