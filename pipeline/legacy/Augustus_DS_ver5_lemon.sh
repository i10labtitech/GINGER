###事前にインストール/パスが必要なツール
#AUGUSTUS

###----事前入力-------
spec= #パラメータのディレクトリ($AUGUSTUS_CONFIG_PATH/speciesの下に同じ名前が無いように気をつける)
dict="Abinitio_predict" #Augustus_prediction以下にできる作業ディレクトリの名前
pre_genome= #マスクしていないゲノム配列
mask_out= #RepeatMaskerの.outファイル
gff= #学習セット(配列の名前はfastaと対応している必要あり)
size=`echo 1000` #学習セットのサイズ
thread=16 #スレッド数
###-------------------
###----事前環境-------
ROOT=$(dirname `which $0`)
script=${ROOT}/script #用いるスクリプトのパス #/data/yuasa/script/abinitio/script
###-------------------
#マスクしたゲノム配列作成
g++ ${script}/simple_low_norepeatmask.cpp -std=c++0x -O3 
./a.out $pre_genome $mask_out ${pre_genome}.masked_without_SR_LC.fa
genome=${pre_genome}.masked_without_SR_LC.fa
rm ./a.out
########################

PWD1=`pwd`

#Augustusの出力ディレクトリの作成・移動
mkdir -p Augustus_prediction
cd Augustus_prediction

#作業ディレクトリの作成・移動
mkdir $dict
cd $dict

#SNAPに渡すファイルのディレクトリーを示したファイルの作成
PWD2=`pwd`
echo $PWD2 > ${PWD1}/PASS_information_for_snap_t${size}.txt

#gff3をgtfに変換
python ${script}/gff2gtf_better.py $gff > ./input.gtf
#学習に用いるCDSの情報のみにする
awk '$3=="CDS"' ./input.gtf > ./input2.gtf
#gtfをgbに変換(遺伝子周辺領域をmax1000bpとして切り出す)
gff2gbSmallDNA.pl input2.gtf $genome 1000 first.gb

#中間ファイル(gtf)の消去
rm input.gtf input2.gtf

#trainingしたパラメータを入れるディレクトリ作成
new_species.pl --species=$spec

#Augustusの入力にふさわしくない遺伝子情報のフィルタリング
#スタートコドンがなかったり、終始コドンが無いようなエラーを取り除く(genericはAugustusのマーカーセット)参考：Scipio
#--stopCodonExcludedFromCDS=trueはCDSにストップコドンを含めているか否かの設定がetrainingでできる。train.errファイルを確認しエラーが多く検出されていれば確認した方が良い。
etraining --species=generic first.gb 2> train.err
fgrep "gene" train.err | cut -f 2 -d " " > bad.etraining-test.lst
filterGenesOut_mRNAname.pl bad.etraining-test.lst first.gb > second.gb

#テスト用のファイル切り出し(テストファイルサイズ100)
randomSplit.pl second.gb 100

line_num=`grep "LOCUS" second.gb.train | wc -l`
num=`echo $(($line_num - $size))`

randomSplit.pl second.gb.train $num
	
#training1
etraining --species=$spec second.gb.train.train
	
#テスト
augustus --species=$spec second.gb.test | tee firsttest.out
grep -A 22 Evaluation firsttest.out > test_result1.txt &
	
#最適化
optimize_augustus.pl --species=$spec second.gb.train.train --cpus=${thread} --kfold=${thread}
	
#training2
etraining --species=$spec second.gb.train.train

#テスト2
augustus --species=$spec second.gb.test | tee second.out
grep -A 22 Evaluation second.out > test_result2.txt &
	
#fastaファイルの分割
python ${script}/pre_assembly_first_step.py $genome $thread
sort -nrk 1 -t "@" list.tMp > list_sort.tMp
wait
python ${script}/distribute_second_step.py list_sort.tMp $thread

for y in $(seq 1 $thread)
do
	python ${script}/tab_separate_third_step.py sPlitfile${y}.tMp > sPlitfile${y}.fasta &
done
wait

#分割したfastaファイルで予測
for k in $(seq 1 $thread)
do
        augustus --noInFrameStop=true --species=${spec} sPlitfile${k}.fasta  > abinitio${k}.gff &
done

wait

#分割して予測したファイルを統合
for l in $(seq 1 $thread)
do
        grep -v "#" abinitio${l}.gff > abinitio${l}_shaped.gff &
done

wait

#予測結果のファイルを統合(IDの付け方はアイソフォームを考慮していない)
cat abinitio*_shaped.gff > abinitio_shaped_cat.gff
python ${script}/re_name_for_AugustusGFF.py abinitio_shaped_cat.gff > Augustus_abinitio_tmp.gtf

#in-frame stopcodon geneを取り除く
${script}/inframe_stopcodon_exclude -i Augustus_abinitio_tmp.gtf -g $pre_genome -d gene -o Augustus_abinitio.gtf

rm abinitio_shaped_cat.gff Augustus_abinitio_tmp.gtf

#散らばったファイルの片付け
mkdir -p each_gff
mv *.gff each_gff
mkdir -p each_fasta
mv *.fasta each_fasta
rm *.tMp
#中間ファイルの削除(必要ならコメント以下をアウト)
rm -r each_gff
rm -r each_fasta
#更新日2020/10/06(湯淺)
#Augstusに--noInFrameStop=true追加
#変更日2019/11/19(湯淺)
#inframeの判定にマスクしてないゲノムを用いる
#変更日：2019/05/15(中村)
#in-frame stopcodonを含む遺伝子の除去
#作成日：2018/09/07(湯淺)
#最終結果のgff：Augustus_abinitio.gtf
#Augustusをインストールする必要あり(構築当時はVer.3.3を使用)
#Python2系を使用
#パラメーターディレクトリの関係もあるのでAugustusは各自インストールしてください。
#コンセプト：GenomeのマルチFASTAファイルを分割することで最後のgffの出力にかかる時間の削減を行う。(出力も結構時間かかる。)
#ファイルの分割はまず、配列が長い順に並べ替えた後、上から順番に指定した数のファイルに分けていく。
#この時の分け方は、ファイルがA,B,C,Dの4つとすると、A→B→C→D→D→C→B→A→A→B.....というような感じで分けてAとDの差が大きくならないようにしている。
