# ./homology.sh [genome] [related species protein] [spalndatabase] [thread] [relate or not] [prefix]

#homology.shを動かすのに必要なscriptの位置
SCRIPT=/data/tkazuki/annotation/homologypipeline/homologyscript/

#コマンドが失敗した時にエラー出力して強制終了する関数
function abort
{
	echo "$@" 1>&2
	exit 1
}

#正しい引数かどうか
if [ $# -ne 6 ]; then
	echo "error::実行するには6個の引数が必要です。" 1>&2
	echo "./homology.sh [genome] [related species protein] [spalndatabase] [thread] [超近縁 -> 1,　近縁 -> 0] [prefix]" 1>&2
	exit 1
fi

#fileが存在するかどうか
[ -f $1 ] || abort "$1 file not exist."
[ -f $2 ] || abort "$2 file not exist."

outputspaln=$6"_spalnresult.gff"
outputspalnfilter=$6"_spalnresult_filter.gff"

#上書き阻止
[ ! -f refer.mfa ] || abort "refer.mfa exist. change name from refer.mfa, or overwrite"
[ ! -f relate.faa ] || abort "relate.faa exist. change name from relate.faa, or overwrite"
[ ! -f $outputspaln ] || abort "$outputspaln exist. change name from $outputspaln, or overwrite"
[ ! -f $outputspalnfilter ] || abort "$outputspalnfilter exist. change name from $outputspalnfilter, or overwrite"
[ ! -d spalndata ] || abort "spalndata directory exist. change name from spalndata, or overwrite"
[ ! -d spalndata2 ] || abort "spalndata2 directory exist. change name from spalndata2, or overwrite"

echo "start homology"
echo "$0 $1 $2 $3 $4 $5 $6"

#spaln動かすためスクリプトの用意
cp $SCRIPT/Makefile .
cp $SCRIPT/catchr.pl .
cp $SCRIPT/makblk.pl .
cp $SCRIPT/makeidx.pl .
cp $SCRIPT/function.hpp .

echo "start spaln"

#spalnを動かす前処理
g++ $SCRIPT/fastarepair.cpp -std=c++0x -O3 -o fastarepair
./fastarepair $1 refer.mfa

g++ $SCRIPT/fastarepair2.cpp -std=c++0x -O3 -o fastarepair
./fastarepair $2 relate.faa

#spalnを動かす
./makeidx.pl -ip refer.mfa
spaln -Q7 -O0 -ospalnresult -M -yS# -T$3 -yB# -yZ# -t$4 -drefer relate.faa

#tmpファイルの削除
rm Makefile catchr.pl makblk.pl makeidx.pl fastarepair

#spaln結果の修正
g++ $SCRIPT/spaln_repair2.cpp -std=c++0x -O3 -o spalnrepair
./spalnrepair spalnresult refer.mfa $outputspaln

#削除
rm spalnresult spalnrepair refer.bkp refer.ent refer.grp refer.idx refer.seq refer.odr

echo "start alignment"
#alignment開始
#alignment準備
g++ $SCRIPT/gffmakecds.cpp -std=c++0x -O3 -o gffmakecds
./gffmakecds $outputspaln refer.mfa prefilter00.fa
g++ $SCRIPT/fastaprotchange.cpp -std=c++0x -O3 -o fastaprotchange
./fastaprotchange prefilter00.fa prefilter01.fa
g++ $SCRIPT/fastarepair2.cpp -std=c++0x -O3 -o fastarepair
./fastarepair prefilter01.fa prefilter.fa
g++ $SCRIPT/fasta_grep.cpp -std=c++0x -O3 -o fastagrep
./fastagrep relate.faa $outputspaln prefilter2.fa

#削除
rm gffmakecds fastaprotchange fastagrep prefilter00.fa refer.mfa fastarepair relate.faa

#分割
mkdir spalndata
mkdir spalndata2
cd spalndata
g++ $SCRIPT/split.cpp -std=c++0x -O3 -o split
./split ../prefilter.fa 2
cd ../spalndata2
g++ $SCRIPT/split.cpp -std=c++0x -O3 -o split
./split ../prefilter2.fa 2
cd ..

#alignment
wcl=`wc -l prefilter.fa | awk '{print $1}'`
wcl2=`expr $wcl / 2 - 1`

touch alignmentresult.txt
g++ $SCRIPT/mafftalignment.cpp -std=c++0x -O3 -o mafftalignment
for i in `seq 0 $wcl2`
do
	cat spalndata/$i spalndata2/$i > premafft.fa
	mafft --auto premafft.fa > mafftresult.txt
	./mafftalignment mafftresult.txt maffttmp.txt
	cat alignmentresult.txt maffttmp.txt > maffttmp2.txt
	mv maffttmp2.txt alignmentresult.txt
done

#削除
rm premafft.fa mafftresult.txt maffttmp.txt mafftalignment prefilter.fa prefilter2.fa
rm -rf spalndata spalndata2

#filter開始
echo "start filtering"

#flameshift遺伝子除去
g++ $SCRIPT/flameshiftgrep.cpp -std=c++0x -O3 -o flameshiftfilter
./flameshiftfilter prefilter01.fa flameshiftlist.txt
g++ $SCRIPT/exclude.cpp -std=c++0x -O3 -o exclude
./exclude flameshiftlist.txt $outputspaln spaln_filter_flameshift.gff

#低identity遺伝子除去
if [ $5 = 0 ]; then
	awk '$2<30' alignmentresult.txt | awk '{print $1}' > alignmentfilterlist.txt
else
	awk '$2<50' alignmentresult.txt | awk '{print $1}' > alignmentfilterlist.txt
fi


./exclude alignmentfilterlist.txt spaln_filter_flameshift.gff $outputspalnfilter

#削除
rm flameshiftlist.txt alignmentfilterlist.txt flameshiftfilter prefilter01.fa 
rm exclude spaln_filter_flameshift.gff function.hpp

echo "finish_homology"
