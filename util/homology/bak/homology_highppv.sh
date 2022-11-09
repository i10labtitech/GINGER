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
if [ $# -ne 3 ]; then
	echo "error::実行するには3個の引数が必要です。" 1>&2
	echo "./homology_highppv.sh [spalnresult_filter.gff3] [alignmentresult.txt] [prefix]" 1>&2
	exit 1
fi

cp $SCRIPT/function.hpp .

#fileが存在するかどうか
[ -f $1 ] || abort "$1 file not exist."
[ -f $2 ] || abort "$2 file not exist."

outputspalnfilter=$3"_spalnresult_over90.gff"

echo "start homology"
echo "$0 $1 $2 $3"

#filter開始
echo "start filtering"

awk '$2<90' $2 | awk '{print $1}' > alignmentfilterlist.txt

g++ $SCRIPT/exclude.cpp -std=c++0x -O3 -o exclude
./exclude alignmentfilterlist.txt $1 $outputspalnfilter

#削除
rm alignmentfilterlist.txt function.hpp exclude

echo "finish_homology"
