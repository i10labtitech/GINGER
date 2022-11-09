###--USAGE--------------------------------------------------------------
#実行コマンドでどのフェーズまで実行するのか指定する
if [ $# -lt 2 ]; then
    echo "$1 ...output file name $2 ...read1 $3 ...read2"
    exit 1
fi
FILE=$1
read1=$2
read2=$3

velveth $FILE 31 -fastq -short -separate $read1 $read2
velvetg $FILE -read_trkg yes 
oases $FILE
