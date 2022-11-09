###--USAGE--------------------------------------------------------------
#実行コマンドでどのフェーズまで実行するのか指定する
if [ $# -lt 4 ]; then
    echo "$1 ...output file name $2 ...read1 $3 ...read2 $4 ... thread $5 ...memory"
    exit 1
fi
FILE=$1
read1=$2
read2=$3
T=$4
memory=$5

# ---- Trinity part ----
#trinityでassemble
mkdir $1
# ---- 1. Trinityでde novo transcriptome
time Trinity --seqType fq --left $2 --right $3 --output $1 --CPU $T --max_memory $memory --full_cleanup
