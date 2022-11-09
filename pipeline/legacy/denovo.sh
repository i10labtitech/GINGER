### Requirements ###
#TransDecoder
#GMAP
#cd-hit
#seqkit

### Utility path ###
ROOT=/data2/fumiya/mytools/Lemon/newLemon_v2.1/Prep/denovo
Filtering=/data2/fumiya/hetero_annotate/denovo_rework/script/filtering.py
Gmap_native=${ROOT}/scripts/gmap_native_to_format_converter.pl
Cufflinks_gff3=${ROOT}/scripts/cufflinks_gtf_to_alignment_gff3.pl
Cufflinks_fasta=${ROOT}/scripts/cufflinks_gtf_genome_to_cdna_fasta.pl
Cdna_align=${ROOT}/scripts/cdna_alignment_orf_to_genome_orf_nolimit.pl
#Cdna_align=${ROOT}/scripts/cdna_alignment_orf_to_genome_orf_nolimit.pl

### Tool path (User-defined) ###
PYTHON=python
GMAP_BUILD=gmap_build
GMAP=gmap
CD_HIT_EST=cd-hit-est
SEQKIT=seqkit
LONGORFS=TransDecoder.LongOrfs
PREDICT=TransDecoder.Predict

### Usage ###
if test $# -ne 5 ; then
    echo "1: genome.fasta <STRING>"
    echo "2: trinity.fasta <STRING>"
    echo "3: oases.fa <STRING>"
    echo "4: prefix in output file name <STRING>"
    echo "5: the number of threads"
else
    genome=$1
    trinity=$2
    oases=$3
    t=$5

    ### Check file existence ###
    [ -f $genome ] || abort "cannot find $genome"
    [ -f $trinity ] || abort "cannot find $trinity"
    [ -f $oases ] || abort "cannot find $oases"
    
    ### sequences ID simplification ###
    touch ${4}_denovo.cds
    $SEQKIT seq -i -m 300 $trinity >> ${4}_denovo.cds
    $SEQKIT seq -i -m 300 $oases >> ${4}_denovo.cds

    ### de novo sequences compression ###
    cd-hit-est -i ${4}_denovo.cds -c 1.0 -o ${4}_denovo.cd.fasta -T $t -M 100000

    ### GMAP alignment ###
    gmap_build -D . -d ${4}_genome_gmap_DB $genome > ${4}_gmap_build.log 2>&1
    gmap -S -t $t -n 1 -D . -d ${4}_genome_gmap_DB ${4}_denovo.cd.fasta > ${4}_denovo.cd.fasta.gmap 2> ${4}_gmap.stderr
    
    ### filtering alignment result ###
    $PYTHON $Filtering ${4}_denovo.cd.fasta.gmap 95 > ${4}_denovo.cd.fasta.gmap.filtered

    ### format conversion for ORF prediction
    $Gmap_native ${4}_denovo.cd.fasta.gmap.filtered GTF > ${4}_denovo.cd.fasta.gtf
    $Cufflinks_gff3 ${4}_denovo.cd.fasta.gtf > ${4}_denovo.cd.fasta.gff3
    $Cufflinks_fasta ${4}_denovo.cd.fasta.gtf $genome >${4}_denovo_genome.fa

    ### ORF prediction ###
    $LONGORFS -t ${4}_denovo_genome.fa -m 90
    $PREDICT --no_refine_starts --single_best_only -t ${4}_denovo_genome.fa
    
    ### format conversion ###
    $Cdna_align ${4}_denovo_genome.fa.transdecoder.gff3 ${4}_denovo.cd.fasta.gff3 ${4}_denovo_genome.fa > ${4}_final.gff3
fi

