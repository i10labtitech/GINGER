# GINGER: Gradual INtegrated GEne Reconstruction

## Description #################################################################

`GINGER` is a tool that is implemented an integrated method for gene
structure prediction in higher eukaryotes. RNA-Seq-based methods,
ab initio-based methods, homology-based methods are performed, and then 
gene structures are reconstructed via dynamic programming with appropriately
weighted and scored exon/intron/intergenic regions. Different prediction
processes and filtering criteria are applied to multi-exon and single-exon
genes. This pipeline is implemented using `Nextflow`.

## Web site ####################################################################

<https://github.com/i10labtitech/GINGER>

## Reference ###################################################################

Taniguchi T, Okuno M, Shinoda T, Kobayashi F, Takahashi K, Yuasa H, 
Nakamura Y, Tanaka H, Kajitani R, Itoh T. 
GINGER: an integrated method for high-accuracy prediction of gene structure 
in higher eukaryotes at the gene and exon level. 
DNA Res. 2023 Aug 1;30(4) 
<https://doi.org/10.1093/dnares/dsad017>

## Synopsis ####################################################################
### Requirements ###############################################################

`conda` (or alternative) packaging command. Only `linux-64` is supported.

### Inputs #####################################################################

* A file containing genome sequences (`FASTA format`)
* A file containing hard-masked genome sequences by RepeatMasker (`FASTA format`)
* A file containing repeat information by RepeatMasker (`***.out`)
* A file containing RNA-Seq reads 1 (`FASTQ format`)
* A file containing RNA-Seq reads 2 (`FASTQ format`)
* Files containing protein sequences of closely related species (`FASTA format`)

### Installation ###############################################################

```
git clone --depth 1 https://github.com/i10labtitech/GINGER
cd GINGER/rattler-build/
 conda create -n ginger-build -c conda-forge rattler-build
 conda activate ginger-build
  rattler-build build
 conda deactivate
 conda create -n ginger-run -c ./output -c bioconda -c conda-forge ginger
cd ../../
conda activate ginger-run
gingerInitCfg #./nextflow.config.template is generated.
```

* [Note] If you have modified something, `git commit` it before `rattler-build build`.

### Settings ###################################################################

```
cp nextflow.config.template nextflow.config
```

Edit `nextflow.config` properly. [Note] File or directory names should be set as full paths.
* Edit the the variables named `INPUT_GENOME, INPUT_MASKEDGENOME, INPUT_REPOUT, INPUT_RNASEQR1, INPUT_RNASEQR2, HOMOLOGY_DATA`.
* Edit the variable named `PDIR` (output directory name).
* Edit the variable named `SCRATCH` (path to scratch directory).
* Edit the variable named `N_THREAD` (number of threads used).
* Edit the variable named `MAX_MEMORY` (maximum memory size used).

### Run ########################################################################

```
nextflow run $(which ginger.nf) -c nextflow.config
ginger_phase0.sh nextflow.config
ginger_phase1_auto.sh nextflow.config > phase1.log
ginger_phase2.sh 100 #minimum CDS length. 100 is just an example.
```

* [Note] An automatically calculated threshold for the score is written 
  like `score threshold = 1.25` in `phase1.log`.

To set a threshold for reconstructed gene structure's scores, 
use `ginger_phase1_manual.sh` with an extra argument for the threshold 
instead of `ginger_phase1_auto.sh`.

```
rm -r ginger_phase1.gff ginger_phase1_result/
ginger_phase1_manual.sh nextflow.config 1.0 #1.0 is just an example.
ginger_phase2.sh 50 #50 is just an example.
```

Finaly, you may summarize the result.

```
ginger_summary.sh nextflow.config
```

### Final output ###############################################################

The final outputs:
* `ginger_phase2.gff` : gene structures by GINGER (GFF3).
  [Note] See http://gmod.org/wiki/GFF3 for details.
* `ginger.pep`        : protein sequences of the gene structurs (FASTA)
* `ginger.cds`        : CDS of the gene structures (FASTA)
* `ginger_stats.tsv`  : statistical information of gene structures

## Obtaining test inputs (if you need) #########################################

```
mkdir sampleDir

## download
curl -L https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Caenorhabditis_elegans/all_assembly_versions/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz -o sampleDir/GCF_000002985.6_WBcel235_genomic.fna.gz
curl -L https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Caenorhabditis_briggsae/all_assembly_versions/GCF_000004555.2_CB4/GCF_000004555.2_CB4_translated_cds.faa.gz -o sampleDir/GCF_000004555.2_CB4_translated_cds.faa.gz
curl -L https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/180/635/GCA_000180635.4_El_Paco_v._4/GCA_000180635.4_El_Paco_v._4_translated_cds.faa.gz -o sampleDir/GCA_000180635.4_El_Paco_v._4_translated_cds.faa.gz
prefetch SRR5849934 -O sampleDir

## preparation
fastq-dump --split-files sampleDir/SRR5849934/SRR5849934.sra --outdir sampleDir
rm -r sampleDir/SRR5849934/
zcat sampleDir/GCF_000002985.6_WBcel235_genomic.fna.gz | seqkit seq -i > sampleDir/GCF_000002985.6_WBcel235_genomic.commentModified.fna
rm sampleDir/GCF_000002985.6_WBcel235_genomic.fna.gz
perl -pe 'tr/[a-z]/N/ unless /^>/' sampleDir/GCF_000002985.6_WBcel235_genomic.commentModified.fna > sampleDir/GCF_000002985.6_WBcel235_genomic.commentModified.masked.fna
touch sampleDir/GCF_000002985.6_WBcel235_genomic.out #dummy RepeatMasker output for test purpose.
gunzip sampleDir/GCF_000004555.2_CB4_translated_cds.faa.gz
gunzip sampleDir/GCA_000180635.4_El_Paco_v._4_translated_cds.faa.gz

## comparison after done
curl -L https://ftp.ncbi.nlm.nih.gov/genomes/refseq/invertebrate/Caenorhabditis_elegans/all_assembly_versions/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.gff.gz -o sampleDir/GCF_000002985.6_WBcel235_genomic.gff.gz
gunzip sampleDir/GCF_000002985.6_WBcel235_genomic.gff.gz
ginger_evaluate.sh sampleDir/GCF_000002985.6_WBcel235_genomic.gff ginger_phase2.gff mRNA CDS mRNA CDS
```

* [Note] On setting `nextflow.config`, delete the items for the third species in `HOMOLOGY_DATA` and set `SPALNDB` as `NematodC`.
