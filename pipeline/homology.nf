#!/usr/bin/env nextflow

/*
Copyright (C) 2018 Itoh Laboratory, Tokyo Institute of Technology

This file is part of GINGER.

GINGER is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

GINGER is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with GINGER; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

if (file(params.PDIR).isDirectory()) {
} else {
    if (file(params.PDIR).isEmpty()) {
        print "mkdir \"${params.PDIR}/\"\n"
        file(params.PDIR).mkdir()
    }
}
if (file(params.PDIR_PREP_HOMOLOGY).isDirectory()) {
    error "The publish directory already exists: \"${params.PDIR_PREP_HOMOLOGY}\""
} else {
    if (file(params.PDIR_PREP).isEmpty()) {
        print "mkdir \"${params.PDIR_PREP}/\"\n"
        file(params.PDIR_PREP).mkdir()
    }
    if (file(params.PDIR_PREP_HOMOLOGY).isEmpty()) {
        print "mkdir \"${params.PDIR_PREP_HOMOLOGY}/\"\n"
        file(params.PDIR_PREP_HOMOLOGY).mkdir()
    }
}
if (file(params.INPUT_GENOME).isFile())    {} else { error "No path : \"${params.INPUT_GENOME}\"" }
if (file(params.INPUT_RNASEQR1).isFile())  {} else { error "No path : \"${params.INPUT_RNASEQR1}\"" }
if (file(params.INPUT_RNASEQR2).isFile())  {} else { error "No path : \"${params.INPUT_RNASEQR2}\"" }

SCRATCH = params.SCRATCH

hDataChPrefix  = Channel.from(params.HOMOLOGY_DATA['PREFIX'])
hDataChProtein = Channel.from(params.HOMOLOGY_DATA['PROTEIN'])
hDataChSpalndb = Channel.from(params.HOMOLOGY_DATA['SPALNDB'])

// --- Homology ---
process homology {

    cpus params.N_THREAD
    scratch "${SCRATCH}"

    input:
    path genomeFasta from params.INPUT_GENOME
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
    cp !{params.MAKBLK} . # used in makeidx.pl

    # Preprocessing
    !{params.UTILPATH_HOMOLOGY}/fastarepair !{genomeFasta} refer.mfa
    !{params.UTILPATH_HOMOLOGY}/fastarepair2 !{protein} relate.faa

    # Spaln
    export PATH=./:\$PATH
    export PATH=`!{params.UTILPATH_HOMOLOGY}/getPath.pl !{params.SPALN}`:\$PATH
    !{params.MAKEIDX} -ip refer.mfa

    # Exon identity calculation
    !{params.SPALN} -Q7 -O4 -ospalnresult_o4 -M -yS# -T!{spalndb} -yB# -yZ -t!{params.N_THREAD} -drefer relate.faa 
    mv spalnresult_o4 !{params.OPREFIX}_spalnresult_alignment.tsv
    
    python !{params.UTILPATH_HOMOLOGY}/o4_to_gff.py !{params.OPREFIX}_spalnresult_alignment.tsv !{prefix} > !{prefix}_spalnresult.gff
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
