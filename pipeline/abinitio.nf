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
// if (file(params.INPUT_TOLEARN1ST).isFile())       {} else { error "No path : \"${params.TOLEARN1ST}\"" }
// if (file("${params.AUGUSTUS_SPEC_DIR}/${params.AUGUSTUS_SPEC}").isDirectory()) {error "\"${params.AUGUSTUS_SPEC}\" is in ${params.AUGUSTUS_SPEC_DIR}!"} 
if (file("${params.AUGUSTUS_WORK_DIR}").isDirectory()) { error "The working directory for Augustus already exists: \"${params.AUGUSTUS_WORK_DIR}\"" }

size = params.AUGUSTUS_TRAINING_SIZE
dict = 'Abinitio_predict'

SCRATCH = params.SCRATCH

// --- Augustus ---
process augustus {

    cpus params.N_THREAD
    scratch "${SCRATCH}"
    stageOutMode 'copy'
    
    publishDir "${params.PDIR_PREP_ABINITIO_AUGUSTUS}", mode: 'copy', overwrite: false
    
    input:
    path trainingData from params.AUGUSTUS_TRAINING_DATA
    path genomeFASTA  from params.INPUT_GENOME
    path genomeRepOUT from params.INPUT_REPOUT
    
    output:
    path "Augustus_prediction/${dict}" into augustusDictDir
    path "Augustus_abinitio.gtf"
    path "second.gb.train.train" into augustusTrainTrain
    
    shell:
    '''
    mkdir -p !{params.AUGUSTUS_WORK_DIR}
    cp -r !{params.AUGUSTUS_CONFIG_DIR} !{params.AUGUSTUS_WORK_DIR}

#    export AUGUSTUS_CONFIG_PATH="!{params.AUGUSTUS_CONFIG_DIR}"
    export AUGUSTUS_CONFIG_PATH="!{params.AUGUSTUS_WORK_DIR}/config"
    export PATH="$PATH:!{params.AUGUSTUS_SCRIPT_DIR}"
    export PATH="$PATH:!{params.AUGUSTUS_BIN_DIR}"

    script=!{params.UTILPATH_ABINITIO}

    ${script}/simple_low_norepeatmask !{genomeFASTA} !{genomeRepOUT} !{genomeFASTA}.masked_without_SR_LC.fa
    genome=!{genomeFASTA}.masked_without_SR_LC.fa
    
    CWD0=`pwd`

    mkdir -p Augustus_prediction
    cd Augustus_prediction

    mkdir !{dict}
    cd !{dict}

    CWD1=`pwd`

    python ${script}/gff2gtf_better.py ${CWD0}/!{trainingData} > ./input.gtf
    awk '$3=="CDS"' ./input.gtf > ./input2.gtf
    perl !{params.AUGUSTUS_SCRIPT_DIR}/gff2gbSmallDNA.pl input2.gtf ${CWD0}/${genome} 1000 first.gb
    
    perl !{params.AUGUSTUS_SCRIPT_DIR}/new_species.pl --species=!{params.AUGUSTUS_SPEC}
    
    !{params.ETRAINING} --species=generic first.gb 2> train.err
    fgrep "gene" train.err | cut -f 2 -d " " > bad.etraining-test.lst
    perl !{params.AUGUSTUS_SCRIPT_DIR}/filterGenesOut_mRNAname.pl bad.etraining-test.lst first.gb > second.gb
    
    perl !{params.AUGUSTUS_SCRIPT_DIR}/randomSplit.pl second.gb 100
    
    line_num=`grep "LOCUS" second.gb.train | wc -l`
    num=`echo $((${line_num} - !{size}))`
    
    perl !{params.AUGUSTUS_SCRIPT_DIR}/randomSplit.pl second.gb.train ${num}
	
    # Training 1
    !{params.ETRAINING} --species=!{params.AUGUSTUS_SPEC} second.gb.train.train
	
    # Test 1
    !{params.AUGUSTUS} --species=!{params.AUGUSTUS_SPEC} second.gb.test | tee firsttest.out
    grep -A 22 Evaluation firsttest.out > test_result1.txt &
	    
    # Optimization
    perl !{params.AUGUSTUS_SCRIPT_DIR}/optimize_augustus.pl --species=!{params.AUGUSTUS_SPEC} second.gb.train.train --cpus=!{params.N_THREAD} --kfold=!{params.N_THREAD}
	
    # Training 2
    !{params.ETRAINING} --species=!{params.AUGUSTUS_SPEC} second.gb.train.train
    
    # Test 2
    !{params.AUGUSTUS} --species=!{params.AUGUSTUS_SPEC} second.gb.test | tee second.out
    grep -A 22 Evaluation second.out > test_result2.txt &
	    
    #
    python ${script}/pre_assembly_first_step.py ${CWD0}/${genome} !{params.N_THREAD}
    sort -nrk 1 -t "@" list.tMp > list_sort.tMp
    wait
    python ${script}/distribute_second_step.py list_sort.tMp !{params.N_THREAD}
    
    for y in $(seq 1 !{params.N_THREAD})
    do
	    python ${script}/tab_separate_third_step.py sPlitfile${y}.tMp > sPlitfile${y}.fasta &
    done
    wait

    #
    for k in $(seq 1 !{params.N_THREAD})
    do
        !{params.AUGUSTUS} --noInFrameStop=true --species=!{params.AUGUSTUS_SPEC} sPlitfile${k}.fasta  > abinitio${k}.gff &
    done
    wait
    
    #
    for l in $(seq 1 !{params.N_THREAD})
    do
        grep -v "#" abinitio${l}.gff > abinitio${l}_shaped.gff &
    done
    wait

    #
    cat abinitio*_shaped.gff > abinitio_shaped_cat.gff
    python ${script}/re_name_for_AugustusGFF.py abinitio_shaped_cat.gff > Augustus_abinitio_tmp.gtf

    #
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
    
    input:
    path trainingData       from params.SNAP_TRAINING_DATA
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
    export PATH="!{params.SNAP_DIR}:$PATH"

    script=!{params.UTILPATH_ABINITIO}

    CWD0=`pwd`

    mkdir -p SNAP
    cd SNAP

    CWD1=`pwd`

    python ${script}/down_size.py ${CWD0}/!{augustusTrainTrain} ${CWD0}/!{trainingData} 
    num_locus=`grep "LOCUS" ${CWD0}/!{augustusTrainTrain} | wc -l`
    gff=down_size${num_locus}.gff

    mkdir snap_work
    cd snap_work
    
    mkdir -p make_data_set
    cd make_data_set
    
    awk '{print $1}' ${CWD1}/${gff} > genome.pre.txt
    python ${script}/same_line_remover.py genome.pre.txt > genome.txt

    python ${script}/rmk_fasta_maker_2018NOV11.py ${CWD0}/!{genomeMaskedFASTA}

    echo "##FASTA" > tmp.txt

    while read line
    do
	    python ${script}/first_gff_picker.py ${CWD1}/${gff} ${line} > ${line}.st.tmp.gff
		python ${script}/second_gff_liner.py ${line}.st.tmp.gff > ${line}.one_line.tmp.gff
		sort -nk 4 ${line}.one_line.tmp.gff > ${line}.sort.tmp.gff
		python ${script}/third_gff_reset.py ${line}.sort.tmp.gff > ${line}.en.tmp.gff
		python ${script}/fasta_pickup_2018NOV11.py ${CWD0}/!{genomeMaskedFASTA} ${line} > ${line}.fasta
        perl ${script}/gff2ann.pl ${line}.en.tmp.gff > ${line}.ann
        mv ${line}.fasta ${line}.dna
	done < ./genome.txt

    cd ../

    mkdir snap_prediction

    cat make_data_set/*.ann > snap_prediction/genome.ann
    cat make_data_set/*.dna > snap_prediction/genome.dna
    wait

    cd snap_prediction/

    !{params.FATHOM} genome.ann genome.dna -gene-stats > gene-stats.log 2>&1
    !{params.FATHOM} genome.ann genome.dna -validate > validate.log 2>&1
    !{params.FATHOM} genome.ann genome.dna -categorize 1000 > categorize.log 2>&1
    !{params.FATHOM} uni.ann uni.dna -export 1000 -plus> uni-plus.log 2>&1
    !{params.FORGE} export.ann export.dna > forge.log 2>&1
    perl !{params.SNAP_DIR}/hmm-assembler.pl !{params.SNAP_SPECIES} . > !{params.SNAP_SPECIES}.hmm
    !{params.SNAP} !{params.SNAP_SPECIES}.hmm ${CWD0}/!{genomeMaskedFASTA} -gff -aa !{params.SNAP_SPECIES}.faa > !{params.SNAP_SPECIES}.gff
    
    python ${script}/SNAP_to_GFF3.py !{params.SNAP_SPECIES}.gff > !{params.SNAP_SPECIES}.lemonfmt.gff
    
    python ${script}/frame_checker_and_inframe_eliminator.snap.lemonfmt.py ${CWD0}/!{genomeMaskedFASTA} !{params.SNAP_SPECIES}.faa !{params.SNAP_SPECIES}.lemonfmt.gff > !{params.SNAP_SPECIES}.lemonfmt.added_frame_infomation.gff
    sed '/^$/d' !{params.SNAP_SPECIES}.lemonfmt.added_frame_infomation.gff > !{params.SNAP_SPECIES}.pre_final.gff
    ${script}/makefasta -f -i !{params.SNAP_SPECIES}.pre_final.gff -g ${CWD0}/!{genomeFASTA} -o !{params.SNAP_SPECIES}.pre_final.fna
    python ${script}/terminal_exon_to_cds_trouble_fix_atSNAP.py !{params.SNAP_SPECIES}.pre_final.fna !{params.SNAP_SPECIES}.faa !{params.SNAP_SPECIES}.pre_final.gff > !{params.SNAP_SPECIES}.final.gff
    ##########################
    cd ../
    mkdir -p final_result
    cd final_result
    cp ../snap_prediction/!{params.SNAP_SPECIES}.gff ./!{params.SNAP_SPECIES}.original.gff
    cp ../snap_prediction/!{params.SNAP_SPECIES}.hmm ./!{params.SNAP_SPECIES}.hmm
    cp ../snap_prediction/!{params.SNAP_SPECIES}.faa ./!{params.SNAP_SPECIES}.faa
    cp ../snap_prediction/!{params.SNAP_SPECIES}.final.gff ./!{params.SNAP_SPECIES}.final.gff
    cp  ./!{params.SNAP_SPECIES}.final.gff ../../../!{params.SNAP_SPECIES}.final.gff
    '''
}
