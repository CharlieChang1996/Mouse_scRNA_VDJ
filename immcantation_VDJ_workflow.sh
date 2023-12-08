#!/bin/bash

# Immcantation VDJ analysis workflow
# Modified by hyzhang Thu Aug 24 11:52:38 UTC8+ 2023

###########################
### PREPARE ENVIRONMENT ###
###########################
#conda activate immcant-env

##########################
### DEFINE DIRECTORIES ###
##########################
# specify main project directory
main='/home/project7/hyzhang_2023/e007_mm_10x_multi/'


# !IMPORTANT!
# The Immcantation pipeline requires some accessory scripts that are not
# included in a nicely packaged conda or python package, but instead must be
# retrieved from their Bitbucket repository.
# All files in "https://bitbucket.org/kleinstein/immcantation/src/master/scripts/"
# should be downloaded to your local "$main/scripts/immcantation" directory.
PATH="$PATH:"$main"/immcantation"  # add immcantation scripts directory to the PATH


#############################################
### RETRIEVE IG BLAST REFERENCE DATABASES ###
#############################################
# Download reference databases
fetch_igblastdb.sh -o $main/data/immcantation/igblast
fetch_imgtdb.sh -o $main/data/immcantation/germlines/imgt

# Build IgBLAST database from IMGT reference sequences
# NOTE! This script uses the "realpath" command, which can be installed on MacOS with `brew install coreutils`
imgt2igblast.sh -i $main/data/immcantation/germlines/imgt -o $main/data/immcantation/igblast
cd $main/datadata/immcantation/igblast/database/
tar -xvf ncbi_human_c_genes.tar
cd $main
#######################################
### RUN IG BLAST ON VDJ FASTA FILES ###
#######################################
AssignGenes.py igblast \
-s $main/*/filtered_contig_*.fasta \
-b $main/data/immcantation/igblast \
--cdb ncbi_human_c_genes \
--organism mouse \
--loci ig \
--format blast


##############################################
### PARSE VDJ FILES INTO CHANGEO DB FORMAT ###
##############################################

# Create filtered VDJ seq database files for each sample
for sample in NP-gB NP-gH
do
    # Create tab-delimited database file to store seq alignment info
    MakeDb.py igblast \
    -i $main''$sample'_vdj_b/filtered_contig_'$sample'_igblast.fmt7' \
    -s $main''$sample'_vdj_b/filtered_contig_'$sample'.fasta' \
    -r $main'/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IG'*'.fasta' \
    -o $main''$sample'_vdj_b/filtered_contig_'$sample'_igblast_db-pass.tab' \
    --10x $main''$sample'_vdj_b/filtered_contig_annotations_'$sample'.csv' \
    --format changeo \
    --extended

    # Filter database to only keep functional sequences
    ParseDb.py select -d $main''$sample'_vdj_b/filtered_contig_'$sample'_igblast_db-pass.tab' \
    -f FUNCTIONAL \
    -u T \
    --outname $sample'_functional'

    # Parse database output into light and heavy chain change-o files
    ParseDb.py select -d $main''$sample'_vdj_b/'$sample'_functional_parse-select.tab' \
    -f LOCUS \
    -u "IGH" \
    --logic all \
    --regex \
    --outname $sample'_heavy'

    ParseDb.py select -d $main''$sample'_vdj_b/'$sample'_functional_parse-select.tab' \
    -f LOCUS \
    -u "IG[LK]" \
    --logic all \
    --regex \
    --outname $sample'_light'
done


#########################################################
### INFER IG GENOTYPE AND ESTIMATE SEQ DIST THRESHOLD ###
#########################################################
# Run R-script to infer sequence genotype and estimate VDJ seq hamming distance threshold
# Results are exported as "IGHV-genotyped_M#.tab", where "M#" is "M1", "M2", etc. for each mouse,
# and a fasta file of the V-segment germline sequences: "IGHV_genotype_M#.fasta".
# A .csv file of estimated threshold values will also be written to the output directory.

# The metadata.csv file should contain metadata about each of the BCR sequencing sample files;
# for an example of how the metadata file should be formatted, see the example metadata.csv
# file in the data/ subdirectory of the repository.
# Rscript $main/scripts/VDJ_analysis/01_VDJ_genotype_and_threshold.R \
# --VDJ_data_path $main'/data/VDJ_OTUs' \
# --metadata_file $main'/data/metadata.csv' \
# --germline_path $main'/data/immcantation/germlines' \
# --density_method 'density' \
# --default_threshold '0.1' \
# --output_path $main'/Analysis'


##################################################
### DEFINE CLONES AND GERMLINE SEQUENCES NP-gB ###
##################################################


    # Define clones (dist = distance threshold)
    # Output file is named "IGHV-genotyped_M#_clone-pass.tab"
    DefineClones.py -d $main'/Analysis/genotyping/IGHV-genotyped_NP-gB.tab' \
    --act set \
    --model ham \
    --norm len \
    --dist 0.13567 \
    --format changeo \
    --outname 'IGHV-genotyped_NP-gB' \
    --log $main'/Analysis/genotyping/IGHV-genotyped_NP-gB_DefineClones.log'

    # Create germline sequences using genotyped sequences from TIgGER
    # Output file is named "IGHV-genotyped_M#_germ-pass.tab"
    CreateGermlines.py -d $main'/Analysis/genotyping/IGHV-genotyped_NP-gB_clone-pass.tab' \
    -g dmask \
    --cloned \
    -r $main'/Analysis/genotyping/IGHV_genotype_NP-gB.fasta' \
    $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHD.fasta \
    $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHJ.fasta \
    --vf V_CALL_GENOTYPED \
    --format changeo \
    --outname 'IGHV-genotyped_NP-gB' \
    --log $main'/Analysis/genotyping/IGHV-genotyped_NP-gB_CreateGermlines.log'
##################################################
### DEFINE CLONES AND GERMLINE SEQUENCES NP-gH ###
##################################################

    # Define clones (dist = distance threshold)
    # Output file is named "IGHV-genotyped_M#_clone-pass.tab"
    DefineClones.py -d $main'/Analysis/genotyping/IGHV-genotyped_NP-gH.tab' \
    --act set \
    --model ham \
    --norm len \
    --dist 0.12727 \
    --format changeo \
    --outname 'IGHV-genotyped_NP-gH' \
    --log $main'/Analysis/genotyping/IGHV-genotyped_NP-gH_DefineClones.log'

    # Create germline sequences using genotyped sequences from TIgGER
    # Output file is named "IGHV-genotyped_M#_germ-pass.tab"
    CreateGermlines.py -d $main'/Analysis/genotyping/IGHV-genotyped_NP-gH_clone-pass.tab' \
    -g dmask \
    --cloned \
    -r $main'/Analysis/genotyping/IGHV_genotype_NP-gH.fasta' \
    $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHD.fasta \
    $main/data/immcantation/germlines/imgt/mouse/vdj/imgt_mouse_IGHJ.fasta \
    --vf V_CALL_GENOTYPED \
    --format changeo \
    --outname 'IGHV-genotyped_NP-gH' \
    --log $main'/Analysis/genotyping/IGHV-genotyped_NP-gH_CreateGermlines.log'




####################################
### QUANTIFY VDJ MUTATION BURDEN ###
####################################
#conda activate R-4.2
# exports results as ChangeO database file "VDJseq_mutation_quant.tab"
Rscript $main/VDJ_analysis/02_VDJ_mutation_quant.R \
--genotyped_path $main'/Analysis/genotyping' \
--output_path $main'/Analysis/mutation'


# Change conda environment to Sauron.v1
# conda deactivate
# source activate Sauron.v1


# #####################################
# ### ADD VDJ DATA TO SEURAT OBJECT ###
# #####################################
# Rscript $main/scripts/VDJ_analysis/03_VDJ_RNAseq_integration.R \
# --Seurat_object_path $main'/analysis/06_cluster/seurat_object.rds' \
# --changeo_db_path $main'/Analysis/mutation/VDJseq_mutation_quant.tab' \
# --output_path $main'/Analysis/seurat_object_VDJannot.rds'


# conda deactivate





