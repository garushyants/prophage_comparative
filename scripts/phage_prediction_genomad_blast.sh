#!/bin/bash
#$ -v SGE_FACILITIES
#$ -P unified
#$ -l osverfull="(7*|8*)"
#$ -l h_rt=259200,h_vmem=48G,mem_free=48G
#$ -pe multicore 8

source ncbi-facilities-load.sh
source ~/.bashrc

IFS=$'\n' alltasks=($(cat test-tasks.txt))
task=${alltasks[$SGE_TASK_ID-1]}
#task=${alltasks[0]}
echo $task
outprefix=$(echo ${task} | sed 's/.fasta.gz//')
pathToGenomesFolder="../test_genomes/"
blastdb="../20231002_34kecoli_blastdb/ecoli34k" #provide a path to the blast database for all E.coli genomes
outblast=${outprefix}_blastn.tsv

covfile=${outprefix}_nucl_cov_megablast.txt
 
 
#######BLAST analysis
maxtargetseq=10000 
identitycutoff=95

zcat ${pathToGenomesFolder}${task} | blastn -db $blastdb -task megablast -evalue 1e-5 -perc_identity $identitycutoff \
-max_target_seqs $maxtargetseq -out $outblast \
-outfmt "7 qacc qlen qstart qend qcovhsp sacc slen sstart send evalue ppos bitscore" -num_threads 8

echo "Blastn search complete"

#######get coverage from tabulated blastn output
hsplengthcutoff=300
evaluecutoff=10e-10
identitycutoffcov=98
mamba activate base
../scripts/blastn_tabular_to_coverage.py -b $outblast -l $hsplengthcutoff -e $evaluecutoff -i $identitycutoffcov
mamba deactivate

echo "Obtained coverage from blast output"

######geNomad
mamba activate genomad
genomad end-to-end --threads 8 --conservative --enable-score-calibration \
 --cleanup ${pathToGenomesFolder}${task} ${outprefix}_cons_fdr ../../genomad_db
mamba deactivate
####checkv
mamba activate vs2
checkv end_to_end ${outprefix}_cons_fdr/${outprefix}_summary/${outprefix}_virus.fna ${outprefix}_cons_fdr_checkv -t 8 -d ../../CheckVDb/checkv-db-v1.5

echo "geNomad and checkv complete"

####analyze in R and save outputs
Rscript ../scripts/vizualize_genomad_blast_predictions_with_phage_info.R -g ${outprefix} -w 4000 -u 0.2

echo "Analyzing results"
