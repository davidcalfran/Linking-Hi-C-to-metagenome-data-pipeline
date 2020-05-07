#!/bin/bash

set -e

#############
# Variables #
#############

MASH_PATH="./mash"
BIN3C_PATH="./bin3C/bin3C.py"

#samtools is assumed to be installed in path
#blastn is assumed to be installed in path
#CheckM is assumed to be installed in path
#Bin3C is supposed to be downloaded and binaries be all in folder callled "bin3C"
#refseq.genomes.k21s1000.msh is supposed to be downloaded and be in current folder. Necessary for mash.

###################
# PARSE ARGUMENTS #
###################

if [[ -z $1 ]]; then
    echo "first argument should be an ASSEMBLY in fasta"
    exit 1
fi

if [[ -z $2 ]]; then
    echo "second argument should be Hi-C PairedEnd file 1"
    exit 1 #0 means script was successful, 1 means script failed
fi

if [[ -z $3 ]]; then
echo "third argument should be Hi-C PairedEnd file 2"
exit 1 #0 means script was successful, 1 means script failed
fi

if [[ -z $4 ]]; then
echo "database for ARGs"
exit 1
fi

if [[ -z $5 ]]; then
echo "database for Plasmids"
exit 1
fi

if [[ -z $6 ]]; then
echo "database for Integrons"
exit 1
fi

if [[ -z $7 ]]; then
echo "give me a number!"
exit 1 #0 means script was successful, 1 means script failed
fi

###Align Hi-C paired-end library to the MEGAHIT assembly using BWA mem

if [[ $start < 2 ]]; then #if booth documents are there
    # STEP1: Generate index file from assembly
    bwa index $1
fi

if [[ $start < 2 ]]; then #if booth documents are there
    # STEP2: aligh Hi-C to Assembly to generate the .sam file
    bwa mem -5SP $1 $2 $3 > hic_mem.sam
fi

if [[ $start < 2 ]]; then
# STEP3: Generate the .bam file
    samtools view -S -h -b -F 2316 hic_mem.sam > hic_mem.bam
fi

if [[ $start < 2 ]]; then
# STEP4: Sort by name the bam file
samtools sort -n hic_mem.bam > hic_mem_name.bam
fi

####Bin our assembly with the aligned Hi-C reads using Bin3C

if [[ $start < 2 ]]; then
# STEP5: Create a contact map for analysis
python2 $BIN3C_PATH mkmap -e MluCI -v $1 hic_mem_name.bam bin3c_out
fi

if [[ $start < 2 ]]; then
# STEP6: Cluster the resulting contact map into genome bins. Only in ubuntu environments with more than 32GB ram
python2 $BIN3C_PATH cluster -v bin3c_out/contact_map.p.gz bin3c_clust
fi

###Binning Quality Control with CheckM

if [[ $start < 2 ]]; then
# STEP 7: Generate a folder with the quality control about the binning process. It can be used for "removing" low quality bins from folder and start over STEP 8
checkm lineage_wf -x fna bin3c_clust/fasta fasta_bins_out
fi

###Assign taxonomy using Mash and get the mash_sorted_clean.txt file

if [[ $start < 2 ]]; then
# STEP 8: Copy all fasta files in current folder
   cp bin3c_clust/fasta/* .
fi

if [[ $start < 2 ]]; then
# STEP 9: Generate a file with hits
for i in `ls CL*.fna`; do $MASH_PATH screen -w -p 4 refseq.genomes.k21s1000.msh $i > $i.tab; done
fi

if [[ $start < 2 ]]; then
# STEP 10: Generate a folder with the tab_sorted mash files
mkdir tab_sorted
fi

if [[ $start < 2 ]]; then
# STEP 11: Sort mash results to have on top the best hit
for i in `ls CL*.fna.tab`; do sort -gr $i > tab_sorted/$i.tab; done
fi

if [[ $start < 2 ]]; then
# STEP 12: Copy all sorted files to current folder
   cp tab_sorted/* .
fi

if [[ $start < 2 ]]; then
# STEP 13: Extract mash sorted hits on a table
for i in `ls CL*.fna.tab.tab`; do head -n 1 $i> $i.tab ; done
fi

if [[ $start < 2 ]]; then
# STEP 14: Extract mash sorted hits on a table.
bash concat.sh
fi

if [[ $start < 2 ]]; then
# STEP 15: Extract names of fasta files (CL....fna) This is ok.
    ls bin3c_clust/fasta/ | sed 's/.fna//g' > names.txt
fi

if [[ $start < 2 ]]; then
# STEP 16: Merge both names file and merged_tophits_mash.tab file in one. Then we select only the columns with the CL... code and the name of the species to have a clean mash results file.
    paste -d ' ' names.txt merged_tophits_mash.tab > MASH_results.txt
fi


### Generate ARGs/Integrons/Plasmids hits on our hits

if [[ $start < 2 ]]; then
# STEP 17: Merge all the bins before generate in a fasta file
 cat bin3c_clust/fasta/* > merged_fasta.fna
fi

if [[ $start < 2 ]]; then
# STEP 18: Index file the databases that are going to be used for getting the hits
makeblastdb -in $4 -dbtype nucl
makeblastdb -in $5 -dbtype nucl
makeblastdb -in $6 -dbtype nucl

fi

if [[ $start < 2 ]]; then
# STEP 19: Perform the blastn with the merged fasta file and the indexed databases

blastn -db $4 -query merged_fasta.fna -evalue 1e-20 -outfmt 7 -out ARG_results.txt
cat ARG_results.txt |awk '/hits found/{getline;print}' | grep -v "#" | awk '{print $1,$2}' > ARG_top_hits.txt
blastn -db $5 -query merged_fasta.fna -evalue 1e-20 -outfmt 7 -out plasmid_results.txt
cat plasmid_results.txt |awk '/hits found/{getline;print}' | grep -v "#" | awk '{print $1,$2}'  > plasmid_top_hits.txt
blastn -db $6 -query merged_fasta.fna -evalue 1e-20 -outfmt 7 -out integrases_results.txt
cat integrases_results.txt |awk '/hits found/{getline;print}' | grep -v "#" | awk '{print $1,$2}' > integrases_top_hits.txt


fi

###Generate the file showing the links between aligned contigs between the assembly and the Hi-C reads

if [[ $start < 2 ]]; then
# STEP 20: Generate alignment to links file necessary for STEP 21 scripts.
samtools view hic_mem_name.bam | awk '{hash[$3"\t"$7]++}END{for (x in hash) {print x"\t"hash[x]/2}}'> contig_links.txt
fi

###Generate sum of interactions between ARGs/Integrases/Plasmids and its respective clusters

if [[ $start < 2 ]]; then
# STEP 21: Generate sum of interaction files
python2 script_ARG_d.py
python2 script_plasmid_d.py
python2 script_integrases_d.py

fi

###Prepare files for heatmap generation

if [[ $start < 2 ]]; then
# STEP 22: Clean files and join sum of interactions with the events
sed 's/,/ /g' result-sum_ARG.txt | tail -n +2 > result-sum_ARG_clean.txt
sed 's/,/ /g' result-sum_plasmid.txt | tail -n +2 > result-sum_plasmid_clean.txt
sed 's/,/ /g' result-sum_integrases.txt | tail -n +2 > result-sum_integrases_clean.txt

#First we remove spaces from script results to generate columns, then we remove header and then we join with top hits of the specific event.
fi

if [[ $start < 2 ]]; then
# STEP 23: Join them with top hits from blast results

join -1 1 -2 1 ARG_top_hits.txt result-sum_ARG_clean.txt  > ARG_hits_suminteractions.txt
join -1 1 -2 1 plasmid_top_hits.txt result-sum_plasmid_clean.txt  > plasmid_hits_suminteractions.txt
join -1 1 -2 1 integrases_top_hits.txt result-sum_integrases_clean.txt  > integrases_hits_suminteractions.txt
fi

if [[ $start < 2 ]]; then
# STEP 24: Remove "_" from ARG_hits_suminteractiosn so the species from MASH can be added and join them
awk '{gsub("_"," ",$1)}1' ARG_hits_suminteractions.txt > ARG_hits_suminteractions_clean.txt
join -1 1 -2 1 ARG_hits_suminteractions_clean.txt MASH_results.txt > ARG_final_results_R.txt
awk '{gsub("_"," ",$1)}1' plasmid_hits_suminteractions.txt > plasmid_hits_suminteractions_clean.txt
join -1 1 -2 1 plasmid_hits_suminteractions_clean.txt MASH_results.txt > plasmid_final_results_R.txt
awk '{gsub("_"," ",$1)}1' integrases_hits_suminteractions.txt > integrases_hits_suminteractions_clean.txt
join -1 1 -2 1 integrases_hits_suminteractions_clean.txt MASH_results.txt > integrases_final_results_R.txt
fi

if [[ $start < 2 ]]; then
# STEP 24: Clean file part I
awk '{print $1,$3,$5,$14,$15}' ARG_final_results_R.txt > ARG_final_results_R_clean.txt
awk '{print $1,$3,$5,$14,$15}' plasmid_final_results_R.txt > plasmid_final_results_R_clean.txt
awk '{print $1,$3,$5,$14,$15}' integrases_final_results_R.txt > integrases_final_results_R_clean.txt
fi

if [[ $start < 2 ]]; then
# STEP 25: Clean file part II: substract column by column and generate the final file
grep "CL" ARG_final_results_R.txt | awk 'BEGIN{FS=".fna.gz"}; {print $2}' | awk 'BEGIN{FS=".1 "}; {print $2}' |awk '{print $1,$2}' > ARG_species.txt #NO SE COMO LIMPIAR CORRECTAMENTE ESTO. PREGUNTAR A ANDREU
awk '{print $4,$5}' plasmid_final_results_R_clean.txt > plasmid_species.txt
awk '{print $4,$5}' integrases_final_results_R_clean.txt > integrases_species.txt
awk '{print $1}' ARG_final_results_R_clean.txt > ARG_cluster.txt
awk '{print $1}' plasmid_final_results_R_clean.txt > plasmid_cluster.txt
awk '{print $1}' integrases_final_results_R_clean.txt > integrases_cluster.txt
grep "CL" ARG_final_results_R.txt | awk 'BEGIN{FS="|"}{split($0,a,"|") ; if (a[length(a)]~"Requires") {print a[length(a)-3]} else {print a[length(a)-2]}}' > ARG_event.txt
awk '{print $2}' ARG_final_results_R_clean.txt | awk 'BEGIN{FS="|"}{split($0,a,"|") ; if (a[length(a)]~"Requires") {print a[length(a)-1]} else {print a[length(a)-0]}}' > ARG_gene.txt
awk '{print $2}' plasmid_final_results_R_clean.txt |sed 's/_.*//' > plasmid_event.txt
awk '{print $2}' integrases_final_results_R_clean.txt > integrases_event.txt
awk '{print $3}' ARG_final_results_R_clean.txt > ARG_sumint.txt
awk '{print $3}' plasmid_final_results_R_clean.txt > plasmid_sumint.txt
awk '{print $3}' integrases_final_results_R_clean.txt > integrases_sumint.txt
fi

#Normalization

if [[ $start < 2 ]]; then
# STEP 26: Clean file part IV: substract number of contigs per cluster and normalize sum of interactions
sed 's/,/ /g' result_CLXXXXOcurrences_ARG.txt | tail -n +2 > Ocurrences_ARG.txt
join -1 1 -2 1 ARG_cluster.txt Ocurrences_ARG.txt | awk {'print $2'} > Ocurrences_ARG_clean.txt
paste ARG_sumint.txt Ocurrences_ARG_clean.txt | awk '{print($1/$2)}' > normalizedsum_ARG.txt

sed 's/,/ /g' result_CLXXXXOcurrences_plasmid.txt | tail -n +2 > Ocurrences_plasmid.txt
join -1 1 -2 1 plasmid_cluster.txt Ocurrences_plasmid.txt | awk {'print $2'} > Ocurrences_plasmid_clean.txt
paste plasmid_sumint.txt Ocurrences_plasmid_clean.txt | awk '{print($1/$2)}' > normalizedsum_plasmid.txt

sed 's/,/ /g' result_CLXXXXOcurrences_integrases.txt | tail -n +2 > Ocurrences_integrases.txt
join -1 1 -2 1 integrases_cluster.txt Ocurrences_integrases.txt | awk {'print $2'} > Ocurrences_integrases_clean.txt
paste integrases_sumint.txt Ocurrences_integrases_clean.txt | awk '{print($1/$2)}' > normalizedsum_integrases.txt
fi


if [[ $start < 2 ]]; then
# STEP 27: Clean file part III: join the parts and sort them by second column (Species)
paste ARG_cluster.txt ARG_species.txt ARG_event.txt ARG_gene.txt ARG_sumint.txt Ocurrences_ARG_clean.txt normalizedsum_ARG.txt | sort -k 2 > ARG_R.txt

paste plasmid_cluster.txt plasmid_species.txt plasmid_event.txt plasmid_sumint.txt Ocurrences_plasmid_clean.txt normalizedsum_plasmid.txt | sort -k 2   > plasmid_R.txt
paste integrases_cluster.txt integrases_species.txt integrases_event.txt integrases_sumint.txt Ocurrences_integrases_clean.txt normalizedsum_integrases.txt | sort -k 2 > integrases_R.txt

fi

if [[ $start < 2 ]]; then
# STEP 28: Add header
echo -e 'CLUSTER\tSPECIES\tARG family\tARG\tSUM\tNumber of Contigs\tNormalizedSum'| cat - ARG_R.txt  > ARG_R_clean.csv
echo -e 'CLUSTER\tSPECIES\tEVENT\tSUM\tNumber of Contigs\tNormalizedSum'| cat - plasmid_R.txt > plasmid_R_clean.csv
echo -e 'CLUSTER\tSPECIES\tEVENT\tSUM\tNumber of Contigs\tNormalizedSum'| cat - integrases_R.txt > integrases_R_clean.csv

fi


if [[ $start < 2 ]]; then
# STEP 28: Run R

Rscript hic.r

fi
