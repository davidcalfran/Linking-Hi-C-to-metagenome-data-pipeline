# Hi-C to metagenome data pipeline

Wastewater treatment plants (WWTP) operate with a natural community resembling that of natural systems (Cydzik-Kwiatkowska and Zielińska, 2016). The incoming wastewater is much more a human microbiome related and the industrial organisms are more related to a community. In the WWTP they meet and have the potential to exchange DNA between urban microbiome and natural microbiome. 

The study of the populations and mobile genetic elements (MGE) from activated sludge as well as other complex environments will try to clarify (i) the emission and fate of genetic fragments from industrial settings (ii) the presence of horizontal gene transfer and tracking of mobile elements between microorganisms present in wastewater plants.

So far, Hi-C is highly implemented in epigenetics and cancer studies on the biomedical field (Burton et al., 2014; Orlando et al., 2018). However, not that many studies have been published on the environmental field where DNA is constantly being released and can be exchange and transferred. However, there is an increase of papers being published in high impact journals (Stalder et al., 2019), thus showing the increasing interest on this field due to its consequences on human health and risk assessments development. 

## Objectives

-	Develop a Hi-C pipeline in bash using already published datasets from wastewater samples.
-	Link antibiotic resistant genes and mobile genetic elements to specific microorganisms. 
-	Obtain which microorganisms can be highly potential candidates to uptake, exchange and transfer targeted genes. 


## Before we start

- bwa is assumed to be installed in path
- samtools is assumed to be installed in path
- blastn is assumed to be installed in path
- CheckM is assumed to be installed in path
- Bin3C binaries are supposed to be downloaded and be saved in current folder as a folder named "bin3C"
- refseq.genomes.k21s1000.msh is supposed to be downloaded and be in current folder. Necessary for mash. It can be downloaded here: https://mash.readthedocs.io/en/latest/tutorials.html important to save is with the following name: refseq.genome.k21s1000.msh
- Python2 is assumed to be installed in path

## Run the code

```sh
bash Hic_MG_pipeline.sh assembly.fa HiC_1.fastq HiC_2.fastq  megares_database_v1.01.fasta enterobacteriaceae.fsa integrase_database.fa 1
```

Where:

- assembly.fa is your metagenome assembly generated with MEGAHIT, SPADES, etc
- HiC_1/2.fastq are your hic sequencing data.
- megares_database_v1.01.fasta is the antibiotic resistance genes database
- enterobacteriaceae.fsa is the plasmid database
- integrase_database.fa is the integrons database
- "1" is the number that starts running the script

# The idea behind the code

## Obtain data

Normally, metagenomics and Hi-C sequencing data are uploaded in NCBI-SRA. It will be needed to download and dump SRA datasets into fastq files. 

For doing that:

```sh
./fastq-dump -I --split-files  file.sra
```

## Quality control

Before starting,it is reccommended to check for data quality of the next sequencing reads obtained. For doing this, a software such as FastQC is reccommended. You can find the software in:

https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

If paired-end files are (i.e) 250 bp and you see that the last 20 nucleotides quality is really low, it would be recccommended to remove the last 20 of both of your paired end files. For doing that, you can use (for example) Trimmomatic:

http://www.usadellab.org/cms/?page=trimmomatic 

## Metagenome Assembly

For metagenome assembly, we used MEGAHIT (there are many others). The parameters were the following:

```
k_max: 99
min_contig_length: 2000
megahit_paramter_preset: null
```

An easy way to perform the assembly would be to install megahit in your computer:

https://github.com/voutcn/megahit

and then run the following command:

```sh
megahit -1 PE1_1.fq.gz -2 PE1_2.fq.gz -o megahit_out
```

One of the outputs on the folder will be the assembly:  MEGAHIT.assembly.fa

## Align Hi-C paired-end library to the MEGAHIT assembly using BWA mem

First of all we need to generate the index on the assembly:

```sh
bwa index MEGAHIT.assembly.fa
```
In the same folder, we will be able to align the Hi-C paired ends to the assembly. 

```sh
bwa mem -5SP [assembly.fasta] [fwd_hic.fastq] [rev_hic.fastq] > hic_mem.sam

samtools view -S -h -b -F 2316  hic_mem.sam > hic_mem.bam
```
We use bwa mem to align Hi-C data, with the -5, -S, and -P options. These options are not documented that well in the online bwa documentation, but they are documented in the usage of bwa mem. -5 is sometimes called “the Hi-C option” as it was designed to help the aligner handle the statistical properties of Hi-C libraries better, mainly by reducing the amount of secondary and alternate mappings the aligner makes as those cause Hi-C data to become ambiguous. The -S and -P options cause the aligner not to try to use assumptions about the reads that might be true for shotgun or mate pair libraries in an effort to rescue more reads.

For the next steps we will need the hic_mem.bam file sorted by name. For doing that:

```sh
samtools sort -n hic_mem.bam > hic_mem_name.bam
```

## Bin our assembly with the aligned Hi-C reads using Bin3C

To extract metagenome-assembled genomes (MAGs) from metagenomic data using Hi-C (binning), we will use bin3C:

https://github.com/cerebis/bin3C

It requires:

1.	The metagenome assembly in fasta format: MEGAHIT.assembly.fa
2.	The Hi-C dataset mapped against the assembly (links) in BAM format and sorted by names: hic_mem_name.bam

These binaries require considerable memory (>32 GB). Because of that, Bin3C had to be installed under a virtual environment in a google cloud virtual machine of more than 32 GB memory. Once installed, two steps were followed:

### Create a contact map for analysis

```sh
python2 ./bin3C mkmap -e MluCI -v MEGAHIT.assembly.fa hic_mem_name.bam bin3c_out
```

Where MluCI is the enzyme that it was used (can be changed). In bin3c_out, there will be the contact map needed for the binning procedure. 

### Cluster the resulting contact map into genome bins

```sh
python2 ./bin3C cluster -v bin3c_out/contact_map.p.gz bin3c_clust
```

this will basically (for our interest) generate a folder containing all the bins(.fna). These fasta files have to be checked for completeness and contamination. For doing that, we performed CheckM. 

## CheckM: binning quality control

This step is crucial in order to proceed with the analysis. It will give us how complete are the bins generated and how much contamination there is. From here, we will filter out the low quality bins: those whose completeness is below 70% and contamination above 15%. For doing it we will use the general CheckM pipeline:

https://github.com/Ecogenomics/CheckM

Where: 

```sh
checkm lineage_wf -x fa fasta_bins_folder fasta_bins_out
```

You can use the --reduced tree option to make it go faster. This will generate a table with all the clusters, number of genomes, number of gene makers, completeness (%) and contamination (%). We reccommend to sort that table by completeness (>70-80%) and contamination (<5-15%) depending on how restrictive you want to be and then re-run CheckM to get the quality graphical representation (otherwise it will not work). 

## Assign taxonomy using Mash

For downloading mash: 

https://mash.readthedocs.io/en/latest/

Mash is a pre-computed database that will assign taxonomy to the binning output. We will need to perform a mash analysis per bin generated and filtered out of Bin3C plus CheckM analysis. For doing it:

Inside the folder with the .fna files, copy the binaries of mash plus the refseq.genomes.k21s1000.ms and run:

```sh
for i in `ls CL*.fna`; do ./mash screen -w -p 4 refseq.genomes.k21s1000.msh $i > $i.tab; done
```

This will generate a .tab file per .fna file. Then we need to sort them out to obtain the best hit.

We put all the $i.tab generated before in a folder “tab” and inside the folder we do the following (maybe necessary to create a tab_sorted folder): 

```sh
for i in `ls CL*.fna.tab`; do sort -gr $i > tab_sorted/$i.tab; done
```

This will give us a tab file per bin with the top 10 (let's say) best hits. Then, we need to extract the top hit per bin and extract it into a table. 

## Extract mash sorted hits on a table

To extract the best hit from each table generated (one per bin), the following was done: 

```sh
for i in `ls CL*.fna.tab.tab`; do head -n 1 $i> extracted/table$i.tab ; done
```

Then we concatenate them (from the folder):

```sh
cat * > merged-file
```

Then we extract the names of the files with data:

```sh
Ls > names.txt
```

And combine files: names.txt with merged-file.txt to get the list of bins (name of bin) plus the best hit assigned by mash. 

## Extract mash sorted hits on a table

For doing this, we needed to download two databases and generate our own one. For the antibiotic resistance genes database, MEGARES database was used (Doster et al., 2020). For the plasmids, the PlasmidFinder database was used (Carattoli et al., 2014). For the integrases, a smaller database was created from INTEGRALL database (Moura et al., 2009), following what it was done in the reference work (Stalder et al., 2019). 

For annotation, BLASTn in command line was decided to be used. In order to select only the best hit per contig in bin and database, the following commands were used. 

First, we needed to generate the database used index with makeblastdb:

```sh
makeblastdb -in database.fa -dbtype nucl
```

Important to mention that the database is built by nucleotides and not aminoacids. Once we had the indexes generated, we proceeded to the annotation. We needed to merge all the bins before generated in a fasta file:


```sh
cat * > merged_fasta.fna  [in fasta folder]
```

The way of obtaining only the best hit per cluster contig; 

```sh
blastn -db megares_database_v1.01.fasta -query merged_fasta.fna -perc_identity 80 -outfmt 7 -out results.txt

cat results.txt |awk '/hits found/{getline;print}' | grep -v "#" > ARG_top_hits.txt
```

This will give us a text file with the best antibiotic resistance gene per cluster and bin:

```
CL0001_0176     Rif|CP00034.1|gene3741|Rifampin|Rifampin-resistant_beta
```

The same can be done with whatever database you want to check. 


## References

1. Burton, J. N., Liachko, I., Dunham, M. J., and Shendure, J. (2014). Species-Level Deconvolution of Metagenome Assemblies with Hi-C–Based Contact Probability Maps. G3&amp;#58; Genes|Genomes|Genetics 4, 1339–1346. doi:10.1534/g3.114.011825.
2. Carattoli, A., Zankari, E., Garciá-Fernández, A., Larsen, M. V., Lund, O., Villa, L., et al. (2014). In Silico detection and typing of plasmids using plasmidfinder and plasmid multilocus sequence typing. Antimicrob. Agents Chemother. 58, 3895–3903. doi:10.1128/AAC.02412-14.
3. Cydzik-Kwiatkowska, A., and Zielińska, M. (2016). Bacterial communities in full-scale wastewater treatment systems. World J. Microbiol. Biotechnol. 32, 1–8. doi:10.1007/s11274-016-2012-9.
4. Demaere, M. Z., and Darling, A. E. (2019). Bin3C: Exploiting Hi-C sequencing data to accurately resolve metagenome-assembled genomes. Genome Biol. 20, 1–16. doi:10.1186/s13059-019-1643-1.
5. Doster, E., Lakin, S. M., Dean, C. J., Wolfe, C., Young, J. G., Boucher, C., et al. (2020). MEGARes 2.0: a database for classification of antimicrobial drug, biocide and metal resistance determinants in metagenomic sequence data. Nucleic Acids Res. 48, D561–D569. doi:10.1093/nar/gkz1010.
6. Fraser, J., Williamson, I., Bickmore, W. A., and Dostie, J. (2015). An Overview of Genome Organization and How We Got There: from FISH to Hi-C. Microbiol. Mol. Biol. Rev. 79, 347–372. doi:10.1128/mmbr.00006-15.
7. Moura, A., Soares, M., Pereira, C., Leitão, N., Henriques, I., and Correia, A. (2009). INTEGRALL: A database and search engine for integrons, integrases and gene cassettes. Bioinformatics 25, 1096–1098. doi:10.1093/bioinformatics/btp105.
8. Orlando, G., Law, P. J., Cornish, A. J., Dobbins, S. E., Chubb, D., Broderick, P., et al. (2018). Promoter capture Hi-C-based identification of recurrent noncoding mutations in colorectal cancer. Nat. Genet. 50, 1375–1380. doi:10.1038/s41588-018-0211-z.
9. Rocha, P. P., Raviram, R., Bonneau, R., and Skok, J. A. (2016). Breaking TADs: insights into hierarchical genome organization, Epigenomics, Future Medicine. 7, 523–526. doi:10.2217/epi.15.25.Breaking.
10. Stalder, T., Press, M. O., Sullivan, S., Liachko, I., and Top, E. M. (2019). Linking the resistome and plasmidome to the microbiome. ISME J., 2437–2446. doi:10.1038/s41396-019-0446-4.
