

# PlantChimera
PlantChimera aims to provide a streamlined, efficient, and accurate method for identifying fusion transcripts in plant transcriptomic datasets. It employs a highly configurable pipeline with stringent alignment, filtering to detect fusion events from RNA-Seq data with high confidence. By integrating multiple layers of filtering, from sequence overlap and gene orientation to entropy and junction distance, the pipeline reduces the likelihood of false positives while retaining biologically significant fusion candidates. The source code can be found at [https://github.com/rebminso/PlantChimera.git](https://github.com/rebminso/PlantChimera.git) and is open source and compatible with multiple platforms. Installation instructions are provided below.

# Depth indented listing of files
``` bash

```

# Getting Started
## 1. Deploy workflow
The repository can be downloaded with the additional dependencies handled by PlantChimera. The repository has a primary directory `PlantChimera` with several subfolders. 

```bash 
    git clone https://github.com/rebminso/PlantChimera.git
    cd PlantChimera
```

## 2. Installation
To simplify using the program, a Dockerfile is provided which will create a Docker image suitable for 
all the required packages. This README covers the use of the Docker version of the tool. 


### (A1.) Pulling the Docker Image from Docker Hub
``` bash 
    docker pull nmrtsahu/plantchimera:latest
```

### (A2.) Running the container interactively

#### a) Mounting the working directory of your local system
```bash
    docker run -it -v /path/to/local/directory:/app/PlantChimera/data nmrtsahu/plantchimera:latest /bin/bash
```
This will open a shell in the container, allowing you to view and interact with the files (e.g., scripts and the config.yaml file and the output) in the mounted directory on your local machine. Users can modify the parameters in the config.yaml file as needed to suit their requirements. Run the PlantChimera.sh  in a shell in the container. To exit the shell without stopping the container, press Ctrl + P followed by Ctrl + Q.

```bash
    PlantChimera.sh -h 
```
#### b) Or running without interacting the shell in the container  
```bash
    docker run -it -v /path/to/local/directory:/app/PlantChimera/data nmrtsahu/plantchimera:latest “ PlantChimera.sh -h” 
```
OR

### (B1.) Pull the Docker image into a Singularity Image File (SIF)

```bash
    singularity pull docker://nmrtsahu/plantchimera:latest
```
This will download the Docker image and convert it into a Singularity Image File (.sif) in the current directory.
```bash
singularity exec --bind /path/to/local/directory:/data plantchimera.sif /bin/bash
./PlantChimera.sh -h
```
## 3. Input Data Setup
a) The input FASTQ files and reference genome are available for download from [here](https://drive.google.com/drive/folders/1Sg8T8NXMD6t7qQO_SjlHSWl4_HRoEl6p?usp=drive_link ). 

b) Once you've downloaded the files, move them to the `PlantChimera/sample` directory. You can do this manually or using the following command:
```bash 
mv /path/to/downloaded/files/*.fastq.gz /path/to/PlantChimera/sample/
```

## 4. Configure workflow
The configuration of the workflow can be altered by modifying the `config/config.yaml`. The imporant directories are `Input`, `Scripts` and `Output`.  


## 5. Run PlantChimera


```bash
    Usage: ./PlantChimera.sh [options]

    Options:
    -r <reference_file>   Path to the reference genome file (.fasta or .fa) (required)
    -i <input_file>       Path to the forward read of paired end sequencing data (required)
    -I <input_file>       Path to the reverse read of paired end sequencing data (required)
    -g <gtf_file>         Path to the genome annotation file (.gtf)(required)
    -T <transcriptome_file> Path to the reference transcriptome file (.fasta or .fa)(required)
    -o <SAMPLE_OUTPUT>    Path to the output folder (required), Note: only enter the name of the sample eg SRR16989272 
    -t <threads>          Number of threads to use (default: 4)
    -s <species>          Species specific output folder (required), Note: folder name should be without space eg. arabidopsis_thaliana or ath
    -p <paralogue_gene>   Path to the paralogue gene file generated from get_paralogues.R Rscript
    -h                    Display this help message

    Description:
    PLantChimera: A novel fusion detection pipeline for plants.

```

## 6. Output

The output of PlantChimera named `PlantChimera_fusions.csv` is a tab delimited file with the following format. Ther are several other intermediatory files generated during the process which might be useful for user.

``` bash
Fusion_Name 	5'Gene ID	5'Breakpoint	3'Gene ID	3'Breakpoint	LeftBreakEntropy	RightBreakEntropy	Sites	t1_region	t2_region	SplicePattern	FusionJunctionSequence	%Homology	SRC	SpliceSiteClass
AT2G34420_AT2G34430 	AT2G34420	chr2:14523214:-	AT2G34430	chr2:14525272:+	1.95	1.98	INTRACHROMOSOMAL	exon_mid	exon_mid	CA_GT	GCTAGAAGTTATCCA-CCACGCTCAGAGCAT	26.7	2	NonCanonicalPattern
AT3G09162_AT3G09160 	AT3G09162	chr3:2807846:-	AT3G09160	chr3:2806948:-	1.94	1.99	INTRACHROMOSOMAL	exon_mid	exon_mid	GT_TC	TGGAATTTGATTCAG-TTAAGGGTTATCGCC	33.3	8	NonCanonicalPattern
AT3G28620_AT3G28630	    AT3G28620	chr3:10728441:+	AT3G28630	chr3:10729044:+	1.99	1.94	INTRACHROMOSOMAL	exon_mid	exon_mid	TT_TA	TCTGTCGCAAACCTG-CCTGAAATTGTTTCA	20.0	4	NonCanonicalPattern
AT4G06534_AT4G06536 	AT4G06534	chr4:3357768:+	AT4G06536	chr4:3360752:+	1.94	1.84	INTRACHROMOSOMAL	exon_mid	exon_mid	TA_AG	CATCTATCTCGATGG-GAATGAAGCTGGTTT	26.7	14NonCanonicalPattern
AT5G28390_AT3G61740 	AT5G28390	chr5:10342026:-	AT3G61740	chr3:22855685:-	1.93	1.92	INTERCHROMOSOMAL	exon_ter	exon_ter	TA_AC	CGTTACAATCCTTAT-AAGTAAGTACATGAG	33.3	2	NonCanonicalPattern
AT5G28390_AT3G61740 	AT5G28390	chr5:10342058:-	AT3G61740	chr3:22855688:-	1.89	1.95	INTERCHROMOSOMAL	exon_mid	exon_mid	AT_TG	AACAAAGTAGCGACA-GACAAGTAAGTACAT	40.0	2	NonCanonicalPattern
AT5G57140_AT5G57150 	AT5G57140	chr5:23151392:+	AT5G57150	chr5:23152184:+	1.91	1.92	INTRACHROMOSOMAL	exon_mid	exon_mid	GT_CA	GACAAAACATATAAG-GTAAGATACAACGGC	26.7	15NonCanonicalPattern
ATCG00065_ATCG00905 	ATCG00065	chrPt:69611:-	ATCG00905	chrPt:98793:-	1.95	1.9	INTRACHROMOSOMAL	exon_ter	exon_ter	TG_CC	ATGTACTCGGGTGTA-ACTATCACCCCCAAA	26.7	9	NonCanonicalPattern
ATCG00065_ATCG01230 	ATCG00065	chrPt:69611:-	ATCG01230	chrPt:139856:+	1.95	1.9	INTRACHROMOSOMAL	exon_ter	exon_ter	TG_CA	ATGTACTCGGGTGTA-CTATCACCCCCAAAA	26.7	6	NonCanonicalPattern
```
Description of each column is provided below:

| **Columns**     | **Description**      | 
| ------------- | ------------  |
| Fusion_Name    |  Name for gene fusion event.|
| 5' Gene ID    | Gene ID of the upstream gene (5' end) involved in the fusion.|
|5'Breakpoint| Genomic location of the breakpoint for the 5' gene. |
|3' Gene ID |Gene ID of the downstream gene (3' end) involved in the fusion. |
|3'Breakpoint |Genomic location of the breakpoint for the 3' gene.|
|LeftBreakEntropy |Entropy score at the 5' breakpoint.|
|RightBreakEntropy|Entropy score at the 3' breakpoint.|
|Sites|Type of chromosomal event (e.g., intra-/inter-chromosomal). |
|FusionJunctionSequence |equence around the fusion junction.|
|SplicePattern| Splice site pattern at the fusion junction.|
|SRC| Supportive read count for the fusion. |
|SpliceSiteClass| Classification of the splice site (e.g., canonical, non-canonical).|



#

## 7. Pratical examples of using PlantChimera

A small input dataset has been provided in `input/sample/` that can be leveraged using modest computational resources, so it should run on any hardware with atlast 4GB ram. 


### Download and configure PlantChimera

``` bash
git clone https://github.com/rebminso/PlantChimera.git
cd PlantChimera
unzip input/sample/*    
docker pull nmrtsahu/plantchimera:latest 
sudo docker run -it -v $PWD:/data nmrtsahu/plantchimera:latest /bin/bash
```
### Run the following command to execute PlantChimera
``` bash
./PlantChimera.sh -r input/arab/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -i input/arab/dump_1.fastq  -I input/arab/dump_2.fastq -g input/arab/Arabidopsis_thaliana.TAIR10.56.gtf -T input/arab/Arabidopsis_thaliana.TAIR10.cdna.all.fa -o sample -t 8 -s ath -p paralogue_gene.txt
```
The output will be genrated in the output folder 


## 8. How to run get_paralogues.R Rscript to get paralogue gene:

To retrieve paralogue genes of Ensembl plants, you can find the dataset names in the provided `ensembl_plants_dataset.txt` file. 

### Instructions
1. Open the `ensembl_plants_dataset.txt` file to identify the appropriate Ensembl dataset name for your plant of interest.
2. Replace `athaliana_eg_gene` in `get_paraologue.R` script with the selected dataset name.
3. Save your changes and run the script to retrieve the paralogue genes.

Following these steps will enable you to access paralogue genes for various plant species using Ensembl datasets.

1. **Run the get_paralogue.R Script**:
   Use the following bash command to run the `get_paraologue.R` script:

   ```bash
   Rscript get_paraologue.R genome_annotation.gtf transcriptome.fasta paralogue_gene.txt
   ```

## 9. Authors

- [Garima Kalakoti](mailto:kalakoti09@gmail.com)
- [Namrata Sahu](mailto.sahunamrata2098@gmail.com)
- [Manish Datt](mailto:manishdatt@gmail.com)

If you have any questions, bug reports, or suggestions, please contact:

**Dr. Shailesh Kumar**  
Staff Scientist, Bioinformatics Laboratory #202  
National Institute of Plant Genome Research (NIPGR), New Delhi  
Email: [shailesh@nipgr.ac.in](mailto:shailesh@nipgr.ac.in)

