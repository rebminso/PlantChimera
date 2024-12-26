```
 ____  _             _    ____ _     _                          
|  _ \| | __ _ _ __ | |_ / ___| |__ (_)_ __ ___   ___ _ __ __ _ 
| |_) | |/ _` | '_ \| __| |   | '_ \| | '_ ` _ \ / _ \ '__/ _` |
|  __/| | (_| | | | | |_| |___| | | | | | | | | |  __/ | | (_| |
|_|   |_|\__,_|_| |_|\__|\____|_| |_|_|_| |_| |_|\___|_|  \__,_|

```
PlantChimera aims to provide a streamlined, efficient, and accurate method for identifying fusion transcripts in plant transcriptomic datasets. It employs a highly configurable pipeline with stringent alignment, filtering to detect fusion events from RNA-Seq data with high confidence. By integrating multiple layers of filtering, from sequence overlap and gene orientation to entropy and junction distance, the pipeline reduces the likelihood of false positives while retaining biologically significant fusion candidates. The source code can be found at [https://github.com/rebminso/PlantChimera.git](https://github.com/rebminso/PlantChimera.git) and is open source and compatible with multiple platforms.

![image](https://github.com/user-attachments/assets/0112bb73-2f02-4a3b-9afb-8048204714bc)






# Depth indented listing of files
``` bash
PlantChimera/
    ├── config.yaml
    ├── Dockerfile
    ├── PlantChimera.sh
    ├── README.md
    ├── Paralogues
    │   ├── ensembl_plants_dataset.txt
    │   ├── get_paralogue.R
    │   └── ath_paralogue_gene.txt
    └── scripts
        ├── blastn_parallel_optimized.py
        ├── chimera_breakpoint_detecter.py
        ├── chimera_filter.py
        ├── chimera_identifier.py
        ├── discordant_extracter.py
        ├── paralogue_remover.py
        └── read_processor.py
```



# Getting Started
## 1. Deploy workflow
The repository can be downloaded with the additional dependencies by following commands. The repository has a primary directory `PlantChimera` with several subfolders. 

```bash 
    git clone https://github.com/rebminso/PlantChimera.git
    cd PlantChimera
```

## 2. Installation
To simplify using the program, a Dockerfile is provided which will create a Docker image suitable for 
all the required packages. This README covers the usage of the Docker version of the tool. 


### (A1.) Pulling the Docker Image from Docker Hub
``` bash 
    docker pull nmrtsahu/plantchimera:v1.0
```

### (A2.) Running the container interactively

#### a) Mounting the working directory of your local system
```bash
    docker run -it -v /path/to/local/directory:/app/PlantChimera/data nmrtsahu/plantchimera:v1.0 /bin/bash
```
This will open a shell in the container, allowing you to view and interact with the files (e.g., scripts and the config.yaml file and the output) in the mounted directory on your local machine. Users can modify the parameters in the config.yaml file as needed to suit their requirements. Run the PlantChimera.sh  in a shell in the container. To exit the shell without stopping the container, press Ctrl + P followed by Ctrl + Q.

```bash
    PlantChimera.sh -h 
```
#### b) Or running without interacting the shell in the container  
```bash
    docker run -it -v /path/to/local/directory:/app/PlantChimera/data nmrtsahu/plantchimera:v1.0 “ PlantChimera.sh -h” 
```
OR

### (B1.) Pull the Docker image into a Singularity Image File (SIF)

```bash
    singularity pull docker://nmrtsahu/plantchimera:v1.0
```
This will download the Docker image and convert it into a Singularity Image File (.sif) in the current directory.
```bash
singularity exec --bind /path/to/local/directory:/data plantchimera.sif /bin/bash
./PlantChimera.sh -h
```
## 3. Input Data Setup
a) The input FASTQ files and reference genome/gtf/transcriptome are available for download from [here](https://drive.google.com/drive/folders/1Sg8T8NXMD6t7qQO_SjlHSWl4_HRoEl6p?usp=drive_link ). 

b) Once you've downloaded the files, move them to the `PlantChimera/sample` directory. You can do this manually or using the following command:
```bash 
mv /path/to/downloaded/files/*.fastq.gz /path/to/PlantChimera/sample/
```

## 4. Configure workflow
The configuration of the workflow can be altered by modifying the `config/config.yaml` if and only if required. The imporant directories are `sample`, `scripts`.  


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
Fusion_Name	    5'Gene ID	5'Breakpoint	5'Transcript_ID	    3'Gene ID	3'Transcript_ID    	3'Breakpoint	LeftBreakEntropy	RightBreakEntropy	FusionJunctionSequence	SRC	Sites	Breakpoint Pattern
AT1G11890::AT1G31910	AT1G11890	chr1:4010892:+	AT1G11890.1	AT1G31910	    AT1G31910.2     chr1:11458686:+	1.94	1.91	    ATCTCCTATCAGGTT_TCCGAATCCGATCTA	4.0	INTRACHROMOSOMAL	E/E
AT1G14790::AT1G76930	AT1G14790	chr1:5096408:-	AT1G14790.1	AT1G76930	    AT1G76930.2	    chr1:28895934:-	1.96	1.8         TGGATGTTTGCACCA_TACCACTCTCCACCA	4.0	INTRACHROMOSOMAL	M/M
AT1G19720::AT1G19715	AT1G19720	chr1:6819716:-	AT1G19720.1	AT1G19715	    AT1G19715.1	    chr1:6819504:-	1.92	1.91	    GGAAGTTGACCAAAG_CAGTGGTCAGCTATC	3.0	INTRACHROMOSOMAL	M/M
AT1G21313::AT5G58470	AT1G21313	chr1:7453979:+	AT1G21313.1	AT5G58470	    AT5G58470.1	    chr5:23640209:-	1.82	1.76	    CTTTGGTGGTGGTGG_AGGTGGAGGTGACAG	4.0	INTERCHROMOSOMAL	M/M
AT1G24996::AT1G24880	AT1G24996	chr1:8800352:-	AT1G24996.1	AT1G24880	    AT1G24880.2	    chr1:8776144:+	1.99	1.95	    TATGGACCTACCTGA_ATCTGCTATCCGCAG	6.0	INTRACHROMOSOMAL	M/E
AT1G25097::AT1G24880	AT1G25097	chr1:8814730:-	AT1G25097.1	AT1G24880	    AT1G24880.2	    chr1:8776144:+	1.99	1.95	    TATGGACCTACCTGA_ATCTGCTATCCGCAG	4.0	INTRACHROMOSOMAL	M/E
AT1G28710::AT1G28765	AT1G28710	chr1:10089444:-	AT1G28710.3	AT1G28765	    AT1G28765.1	    chr1:10107645:-	1.92	1.74	    CATACCACTAGACCA_TTTTTATAATGCTGT	4.0	INTRACHROMOSOMAL	M/E
AT1G30230::AT2G42660	AT1G30230	chr1:10639510:+	AT1G30230.1	AT2G42660	    AT2G42660.1	    chr2:17767874:-	1.81	1.84	    TGTGTATGCTTATAG_AGGTAATCTTCATCG	3.0	INTERCHROMOSOMAL	M/E
```
Description of each column is provided below:

| **Columns**     | **Description**      | 
| ------------- | ------------  |
| Fusion_Name    |  Name for gene fusion event.|
|5' Gene ID    | Gene ID of the upstream gene (5' end) involved in the fusion.|
|5'Breakpoint| Genomic location of the breakpoint for the 5' gene. |
|3' Gene ID |Gene ID of the downstream gene (3' end) involved in the fusion. |
|3'Breakpoint |Genomic location of the breakpoint for the 3' gene.|
|LeftBreakEntropy |Entropy score at the 5' breakpoint, indicating its randomness or complexity..|
|RightBreakEntropy|Entropy score at the 3' breakpoint, indicating its randomness or complexity..|
|Sites|Type of chromosomal event (e.g., intra-/inter-chromosomal). |
|FusionJunctionSequence |The nucleotide sequence at the fusion junction.|
|SRC| Supportive read count for the fusion. |
|Breakpoint Pattern| Location of the fusion transcript breakpoints at the exon boundaries [E] or within Exons[M].|



#

## 7. Pratical examples of using PlantChimera

A small input dataset can be downloaded from [here](https://drive.google.com/drive/folders/1Sg8T8NXMD6t7qQO_SjlHSWl4_HRoEl6p?usp=drive_link ) that can be leveraged using modest computational resources, so it should run on any hardware with atlast 4GB ram and make sure to transfer the downloaded files to the `sample/` directory 

### Download and configure PlantChimera

``` bash
git clone https://github.com/rebminso/PlantChimera.git
cd PlantChimera
mkdir sample
```

#### - Transfer the downloaded [files](https://drive.google.com/drive/folders/1Sg8T8NXMD6t7qQO_SjlHSWl4_HRoEl6p?usp=drive_link ) to the `/path/to/PlantChimera/sample/` directory. 

``` bash
unzip sample/*    
docker pull nmrtsahu/plantchimera:latest 
sudo docker run -it -v $PWD:/data nmrtsahu/plantchimera:v1.0 /bin/bash
```
### Run the following command to execute PlantChimera
``` bash
./PlantChimera.sh -r sample/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -i sample/sample_1.fastq  -I sample/sample_1.fastq -g sample/Arabidopsis_thaliana.TAIR10.56.gtf -T sample/Arabidopsis_thaliana.TAIR10.cdna.all.fa -o output -t 8 -s Ath -p Paralogues/ath_paralogue_gene.txt
```
The results will be generated in the `output/` folder 


## 8. How to run get_paralogues.R Rscript to get paralogue gene:

To retrieve paralogue genes of Ensembl plants, you can find the dataset names in the provided `Paralogues/ensembl_plants_dataset.txt` file. 

### Instructions
1. Download genome and transcriptome of the plant of interest in `Paralogues/` directory. 
2. Open the `Paralogues/ensembl_plants_dataset.txt` file to identify the appropriate Ensembl dataset name for your plant of interest.
3. Replace `athaliana_eg_gene` in `get_paraologue.R` script with the selected dataset name. 
4. Save your changes and run the script to retrieve the paralogues genes.

Following these steps will enable you to access paralogue genes for various plant species using Ensembl datasets.

1. **Run the get_paralogue.R Script**:
   Use the following bash command to run the `get_paraologue.R` script:

   ```bash
   Rscript paralogue/get_paraologue.R paralogue/genome.gtf paralogue/transcriptome.fasta paralogue_gene.txt
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

