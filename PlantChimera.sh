#!/bin/bash

REFERENCE_FILE=""
INPUT1_FILE=""
INPUT2_FILE=""
ANNOT_FILE=""
TRANS_FILE=""
SAMPLE_OUTPUT=""
PARALOGUE_FILE=""
THREADS=4  # Default value for threads
SPECIES=""

# Function to display usage
usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Options:"
    echo "  -r <reference_file>   Path to the reference genome file (.fasta or .fa) (required)"
    echo "  -i <input_file>       Path to the forward read of paired end sequencing data (required)"
    echo "  -I <input_file>       Path to the reve rse read of paired end sequencing data (required)"
    echo "  -g <gtf_file>         Path to the genome annotation file (.gtf)(required)"
    echo "  -T <transcriptome_file> Path to the reference transcriptome file (.fasta or .fa)(required)"
    echo "  -o <SAMPLE_OUTPUT>    Path to the output folder (required), Note: only enter the name of the sample eg SRR16989272 "
    echo "  -t <threads>          Number of threads to use (default: 4)"
    echo "  -s <species>          Species specific output folder (required), Note: folder name should be without space eg. arabidopsis_thaliana or ath"
    echo "  -p <paralogue_gene>   Path to the paralogue gene file"
    echo "  -h                    Display this help message"
    echo ""
    echo "Description:"
    echo " PLantChimera: Fusion detection novel pipeline for plants."
    echo ""
    exit 1
}

while getopts ":r:i:I:g:T:o:t:s:p:h" opt; do
    case ${opt} in
        r )
            REFERENCE_FILE=$OPTARG
            ;;
        i )
            INPUT1_FILE=$OPTARG
            ;;
        I )
            INPUT2_FILE=$OPTARG
            ;;
        g )
            ANNOT_FILE=$OPTARG
            ;;
        T )
            TRANS_FILE=$OPTARG
            ;;
        o )
            SAMPLE_OUTPUT=$OPTARG
            ;;
        t )
            THREADS=$OPTARG
            ;;
        s )
            SPECIES=$OPTARG
            ;;
        p)
            PARALOGUE_FILE=$OPTARG
            ;;
        h )
            usage
            ;;
        \? )
            usage
            ;;
        : )
            echo "Invalid option: -$OPTARG requires an argument" >&2
            usage
            ;;
    esac
done
shift $((OPTIND -1))

# Check if all required options are provided
if [ -z "$REFERENCE_FILE" ] || [ -z "$INPUT1_FILE" ] || [ -z "$INPUT2_FILE" ] || [ -z "$ANNOT_FILE" ] || [ -z "$TRANS_FILE" ] || [ -z "$SAMPLE_OUTPUT" ] || [ -z "$SPECIES" ]; then
    usage
fi



# Define directories dynamically based on species name
OUTPUT="$PWD/output/$SPECIES"           # For universal output like index, paralogous gene, transcript info
OUTPUT_DIR="$PWD/output/$SPECIES/$SAMPLE_OUTPUT"  # For sample-specific folder
SCRIPT_DIR="$(dirname "$(realpath "$0")")/scripts"

# Create directories if they do not exist
mkdir -p "$OUTPUT"
mkdir -p "$OUTPUT_DIR"

echo "Sample-specific output directory: $OUTPUT_DIR"

# read config file 
config_dir="$(dirname "$(realpath "$0")")/"
config_file="$config_dir/config.yaml"
# Accessing top-level value
with_paralogous_genes=$(yq eval '.with_paralogous_genes' "$config_file")
# Accessing values under the parameters key
wordsize=$(yq eval '.parameters.blast.wordsize' "$config_file")
perc_identity=$(yq eval '.parameters.blast.perc_identity' "$config_file")
num_threads=$(yq eval '.parameters.blast.num_threads' "$config_file")
max_parallel=$(yq eval '.parameters.blast.max_parallel' "$config_file")
AnchorLength=$(yq eval '.parameters.fusion_filter.AnchorLength' "$config_file")
ReadLength=$(yq eval '.parameters.fusion_filter.ReadLength' "$config_file")
JunctionSeq=$(yq eval '.parameters.fusion_filter.JunctionSeq' "$config_file")
overlap_query=$(yq eval '.parameters.fusion_fil_devmode.overlap_query' "$config_file")
gap=$(yq eval '.parameters.fusion_fil_devmode.gap' "$config_file")
subject_query=$(yq eval '.parameters.fusion_fil_devmode.subject_query' "$config_file")
gap_subject=$(yq eval '.parameters.fusion_fil_devmode.gap_subject' "$config_file")
JunctionDist=$(yq eval '.parameters.fusion_fil_devmode.JunctionDist' "$config_file")
ShanEnt_Seq=$(yq eval '.parameters.fusion_fil_devmode.ShanEnt_Seq' "$config_file")
BptDinucleotide=$(yq eval '.parameters.fusion_fil_devmode.BptDinucleotide' "$config_file")
TopHits=$(yq eval '.parameters.fusion_fil_devmode.TopHits' "$config_file")
Promiscioushits=$(yq eval '.parameters.fusion_fil_devmode.Promiscioushits' "$config_file")
split=$(yq eval '.parameters.fusion_fil_devmode.split' "$config_file")
span=$(yq eval '.parameters.fusion_fil_devmode.span' "$config_file")




#-----STEP 1 - INDEXING----- 
# Create index directory within the output folder
INDEX_DIR="$OUTPUT/index"
mkdir -p "$INDEX_DIR"
# Define the index prefix inside the index directory
INDEX_PREFIX="$INDEX_DIR/index"

# Extract the base name of the input file
#$BASE_NAME=$(basename "$INPUT1_FILE" | sed -E 's/(_[12]\.fastq|\.[12]\.fastq)//')

# Check if index files exist
INDEX_FILES=("${INDEX_PREFIX}.amb" "${INDEX_PREFIX}.ann" "${INDEX_PREFIX}.bwt" "${INDEX_PREFIX}.pac" "${INDEX_PREFIX}.sa")

INDEX_EXISTS=true
for file in "${INDEX_FILES[@]}"; do
    if [ ! -f "$file" ]; then
        INDEX_EXISTS=false
        break
    fi
done
echo ""
echo "[ Step 1/10 ] : Indexing"
start_time=$(date +%s)
if [ "$INDEX_EXISTS" = true ]; then
    echo "Index files already exist. Skipping indexing."
else
    echo "Index files not found. Running BWA index..."
    bwa index -p "$INDEX_PREFIX" -a bwtsw "$TRANS_FILE"
    if [ $? -ne 0 ]; then
        echo "Error in indexing the transcriptome."
        exit 1
    fi
    echo "Indexing Process is completed."
fi
end_time=$(date +%s)
echo "Indexing took $((end_time - start_time)) seconds."
echo ""

#-----STEP 2 - ALIGNMENT AND PROCESSING-----
echo "[ Step  2/10 ] : Alignment and processing"

bam_file="$OUTPUT_DIR/${SAMPLE_OUTPUT}_sorted.bam"

if [ -f "$bam_file" ]; then
    echo "BAM file already exists. Skipping alignment and sorting."
    echo ""
else
    start_time=$(date +%s)

    if [ "$ReadLength" -gt 70 ]; then
        echo "Read length > 70. Running BWA MEM..."
        bwa mem -t "$THREADS" "$INDEX_PREFIX" "$INPUT1_FILE" "$INPUT2_FILE" > "$OUTPUT_DIR/${SAMPLE_OUTPUT}.sam"
    else
        echo "Read length <= 70. Running BWA ALN..."
        bwa aln -t "$THREADS" "$INDEX_PREFIX" "$INPUT1_FILE" > "$OUTPUT_DIR/${SAMPLE_OUTPUT}_1.sai"
        bwa aln -t "$THREADS" "$INDEX_PREFIX" "$INPUT2_FILE" > "$OUTPUT_DIR/${SAMPLE_OUTPUT}_2.sai"
        echo "Running BWA SAMPE..."
        bwa sampe "$INDEX_PREFIX" "$OUTPUT_DIR/${SAMPLE_OUTPUT}_1.sai" "$OUTPUT_DIR/${SAMPLE_OUTPUT}_2.sai" "$INPUT1_FILE" "$INPUT2_FILE" > "$OUTPUT_DIR/${SAMPLE_OUTPUT}.sam"
    fi

    if [ $? -ne 0 ]; then
        echo "Error running alignment."
        exit 1
    fi

    echo "Alignment is completed."
    end_time=$(date +%s)
    echo "Alignment took $((end_time - start_time)) seconds."
    echo ""

    #SAM to BAM processing
    start_time=$(date +%s)
    echo "Processing SAM file with Samtools..."
    samtools view -bS "$OUTPUT_DIR/${SAMPLE_OUTPUT}.sam" -o "$OUTPUT_DIR/temp.bam"

    if [ $? -ne 0 ]; then
        echo "Error converting SAM to BAM."
        exit 1
    fi

    #bam sorting 
    samtools sort -n "$OUTPUT_DIR/temp.bam" -o "$bam_file"
    if [ $? -ne 0 ]; then
        echo "Error sorting BAM file."
        exit 1
    fi

    rm "$OUTPUT_DIR/temp.bam"
    end_time=$(date +%s)
    echo "SAM to BAM processing took $((end_time - start_time)) seconds."
    echo ""
fi

#-----STEP 3 - BAM TO BED CONVERSION-----
echo "[ Step 3/10 ] : BAM to BED Conversion"

start_time=$(date +%s)
bed_file="$OUTPUT_DIR/${SAMPLE_OUTPUT}_PE.bedpe"

if [ -f "$bed_file" ]; then
    echo "$bed_file file alredy exists. Skipping the BAM to BED conversion step"
    echo ""
else
    echo "BAM file to BED conversion is running..."
    bamToBed -bedpe -i "$OUTPUT_DIR/${SAMPLE_OUTPUT}_sorted.bam" > "$OUTPUT_DIR/${SAMPLE_OUTPUT}_PE.bedpe"
    echo "BAM file to BED conversion is completed."
    end_time=$(date +%s)
    echo "BAM to BED conversion took $((end_time - start_time)) seconds."
    echo ""
    if [ $? -ne 0 ]; then
        echo "Error running BAM to BED conversion."
        exit 1
    fi    
fi

#-----STEP 4 - SA_FILE generation step-----
echo "[ Step 4/10 ] : Extracting discordant reads "
# Timing the python3 script discordant_extracter.py
sa_file="$OUTPUT_DIR/${SAMPLE_OUTPUT}_SAfile.txt"

if [ -f "$sa_file" ]; then
    echo "$OUTPUT_DIR/${SAMPLE_OUTPUT}_SAfile.txt file already exists. Skipping SA_file generation step."
    echo ""

else
    start_time=$(date +%s)
    python3 $SCRIPT_DIR/discordant_extracter.py "$OUTPUT_DIR/${SAMPLE_OUTPUT}_PE.bedpe" "$OUTPUT_DIR/${SAMPLE_OUTPUT}_SAfile.txt"
    if [ $? -ne 0 ]; then
        echo "Error running python3 file discordant_extracter.py."
        exit 1
    fi
    #echo "SA_file.txt generated in output directory $OUTPUT_DIR."
    end_time=$(date +%s)
    echo "discordant_extracter.py took $((end_time - start_time)) seconds."
    echo ""
fi

#----- STEP 5 - removal of paralogues gene fusion pair -----
echo "[ Step 5/10 ] : removal of paralogues gene fusion pair"
start_time=$(date +%s)
paralogue_remover_out="$OUTPUT_DIR/${SAMPLE_OUTPUT}_newSAfile.txt"

if [ ! -f "$paralogue_remover_out" ]; then
    echo "Removing duplicates gene pairs and paralogues gene pairs from SAfile..."
    python3 "$SCRIPT_DIR/paralogue_remover.py" "$OUTPUT_DIR/${SAMPLE_OUTPUT}_SAfile.txt" "$ANNOT_FILE" "$OUTPUT/gtf_transcript_info.txt" "$OUTPUT_DIR/${SAMPLE_OUTPUT}_newSAfile.txt" "$PARALOGUE_FILE"
    end_time=$(date +%s)
    echo "paralogue_remover.py took $((end_time - start_time)) seconds."
    echo ""

    if [ $? -ne 0 ]; then
        echo "Error running python3 file paralogue_remover.py."
        exit 1
    fi

else
    echo "$paralogue_remover_out already exists. Skipping the removal of homologous genes and same gene fusions."
    echo ""
fi

#----- STEP 6 - Blastn running -----
mkdir -p "$OUTPUT_DIR/temp"
echo "[ Step 6/10 ] : Blastn running ..."

blast_out="$OUTPUT_DIR/${SAMPLE_OUTPUT}_SAfile_blast.txt"
start_time=$(date +%s) 

if [ ! -f "$blast_out" ]; then
    echo "Running blast_parallel.py..."

    python3 $SCRIPT_DIR/blastn_parallel_optimized.py --input "$INPUT1_FILE,$INPUT2_FILE" --transcript "$TRANS_FILE" --safile "$OUTPUT_DIR/${SAMPLE_OUTPUT}_newSAfile.txt" --temp_path "$OUTPUT_DIR/temp" --output "$blast_out" --num_threads $num_threads --max_parallel $max_parallel --wordsize $wordsize --perIdentity $perc_identity
    end_time=$(date +%s)
    echo "blast_parallel.py took $((end_time - start_time)) seconds."
    echo ""
    rm -r "$OUTPUT_DIR/temp"

    if [ $? -ne 0 ]; then
        echo "Error running blast_parallel.py."
        exit 1
    fi
else 
    echo "Blast output is already in $blast_out. Skipping running the Blastn step."
    echo ""
fi

#----- STEP 7 - read_processor.py.py -----
echo "[ Step 7/10 ] : Running read_processor.py"

read_processor_out="$OUTPUT_DIR/${SAMPLE_OUTPUT}_df_blast_output.txt"

if [ ! -f $read_processor_out ]; then 
    
    start_time=$(date +%s)

    python3 $SCRIPT_DIR/read_processor.py "$OUTPUT_DIR/${SAMPLE_OUTPUT}_SAfile_blast.txt" "$OUTPUT_DIR/${SAMPLE_OUTPUT}_newSAfile.txt" "$OUTPUT_DIR/${SAMPLE_OUTPUT}_fusiondf.csv" "$OUTPUT_DIR/${SAMPLE_OUTPUT}_df_blast_output.txt"

    if [ $? -ne 0 ]; then
        echo "Error running read_processor.py."
        exit 1
    fi

    end_time=$(date +%s)

    echo "read_processor.py took $((end_time - start_time)) seconds"
    echo ""
else
    echo "$read_processor_out already exist. Skipping read_processor.py running step"
    echo ""
fi

#----- STEP 8 - chimera_breakpoint_detecter.py -----
start_time=$(date +%s)

bp_detector_out=$OUTPUT_DIR/${SAMPLE_OUTPUT}_allsplit.txt
echo "[ Step 8/10 ] : chimera_breakpoint_detecter.py"

if [ ! -f $bp_detector_out ]; then

    python $SCRIPT_DIR/chimera_breakpoint_detecter.py "$OUTPUT_DIR/${SAMPLE_OUTPUT}_df_blast_output.txt" $overlap_query $gap "$OUTPUT_DIR/${SAMPLE_OUTPUT}_allsplit.txt" $AnchorLength $ReadLength $subject_query $gap_subject  

    if [ $? -ne 0 ]; then
        echo "Error running chimera_breakpoint_detecter.py"
        exit 1
    fi


    end_time=$(date +%s)
    echo "chimera_breakpoint_detecter.py took $((end_time - start_time)) seconds."
    echo ""
else 
    echo "$bp_detector_out already exists. Skipping the step "
    echo ""
fi

#--------STEP 9 : chimera_iden_modified.py----------
echo "[ Step 9/10 ] : chimera_identifier.py"

chimera_iden=$OUTPUT_DIR/${SAMPLE_OUTPUT}_chimera_identified.csv

if [ ! -f $chimera_iden ]; then 

    start_time=$(date +%s)  
    python $SCRIPT_DIR/chimera_iden_modified.py $ANNOT_FILE $OUTPUT/transcript_info.csv $OUTPUT/exons_position_adjusted.csv $OUTPUT_DIR/${SAMPLE_OUTPUT}_allsplit.txt $OUTPUT_DIR/${SAMPLE_OUTPUT}_chimera_identified.csv $OUTPUT_DIR/${SAMPLE_OUTPUT}_fusiondf.csv
    if [ $? -ne 0 ]; then
        echo "Error running chimera_identifier.py."
        exit 1nan
    fi
    end_time=$(date +%s)
    echo "chimera_identifier.py took $((end_time - start_time)) seconds."
    echo ""
else
    echo ""${SAMPLE_OUTPUT}_chimera_identified.csv" already exists in $OUTPUT_DIR directory. Skipping the chimera identification step"
    echo ""
fi

#--------STEP 10 : chimera_filter.py----------
echo "[ Step 10/10 ] : chimera_filter.py"

start_time=$(date +%s)
python $SCRIPT_DIR/chimera_filter.py $OUTPUT_DIR/${SAMPLE_OUTPUT}_chimera_identified.csv $JunctionDist $ShanEnt_Seq $REFERENCE_FILE $JunctionSeq $BptDinucleotide $TopHits $Promiscioushits $split $span $OUTPUT_DIR/${SAMPLE_OUTPUT}_PlantChimera_fusions.csv
if [ $? -ne 0 ]; then
    echo "Error running chimera_filter.py."
    exit 1
fi

#remove intermediate files 
rm -r $OUTPUT_DIR/${SAMPLE_OUTPUT}_fusiondf.csv
rm -r $OUTPUT_DIR/${SAMPLE_OUTPUT}_chimera_identified.csv
rm -r $OUTPUT_DIR/${SAMPLE_OUTPUT}_df_blast_output.txt
rm -r $OUTPUT_DIR/${SAMPLE_OUTPUT}_allsplit.txt
rm -r $OUTPUT_DIR/${SAMPLE_OUTPUT}_SAfile_blast.txt
rm -r $OUTPUT_DIR/${SAMPLE_OUTPUT}_newSAfile.txt
rm -r $OUTPUT_DIR/${SAMPLE_OUTPUT}_SAfile.txt
rm -r $OUTPUT_DIR/${SAMPLE_OUTPUT}_PE.bedpe
#rm -r $OUTPUT_DIR/${SAMPLE_OUTPUT}.sam
#rm -r $OUTPUT_DIR/*.fasta $OUTPUT_DIR/*.bed 

echo "chimera_filter.py took $((end_time - start_time)) seconds."
echo ""
echo "Script completed successfully. Output file saved in "$OUTPUT_DIR/${SAMPLE_OUTPUT}_PlantChimera_fusions.csv" directory"

