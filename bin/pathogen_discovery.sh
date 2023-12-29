#!/bin/bash


###### Defining workflow parameters ######

OPTIONS=$(getopt -o a:b:c:d:e:f:g:h: --long sample_num:,read_type:,input_dir:,script_dir:,submit_type:,kraken_db:,read_1:,output_dir: -- "$@")

if [ $? -ne 0 ]; then
    echo "Error: Invalid option"
    exit 1
fi

eval set -- "$OPTIONS"

sample_num=""
read_type=""
input_dir=""
script_dir=""
submit_type=""
kraken_db=""
read_1=""
output_dir=""

while true; do
    case "$1" in
        --sample_num)
            sample_num="$2"
            shift 2
            ;;
        --read_type)
            read_type="$2"
            shift 2
            ;;
        --input_dir)
            input_dir="$2"
            shift 2
            ;;
        --script_dir)
            script_dir="$2"
            shift 2
            ;;
        --submit_type)
            submit_type="$2"
            shift 2
            ;;
        --kraken_db)
            kraken_db="$2"
            shift 2
            ;;
        --read_1)
            read_1="$2"
            shift 2
            ;;
        --output_dir)
            output_dir="$2"
            shift 2
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Internal error!"
            exit 1
            ;;
    esac
done

# set path to kraken db if specified
if [[ $read_type == "illumina-paired" ]]; then

	if [[ $read_1 == *"_r1"* ]]; then

		read_2="${read_1/_r1/_r2}"

	else

		read_2="${read_1/_R1/_R2}"

	fi

	read_1=$input_dir/$read_1
	read_2=$input_dir/$read_2

fi

if [[ $submit_type == "qsub" ]]; then

	blast_threads=""

else

	blast_threads="-num_threads 8"

fi


###### Defining Input and Reference directories ######

PROJECT_DIR=$script_dir/..
REF_DIR=$PROJECT_DIR/data/reference
ADAPTER_REF=/scicomp/reference/adapters/sequencing-adapters.fasta
BLASTDB=/scicomp/reference/ncbi-blast-taxdb


###### Defining output directories ######

if [[ -z $output_dir ]]; then

	OUTPUT_DIR=$PROJECT_DIR/results/$read_type

else

	# setting the output dir to user's input and removing a trailing forward slash if it was included
	OUTPUT_DIR="${output_dir%/}"

fi

BLAST_OUT=$OUTPUT_DIR/blast
CHOPPER_OUT=$OUTPUT_DIR/chopper
FASTP_OUT=$OUTPUT_DIR/fastp
FASTQC_OUT=$OUTPUT_DIR/fastqc
FLYE_OUT=$OUTPUT_DIR/flye
KRAKEN2_OUT=$OUTPUT_DIR/kraken2
REFORMAT_OUT=$OUTPUT_DIR/reformat
SPADES_OUT=$OUTPUT_DIR/spades
QUAST_OUT=$OUTPUT_DIR/quast

if [[ ! -d $OUTPUT_DIR ]]; then

	mkdir $OUTPUT_DIR

fi

output_dirs="blast chopper fastp fastqc flye kraken2 reformat spades quast"

# dynamically generating output dirs based on tools used
for i in $output_dirs; do

	if [[ ! -d "$OUTPUT_DIR/$i" ]]; then

		mkdir $OUTPUT_DIR/$i
	
	fi

	if [[ ! -d "$OUTPUT_DIR/$i/$sample_num" ]]; then

		mkdir $OUTPUT_DIR/$i/$sample_num

	fi

done

# redefining output paths
BLAST_OUT=$BLAST_OUT/$sample_num
CHOPPER_OUT=$CHOPPER_OUT/$sample_num
FASTQC_OUT=$FASTQC_OUT/$sample_num
FASTP_OUT=$FASTP_OUT/$sample_num
FLYE_OUT=$FLYE_OUT/$sample_num
KRAKEN2_OUT=$KRAKEN2_OUT/$sample_num
SPADES_OUT=$SPADES_OUT/$sample_num
QUAST_OUT=$QUAST_OUT/$sample_num

###### Defining Singularity Containers ######

BBMAP_CONTAINER="/scicomp/reference/singularity/containers/bbmap/bbmap%3A39.01--h92535d8_1"
BLAST_CONTAINER="/scicomp/reference/singularity/containers/blast/blast%3A2.14.1--pl5321h6f7f691_0"
CHOPPER_CONTAINER="/scicomp/reference/singularity/containers/chopper/chopper%3A0.7.0--hdcf5f25_0"
FASTP_CONTAINER="/scicomp/reference/singularity/containers/fastp/fastp%3A0.23.4--hadf994f_2"
FASTQC_CONTAINER="/scicomp/reference/singularity/containers/fastqc/fastqc%3A0.12.1--hdfd78af_0"
FLYE_CONTAINER="/scicomp/reference/singularity/containers/flye/flye%3A2.9.3--py39hd65a603_0"
KRAKEN2_CONTAINER="/scicomp/reference-pure/singularity/containers/kraken2/kraken2%3A2.1.3--pl5321hdcf5f25_0"
QUAST_CONTAINER="/scicomp/reference-pure/singularity/containers/quast/quast%3A5.2.0--py39pl5321h4e691d4_3"
SPADES_CONTAINER="/scicomp/reference-pure/singularity/containers/spades/spades%3A3.15.5--h95f258a_1"


#------ START PATHOGEN DISCOVERY WORKFLOW ------#


#File names format: fastq_pass/barcode##/*.fastq.gz

echo -e "\nSample number is: " $sample_num

echo -e "\nTrimming Adapter Sequences\n"

if [[ $read_type == "nanopore" || $read_type == "pacbio" ]]; then

	# Remove amplicon primers from reads. Reads should already be demultiplexed and barcode removed.
	gunzip -c $input_dir/$sample_num/*.fastq.gz | singularity exec $CHOPPER_CONTAINER chopper --headcrop 29 --tailcrop 29 | gzip > $CHOPPER_OUT/"$sample_num"_filtered_reads.fastq.gz

	TRIM_OUT=$CHOPPER_OUT/"$sample_num"_filtered_reads.fastq.gz

	# Subsample cleaned and chopped reads:
	singularity exec $BBMAP_CONTAINER reformat.sh \
		in=$TRIM_OUT \
		out=$REFORMAT_OUT/"$sample_num"_filtered_reads_SS.fastq.gz \
		samplerate=1 \
		ow=t \
		ignorebadquality=t \
		qin=33 \

	fastqc_args="$REFORMAT_OUT/"$sample_num"_filtered_reads_SS.fastq.gz"
	kraken2_args="$KRAKEN2_OUT/unclass-$sample_num.fastq $REFORMAT_OUT/"$sample_num"_filtered_reads_SS.fastq.gz"
	kraken2_std_args="$REFORMAT_OUT/"$sample_num"_filtered_reads_SS.fastq.gz"

else

	singularity exec $FASTP_CONTAINER fastp \
		-i $read_1 \
		-I $read_2 \
		-o $FASTP_OUT/"$sample_num"-R1-trimmed.fastq.gz \
		-O $FASTP_OUT/"$sample_num"-R2-trimmed.fastq.gz \
		-h $FASTP_OUT/"$sample_num"-fastp.html \
		-j $FASTP_OUT/"$sample_num"-fastp.json

	read_1_trimmed=$FASTP_OUT/"$sample_num"-R1-trimmed.fastq.gz
	read_2_trimmed=$FASTP_OUT/"$sample_num"-R2-trimmed.fastq.gz

	fastqc_args="$read_1_trimmed $read_2_trimmed"
	kraken2_args="$KRAKEN2_OUT/unclass-$sample_num-r#.fastq --paired $read_1_trimmed $read_2_trimmed"
	kraken2_std_args="--paired $read_1_trimmed $read_2_trimmed"

fi

echo -e "\nTrimming Complete\n"

# Perform QC on trimmed reads:
singularity exec $FASTQC_CONTAINER fastqc \
    --quiet \
    --threads 6 \
	--outdir $FASTQC_OUT \
    $fastqc_args

# Removing human reads with Kraken2:
singularity exec $KRAKEN2_CONTAINER kraken2 --threads 4 --db $kraken_db --unclassified-out $kraken2_args > $KRAKEN2_OUT/log.txt

# Setting up kraken2 standard DB params
if [[ $read_type == "illumina-paired" ]]; then

	kraken2_std_args="--paired $KRAKEN2_OUT/unclass-$sample_num-r_1.fastq $KRAKEN2_OUT/unclass-$sample_num-r_2.fastq"

else

	kraken2_std_args="$KRAKEN2_OUT/unclass-$sample_num.fastq"

fi

# Generating classification report for remaining reads with kraken2 standard DB
singularity exec $KRAKEN2_CONTAINER kraken2 --threads 4 --db /scicomp/reference/kraken_std8/ --report $KRAKEN2_OUT/classification_report.txt $kraken2_std_args > $KRAKEN2_OUT/log.txt

# De novo assembling unclassified reads (spades for illumina / flye for pacbio and nanopore):
echo -e "\nDe Novo assembly"

# predefining spades assembler as function

function spades_command {

	singularity exec $SPADES_CONTAINER spades.py \
		-k auto \
		-t 8 \
		-o $SPADES_OUT \
		$1

}


if [[ $read_type == "illumina-paired" ]]; then

	# running short reads through spades assembler

	spades "-1 $KRAKEN2_OUT/unclass-$sample_num-r_1.fastq -2 $KRAKEN2_OUT/unclass-$sample_num-r_2.fastq"
	
	ASSEMBLY_OUT=$SPADES_OUT/contigs.fasta

else
	
	if [[ $read_type == "nanopore" ]]; then

		flye_args="--nano-hq"

	else

		flye_args="--pacbio-hifi"

	fi

	# running long reads through flye assembler

	singularity exec $FLYE_CONTAINER flye \
		$flye_args \
		$KRAKEN2_OUT/unclass-$sample_num.fastq \
		--out-dir $FLYE_OUT

	if [[ ! -f "$FLYE_OUT/assembly.fasta" ]]; then

		# if long read assembler fails, brute force assembly with spades

		echo -e "\n!! WARNING !!: Long read assembly failed - read quality/length may be insufficient for a long read assembler. Attempting to assemble as short, unpaired reads with spades\nPlease evaluate assembly quality with QUAST results to determine the validity of the assembled contigs"

		spades_command "-s $KRAKEN2_OUT/unclass-$sample_num.fastq"

		ASSEMBLY_OUT=$SPADES_OUT/contigs.fasta

	else
	
		ASSEMBLY_OUT=$FLYE_OUT/assembly.fasta

	fi

fi

# Quast for assembly metrics:
echo -e "\nRunning quast for de novo assembly metrics:\n"

singularity exec $QUAST_CONTAINER quast.py --output-dir $QUAST_OUT $ASSEMBLY_OUT

# Identify contigs by nucleotide blast:
echo -e "\nRunning contigs through BLASTN\n"

function blast_command {

	blast_args="6 qseqid sacc pident length mismatch evalue bitscore stitle"

	singularity exec --env BLASTDB=/scicomp/reference/ncbi-blast-taxdb $BLAST_CONTAINER /usr/local/bin/blastn \
		-db "/scicomp/reference/ncbi-blast-databases/nt" \
		-query $1 \
		-outfmt "$blast_args" \
		-out $BLAST_OUT/$sample_num'_blast_hits.tsv' \
		-qcov_hsp_perc 70 \
		-perc_identity 70 \
		-evalue 0.0001 \
		-max_target_seqs 1 $blast_threads

	# formatting header string for blast output
	headers="${blast_args//6 /}"
	headers="${headers// /, }"

	echo -e "$headers\n$(cat $BLAST_OUT/$sample_num'_blast_hits.tsv')" > $BLAST_OUT/$sample_num'_blast_hits.tsv'

	if [[ $ASSEMBLY_OUT == *"spades"* && $read_type != "illumina-paired" ]]; then

		echo -e "# !!! WARNING !!!: BLAST Hits are a result of contigs generated from long reads forced through a spades short read assembly process. Please take assembly quality into consideration to validate BLAST hits\n\n$(cat $BLAST_OUT/$sample_num'_blast_hits.tsv')" > $BLAST_OUT/$sample_num'_blast_hits.tsv'

	fi

}

blast_command $ASSEMBLY_OUT

# if the long read contigs did not produce any BLAST hits, then we will run the unclass reads through spades (instead of flye) to generate a different set of contigs and potentially get a BLAST hit

if [[ $(grep -e "0 hits found" $BLAST_OUT/$sample_num'_blast_hits.tsv') && $ASSEMBLY_OUT == *"flye"* && $read_type != "illumina-paired" ]]; then

	echo -e "\nNo BLAST hits found with long read assembly - Rerunning unclassified reads through spades assembler"
	echo -e "\n!!! WARNING !!!: Attempting to assemble long reads as short, unpaired reads with spades\nPlease evaluate assembly quality with QUAST results to determine the validity of the assembled contigs"

	spades_command "-s $KRAKEN2_OUT/unclass-$sample_num.fastq"

	ASSEMBLY_OUT=$SPADES_OUT/contigs.fasta

	echo -e "\nRunning quast on forced spades assembly:\n"

	singularity exec $QUAST_CONTAINER quast.py --output-dir $QUAST_OUT $ASSEMBLY_OUT

	echo -e "\nRunning forced, spades contigs through BLASTN\n"

	blast_command $ASSEMBLY_OUT

fi

echo -e "\nScript finished\n"
