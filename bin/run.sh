#!/bin/bash

SCRIPTDIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )

# Defining script parameters
while getopts ":qt:hk:i:e:o:" opt; do
  case $opt in
    q) qsub=true ;;
    t) read_type="$OPTARG" ;;
    h) help=true ;;
    k) kraken_db="$OPTARG" ;;
    i) input_dir="$OPTARG" ;;
    e) email="$OPTARG" ;;
    o) output_dir="$OPTARG" ;;
    \?) echo "Invalid option -$OPTARG" >&2
        exit 1
    ;;
    :) echo "Option -$OPTARG requires an argument." >&2
       exit 1
    ;;
  esac
done


if [[ $help || -z $input_dir || -z $read_type ]]; then

  echo -e "\nUsage: $0 [OPTION]..."
  echo -e "Submits Pathogen Discovery script to cluster or runs the script iteratively over the sample directories on your local machine/host"
  echo -e "\nArgument list:\n"
  echo -e "-h\t\tHELP"
  echo -e "-t <STRING>\tread-type: illumina | pacbio | nanopore\t**REQUIRED**"
  echo -e "-i <STRING>\tpath to input directory\t**REQUIRED*"
  echo -e "-k <STRING>\tpath to kraken reference database - **REQUIRED*"
  echo -e "-q\t\tsubmit pathogen discovery script for each sample to the cluster with qsub \t(if -q is NOT specified, the script will run iteratively over each sample on a local, interactive node)"
  echo -e "-e <STRING>\temail for qsub notifications"
  echo -e "-o <STRING>\tpath to output directory - default is 'results' in this project directory"
  echo
  exit 1

elif [[ $read_type != "illumina-paired" && $read_type != "pacbio" && $read_type != "nanopore" ]]; then

  echo -e "ERROR: not a valid read type"
  echo -e "-t <STRING>\tread-type: illumina-paired | pacbio | nanopore\n"
  exit 1
  
fi

if [[ -z $kraken_db ]]; then

  kraken_db="/scicomp/reference/kraken-human/"

fi


function qsub_command {

  qsub \
    -N Pathogen_Discovery_$1 \
    -M $email \
    -m abe \
    -q highmem.q \
    -pe smp 8 \
    -l h_vmem=120G \
    -cwd \
    -o $SCRIPTDIR/../qsub_submissions \
    -e $SCRIPTDIR/../qsub_submissions/$1'_error.log' \
    $SCRIPTDIR/pathogen_discovery.sh --sample_num "$1" --read_type $read_type --input_dir "$2" --script_dir $SCRIPTDIR --kraken_db $kraken_db --submit_type "qsub" --read_1 "$3" --output_dir $output_dir

}

# setting up iteration logic for paired-reads
if [[ $read_type == "illumina-paired" ]]; then

  count=1
  
  for i in $(cat $SCRIPTDIR/../data/sample_names.txt); do

    for read in $(ls "$input_dir/$i"); do

      if [ $(( $count % 2 )) -eq 1 ]; then

        label="${read%%_r1*}"
        label="${label%%_R1*}"
        
        ((count++))

        if [[ $qsub ]]; then

          qsub_command $label $input_dir/$i $read

        else

          bash $SCRIPTDIR/pathogen_discovery.sh --sample_num $label --read_type $read_type --input_dir $input_dir/$i --script_dir $SCRIPTDIR --kraken_db "$kraken_db" --submit_type "interactive" --read_1 $read --output_dir $output_dir

        fi

      else

        ((count++))
        continue

      fi

    done

  done

else

  for i in $(cat $SCRIPTDIR/../data/sample_names.txt); do

    if [[ $qsub ]]; then

      qsub_command $i $input_dir
    
    else

      bash $SCRIPTDIR/pathogen_discovery.sh --sample_num $i --read_type $read_type --input_dir $input_dir --script_dir $SCRIPTDIR --kraken_db "$kraken_db" --submit_type "interactive" --output_dir $output_dir

    fi

  done

fi

