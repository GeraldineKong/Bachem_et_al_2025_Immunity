#!/bin/sh

flag_manifest='manifest.tsv'
flag_outdir='kraken2_res'
flag_kraken_db=kraken2_pluspf/
flag_threads=30
inputerror='false'

echo ''

while getopts ':hi:o:d:t:' opt
do
    case ${opt} in
	h)
	    echo -e '\tPerforms Kraken2 taxonomic classification and abundance estimation with Bracken'
        echo -e ''
        echo -e '\tOptions:'
        echo -e ''
	    echo -e '\t -i Manifest File (headerless TSV file with sample name and absolute paths to demultiplexed forward and reverse reads) [default: manifest.tsv]'
		echo -e '\t -o Output Directory [default: kraken2_res/]'
        echo -e '\t -d Path to Kraken2 database [default: pluspf]'
		echo -e '\t -t Threads/CPUs To Use [default: 30]'
	    echo ''
	    exit 0
	    ;;
	i) flag_manifest=$OPTARG ;;
	o) flag_outdir=$OPTARG ;;
    d) flag_kraken_db=$OPTARG ;;
	t) flag_threads=$OPTARG ;;
	:) echo -e '\t Error: Use -h for full options list\n'
	   exit 1
    esac
done

# CONFIRM THAT MANIFEST FILE EXISTS
if [ ! -f $flag_manifest ]
then
	echo 'Manifest File Not Found'
    inputerror='true'
fi

# CONFIRM THAT FILES LISTED IN MANIFEST EXIST
while read sample read1 read2
do
	if [ ! -f $read1 ]
	then
		echo 'File '$read1' does not exist'
		inputerror='true'
	fi

	if [ ! -f $read2 ]
	then
		echo 'File '$read2' does not exist'
		inputerror='true'
	fi
done < $flag_manifest

# IF ANY OF THE ABOVE CONDITIONS ARE VIOLATED, PRINT MESSAGES AND EXIT
if [ $inputerror = 'true' ]
then
    echo ''
    exit 1
fi

# IF OUTPUT DIRECTORY ALREADY EXISTS
# WARN THAT CONTENTS WILL BE OVERWRITTEN
# WAIT 5 SECONDS TO GIVE USER TIME TO CANCEL
# CONTINUE
if [ ! -d $flag_outdir ]
then
	mkdir -p $flag_outdir
else
    echo ''
	echo 'Output Directory '$flag_outdir' Already Exists: Contents Will Be Overwritten'
    sleep 5
fi

# MAIN LOOP
while read sample read1 read2
do
  ## Kraken2 assignment
  kraken2 --paired $read1 $read2 \
    --threads $flag_threads \
    --gzip-compressed \
    --db $flag_kraken_db \
    --report $flag_outdir/${sample}_report.txt \
    --output $flag_outdir/${sample}_kraken2.out

  ## Bracken2 estimation
  /home/kongg1/tools/Bracken-2.7/bracken -d $flag_kraken_db \
                                         -i $flag_outdir/${sample}_report.txt \
                                         -o $flag_outdir/${sample}_bracken_species.txt \
                                         -r 100 \
                                         -l S \
                                         -t 10

  ## Kraken2 style report to MPA style
  python ~/tools/Bracken-2.7/KrakenTools-master/kreport2mpa.py -r $flag_outdir/${sample}_bracken_species.txt \
                                                               -o $flag_outdir/${sample}_bracken_mpa.txt \
                                                               --display-header

done < $flag_manifest

# Combining all reports into a table
cd $flag_outdir
python ~/tools/Bracken-2.7/KrakenTools-master/combine_mpa.py -i *_bracken_mpa.txt -o combined_bracken_mpa_report.txt
