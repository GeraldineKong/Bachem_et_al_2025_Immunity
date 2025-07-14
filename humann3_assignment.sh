#!/bin/bash

flag_manifest='manifest.txt'
flag_merge_out='merged_fastq/'
flag_outdir='humann3_res/'
flag_metaphlan='MetaPhlAn-master/metaphlan'
flag_chocophlan_db='humann_database/chocophlan'
flag_threads=10
inputerror='false'

while getopts ':hm:m2:o:p:d:t:' opt; do
  case ${opt} in
    h)
      echo "Usage: $0 [-m manifest] [-m2 merged_dir] [-o output_dir] [-p metaphlan_path] [-d chocophlan_db] [-t threads]"
      echo ''
      echo "Options:"
      echo "  -m  Manifest file (3-column TSV with sample names, absolute path to read1, absolute path to read2) [default: manifest.txt]"
      echo "  -m2 Output directory for merged reads"
      echo "  -o  Output directory for HUMAnN3 results"
      echo "  -p  Path to MetaPhlAn executable"
      echo "  -d  Path to ChocoPhlAn nucleotide database"
      echo "  -t  Number of threads [default: 10]"
      exit 0
      ;;
    m) flag_manifest=$OPTARG ;;
    m2) flag_merge_out=$OPTARG ;;
    o) flag_outdir=$OPTARG ;;
    p) flag_metaphlan=$OPTARG ;;
    d) flag_chocophlan_db=$OPTARG ;;
    t) flag_threads=$OPTARG ;;
    :) echo "Error: Missing argument for -$OPTARG" >&2; exit 1 ;;
    \?) echo "Error: Invalid option -$OPTARG" >&2; exit 1 ;;
  esac
done

# Check manifest
if [ ! -f "$flag_manifest" ]; then
  echo "Error: Manifest file '$flag_manifest' not found."
  inputerror='true'
fi

# Check directories
mkdir -p "$flag_merge_out"
mkdir -p "$flag_outdir"

# Check sample files
while read sample read1 read2; do
  [ ! -f "$read1" ] && echo "Missing file: $read1" && inputerror='true'
  [ ! -f "$read2" ] && echo "Missing file: $read2" && inputerror='true'
done < "$flag_manifest"

if [ "$inputerror" = "true" ]; then
  echo "Aborting due to input errors."
  exit 1
fi

# Run HUMAnN3
while read sample read1 read2; do
  echo "Processing $sample..."

  # Merge reads
  merged="$flag_merge_out/${sample}_merged.fastq"
  zcat "$read1" \
       "$read2" > "$merged"

  # Run HUMAnN3
  humann --input "$merged" \
         --output "$flag_outdir" \
         --threads "$flag_threads" \
         --metaphlan "$flag_metaphlan" \
         --nucleotide-database "$flag_chocophlan_db" \
         --verbose \
         --remove-temp-output

done < "$flag_manifest"

# Clean up intermediate files
rm -rf "$flag_merge_out"