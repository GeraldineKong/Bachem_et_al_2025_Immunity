#!/usr/bin/env nextflow


params.manifest = "manifest.txt" // TSV file with 3 columns: sample_name, absolute_path_R1, absolute_path_R2
params.outdir = "metagenome_assembled_genome"
params.threads = 16

params.checkm2_db = "checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"
params.gunc_db = "gunc_db/gunc_db_progenomes2.1.dmnd"

Channel
    .fromPath(params.manifest)
    .splitCsv(header: false, sep: '\t')
    .map { row -> tuple(row[0], file(row[1]), file(row[2])) }
    .set { sample_reads }

process AssembleSample {
    tag "$sample"

    input:
    tuple val(sample), path(r1), path(r2)

    output:
    tuple(val(sample), path("assembly_${sample}/contigs.fasta")), emit: assemblies

    script:
    """
    if [ ! -f "$r1" ] || [ ! -f "$r2" ]; then
        echo "Missing input FASTQ files for $sample"
        exit 1
    fi

    /home/kongg1/tools/SPAdes-3.15.5-Linux/bin/metaspades.py \
        -1 $r1 -2 $r2 \
        -o assembly_${sample} \
        -t ${params.threads} \
        --checkpoints last
    """
}

process AlignToAssembly {
    tag "$sample"

    input:
    tuple val(sample), path(assembly), path(r1), path(r2) from aligned_input

    output:
    tuple(val(sample), path(assembly), path("${sample}.bam")), emit: sample_bams

    script:
    """
    bowtie2-build $assembly ${sample}_index

    bowtie2 -x ${sample}_index -1 $r1 -2 $r2 -p ${params.threads} | \
      samtools sort -@ ${params.threads} -O BAM -o ${sample}.bam
    """
}

process Binning {
    tag "$sample"

    input:
    tuple val(sample), path(assembly), path(bam)

    output:
    path "bins_${sample}", emit: sample_bins

    script:
    """
    vamb --outdir bins_${sample} \
         --fasta $assembly \
         --bamfiles $bam \
         -p ${params.threads} -m 1000 --minfasta 1000
    """
}

process CombineBins {
    tag "combine_all_bins"

    input:
    path sample_bins.collect()

    output:
    path "combined_bins", emit: combined_bins

    script:
    """
    mkdir -p combined_bins
    find bins_* -name "*.fa" -exec cp {} combined_bins/ \;
    """
}

process CheckM2 {
    tag "checkm2"

    input:
    path combined_bins

    output:
    path "checkm2_out", emit: checkm_out

    script:
    """
    checkm2 predict \
      --threads ${params.threads} \
      --input combined_bins \
      --output-directory checkm2_out \
      --remove_intermediates -x .fa \
      --database_path ${params.checkm2_db}
    """
}

process RunGUNC {
    tag "gunc"

    input:
    path combined_bins

    output:
    path "gunc_out", emit: gunc_out

    script:
    """
    gunc run --input_dir combined_bins \
             --db_file ${params.gunc_db} \
             --threads ${params.threads} \
             --out_dir gunc_out
    """
}

process FilterGoodBins {
    tag "filter_good_bins"

    input:
    path checkm2_dir
    path gunc_dir
    path bin_dir

    output:
    path "good_bins", emit: good_bins

    script:
    """
    mkdir -p good_bins

    # Merge Checkm2 and GUNC results, filter for good bins
    python3 bin/gunc_checkm2_merge.py \
      --checkm2 $checkm2_dir/quality_report.tsv \
      --gunc $gunc_dir/GUNC*.tsv \
      --output passing_bins.tsv \
      --comp 90 \
      --cont 5

    # Extract genome names and copy .fa files from bin_dir to good_bins/
    tail -n +2 passing_bins.tsv | cut -f1 | while read bin; do
        cp $bin_dir/\${bin}.fa good_bins/
    done
    """
}

process Dereplicate {
    tag "dRep"

    input:
    path good_bins

    output:
    path "drep_out/dereplicated_genomes", emit: drep_genomes_dir

    script:
    """
    dRep dereplicate drep_out \
      --ignoreGenomeQuality \
      -p ${params.threads} \
      -g $good_bins/*.fa
    """
}

process CombineDereplicatedMAGs {
    tag "combine_and_index_MAGs"

    input:
    path drep_genomes_dir

    output:
    path "mag_index.*"

    script:
    """
    cat $drep_genomes_dir/*.fa > dereplicated_MAGs.fa

    bowtie2-build dereplicated_MAGs.fa mag_index
    """
}

process QuantMAGs {
    tag "$sample"

    input:
    tuple val(sample), path(r1), path(r2)
    path "mag_index.*"

    output:
    path "quant_results/${sample}.tsv", emit: quant_profiles

    script:
    """
    mkdir -p bams quant_results

    bowtie2 -x mag_index -1 $r1 -2 $r2 -p ${params.threads} | \
        samtools sort -@ ${params.threads} -O BAM -o bams/${sample}.bam

    samtools view bams/${sample}.bam | \
        woltka classify -i - -o quant_results/${sample}.tsv --to-tsv --unassigned
    """
}

workflow {
    // Assemble and align
    assembly_dirs = AssembleSample(sample_reads)
    aligned_input = assemblies
        .map { sample, assembly -> tuple(sample, [assembly]) }
        .join(sample_reads.map { sample, r1, r2 -> tuple(sample, [r1, r2]) })
        .map { sample, a, r -> tuple(sample, a[0], r[0], r[1]) }
    sample_bams = AlignToAssembly(aligned_input)

    // Binning 
    sample_bins = Binning(sample_bams)

    // Combine bins and QC
    combined_bins = CombineBins(sample_bins)
    checkm = CheckM2(combined_bins)
    gunc   = RunGUNC(combined_bins)

    // Filter good bins
    good_bins = FilterGoodBins(checkm, gunc, combined_bins)

    // Dereplicate
    drep_out = Dereplicate(good_bins)

    // Combine dereplicated MAGs and create index
    mag_indexed = CombineDereplicatedMAGs(drep_out)

    // Generate abundance table
    quant_input = sample_reads.map { sample, r1, r2 -> tuple(sample, r1, r2) }.combine(mag_indexed)
    quant_profiles = QuantMAGs(quant_input)

}
