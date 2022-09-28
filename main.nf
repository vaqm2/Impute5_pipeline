#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def help_msg() {
    log.info """
    Nextflow pipeline to impute data in VCF format using Impute5
    Author: Vivek Appadurai | vivek.appadurai@regionh.dk

    USAGE: nextflow run main.nf

    OPTIONS:

    --ref <ref.vcf.gz> [Gzipp'ed Tab indexed VCF file of reference haplotypes]
    --gt <target.vcf.gz> [Gzipp'ed Tab indexed VCF file of phased target haplotypes]
    --map <chr21.gmap.txt.gz> [Recombination Map]
    --chr <21> [Chromosome to Impute]
    --out <out_prefix> [Prefix to output file]
    --help <Prints this message>
    """
}

if(params.help) {
    help_msg()
    exit 0
}

log.info """
==================================================================================
I M P U T E 5 - P I P E L I N E - N F
==================================================================================
Reference Haplotypes: $params.ref
Target Haplotypes:    $params.gt
Recombination Map:    $params.map
Chromosome:           $params.chr
Output Prefix:        $params.out
==================================================================================
"""

ref_index = params.ref + ".tbi"
gt_index = params.gt + ".tbi"

process chunk_regions {
    label 'low_mem'

    input:
        tuple path(reference_haplotypes),
        path(reference_haplotypes_index),
        path(target_haplotypes),
        path(target_haplotypes_index),
        val(chromosome),
        val(out_prefix),
        path(impute5_chunker)

    output:
        path("${out_prefix}_chr${chromosome}.chunks.csv")

    script:
    """
    ./$impute5_chunker --h $reference_haplotypes \
        --g $target_haplotypes \
        --r $chromosome \
        --o ${out_prefix}_chr${chromosome}.chunks.txt

    echo -e "chr,bufferRegion,imputeRegion" > ${out_prefix}_chr${chromosome}.chunks.csv
    awk \'{print \$2","\$3","\$4}\' ${out_prefix}_chr${chromosome}.chunks.txt >> ${out_prefix}_chr${chromosome}.chunks.csv
    """
}

process impute_chunks {
    label 'big_mem'

    input:
        tuple val(chromosome),
        val(buffer_region),
        val(impute_region),
        path(reference_haplotypes),
        path(reference_haplotypes_index),
        path(target_haplotypes),
        path(target_haplotyoes_index),
        path(map),
        val(out_prefix),
        path(impute5)

    output:
        path("${out_prefix}_${impute_region}.vcf.gz")

    script:
    """
    ./impute5_1.1.5_static --h $reference_haplotypes \
        --g $target_haplotypes \
        --m $map \
        --r $impute_region \
        --buffer-region $buffer_region \
        --out-gp-field \
        --out-ap-field \
        --o ${out_prefix}_${impute_region}.vcf.gz
    """
}

process ligate_chunks {
    label 'low_mem'
    publishDir launchDir

    input:
        tuple path(vcf_files),
        val(out_prefix),
        val(chromosome),
        path(bcftools),
        path(tabix)

    output:
        tuple path("${out_prefix}_chr${chromosome}.vcf.gz"),
        path("${out_prefix}_chr${chromosome}.vcf.gz.tbi")

    script:
    """
    ls ${out_prefix}_${chromosome}*.vcf.gz > chunks_to_ligate.txt
    ./bcftools concat --file-list chunks_to_ligate.txt -Oz -o ${out_prefix}_chr${chromosome}.vcf.gz
    ./tabix -p vcf ${out_prefix}_chr${chromosome}.vcf.gz
    """
}

workflow {
    Channel.fromPath(params.ref) \
    | combine(Channel.fromPath(ref_index)) \
    | combine(Channel.fromPath(params.gt)) \
    | combine(Channel.fromPath(gt_index)) \
    | combine(Channel.of(params.chr)) \
    | combine(Channel.of(params.out)) \
    | combine(Channel.fromPath(params.impute5_chunker_path)) \
    | chunk_regions \
    | splitCsv(header:true) \
    | map{ row -> tuple(row.chr, row.bufferRegion, row.imputeRegion) } \
    | combine(Channel.fromPath(params.ref)) \
    | combine(Channel.fromPath(ref_index)) \
    | combine(Channel.fromPath(params.gt)) \
    | combine(Channel.fromPath(gt_index)) \
    | combine(Channel.fromPath(params.map)) \
    | combine(Channel.of(params.out)) \
    | combine(Channel.fromPath(params.impute5_path)) \
    | impute_chunks \
    | collect \
    | toList \
    | combine(Channel.of(params.out)) \
    | combine(Channel.of(params.chr)) \
    | combine(Channel.fromPath(params.bcftools_path)) \
    | combine(Channel.fromPath(params.tabix_path)) \
    | ligate_chunks
}