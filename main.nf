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

process chunk_regions {
    label 'low_mem'

    input:
        tuple path(reference_haplotypes), 
        path(target_haploptypes),
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
    awk '{print \$2","\$3","\$4"}' ${out_prefix}_chr${chromosome}.chunks.txt >> ${out_prefix}_chr${chromosome}.chunks.csv
    """
}

process impute_chunks {
    label 'high_mem'

    input:
        tuple val(chromosome),
        val(buffer_region),
        val(impute_region),
        path(reference_haplotypes), 
        path(target_haploptypes),
        path(map),
        val(out_prefix),
        path(impute5)

    output:
        path("${out_prefix}_${impute_region}.vcf.gz")

    script:
    """
    ./impute5 --h $reference_haplotypes \
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
        tuple path(chunks_to_ligate),
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
    Channel.of(params.ref) \
    | combine(Channel.fromPath(params.gt)) \
    | combine(Channel.fromPath(params.chr)) \
    | combine(Channel.of(params.out)) \
    | combine(Channel.fromPath(params.impute5_chunker)) \
    | chunk_regions \
    | splitCsv(header:true) \
    | map{row -> tuple(row.chr, row.bufferRegion, row.imputeRegion)} \
    | combine(Channel.fromPath(params.ref)) \
    | combine(Channel.fromPath(params.gt)) \
    | combine(Channel.fromPath(params.map)) \
    | combine(Channel.of(params.out)) \
    | combine(Channel.fromPath(params.impute5_path)) \
    | impute_chunks \
    | collect \
    | combine(Channel.of(params.out)) \
    | combine(Channel.of(params.chr)) \
    | combine(Channel.fromPath(bcftools_path)) \
    | combine(Channel.fromPath(tabix_path)) \
    | ligate_chunks
}