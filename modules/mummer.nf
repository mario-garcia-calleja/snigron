process DNADIFF {
    tag "${sample}"
    label 'process_low'
    publishDir "${params.outdir}/dnadiff", mode: params.publish_mode
    conda "bioconda::mummer=3.23"

    input:
    tuple val(sample), path(assembly)
    path reference

    output:
    tuple val(sample), path("${sample}.snps"), emit: snps
    path "${sample}.report",                   emit: report
    path "${sample}.{1coords,delta}"

    script:
    """
    dnadiff -p ${sample} ${reference} ${assembly}
    """

    stub:
    "touch ${sample}.snps ${sample}.report ${sample}.1coords ${sample}.delta"
}
