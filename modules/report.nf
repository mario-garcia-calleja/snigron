process VARIANT_REPORT {
    label 'process_low'
    publishDir "${params.outdir}/report", mode: params.publish_mode
    conda "${projectDir}/environment.yml"

    input:
    path snps
    path regions
    path comparisons

    output:
    path 'snigron_variants_vs_reference.tsv'
    path 'snigron_comparisons.tsv'
    path 'snigron_summary.tsv'
    path 'snigron_report.xlsx'

    script:
    """
    snigron_variants.py \\
        --snps ${snps} \\
        --regions ${regions} \\
        --comparisons ${comparisons} \\
        --flank ${params.flank}
    """

    stub:
    """
    touch snigron_variants_vs_reference.tsv snigron_comparisons.tsv
    touch snigron_summary.tsv snigron_report.xlsx
    """
}
