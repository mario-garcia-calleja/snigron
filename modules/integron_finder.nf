process INTEGRON_FINDER {
    tag "${sample}"
    label 'process_medium'
    publishDir "${params.outdir}/integrons", mode: params.publish_mode
    conda "bioconda::integron_finder=2.0.6"

    input:
    tuple val(sample), path(assembly)
    path func_db                 // directory of pressed HMMs (e.g. Resfams); empty list = none

    output:
    path "Results_Integron_Finder_*", emit: results

    script:
    def func = func_db ? "--func-annot --path-func-annot ${func_db}" : ''
    """
    integron_finder ${assembly} \\
        --local-max \\
        --cpu ${task.cpus} \\
        ${func}
    """

    stub:
    "mkdir -p Results_Integron_Finder_${sample} && touch Results_Integron_Finder_${sample}/${sample}.summary"
}
