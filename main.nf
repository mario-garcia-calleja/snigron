#!/usr/bin/env nextflow
/*
 * SNIGRON: SNPs, indels and integrons in bacterial assemblies
 *
 *   assembly FASTA -> dnadiff vs reference -> variants in reference coordinates
 *                  -> per-comparison acquired / lost variants + region annotation
 *   assembly FASTA -> IntegronFinder (optionally with functional HMM annotation)
 */

include { DNADIFF         } from './modules/mummer'
include { INTEGRON_FINDER } from './modules/integron_finder'
include { VARIANT_REPORT  } from './modules/report'

def requireFile(String label, value) {
    if (!value) error "Missing required parameter --${label}"
    def f = file(value)
    if (!f.exists()) error "--${label}: file not found: ${value}"
    return f
}

// Relative paths in a samplesheet are resolved against the samplesheet's own folder
def sheetFile(sheet, String path) {
    def f = (path.startsWith('/') || path =~ /^[a-zA-Z0-9]+:\/\//) ? file(path) : file(sheet.parent.resolve(path).toString())
    if (!f.exists()) error "File listed in ${sheet.name} not found: ${path}"
    return f
}

workflow {
    log.info """
    S N I G R O N   ~   variants and integrons in bacterial assemblies
    ================================================================
    samplesheet : ${params.samplesheet}
    comparisons : ${params.comparisons}
    reference   : ${params.reference}
    regions     : ${params.regions}
    func HMMs   : ${params.func_hmm_dir ?: '- (integrons without functional annotation)'}
    outdir      : ${params.outdir}
    """.stripIndent()

    def reference   = requireFile('reference', params.reference)
    def regions     = requireFile('regions', params.regions)
    def samplesheet = requireFile('samplesheet', params.samplesheet)
    def comparisons = requireFile('comparisons', params.comparisons)

    assemblies = Channel.fromPath(samplesheet)
        .splitCsv(header: true, strip: true)
        .map { r ->
            if (!r.sample || !r.fasta) error "Samplesheet needs 'sample' and 'fasta' columns: ${r}"
            tuple(r.sample, sheetFile(samplesheet, r.fasta))
        }

    // Fail early if a comparison references a sample that is not in the samplesheet
    known = assemblies.map { it[0] }.collect().map { it as Set }
    Channel.fromPath(comparisons)
        .splitCsv(header: true, strip: true)
        .combine(known)
        .subscribe { row, names ->
            [row.control, row.treated].each { s ->
                if (!(s in names)) error "Comparison '${row.comparison}': unknown sample '${s}'"
            }
        }

    DNADIFF(assemblies, reference)
    VARIANT_REPORT(DNADIFF.out.snps.map { it[1] }.collect(), regions, comparisons)

    if (!params.skip_integrons) {
        def func_db = params.func_hmm_dir ? file(params.func_hmm_dir, checkIfExists: true) : []
        INTEGRON_FINDER(assemblies, func_db)
    }
}
