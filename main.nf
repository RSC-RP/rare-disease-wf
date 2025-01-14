nextflow.enable.dsl=2

include { FILTER_VCF } from './workflows/filter_vcf'
include { CALL_TRIOS } from './workflows/call_trios'

//include { preprocessvcf } from './modules/local/preprocessvcf/main'

workflow {
    FILTER_VCF()
}

workflow calltrios {
    CALL_TRIOS()
}
