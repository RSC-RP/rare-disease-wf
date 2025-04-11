process SLIVAR_DUODEL {
    container "https://depot.galaxyproject.org/singularity/slivar:0.3.1--h4e814b3_1"
    publishDir "$params.outdir", mode: 'copy', overwrite: true
    label 'process_single'

    input:
    path(invcf)
    path(ped)
    path(selfchain)

    output:
    path(outbed)

    script:
    outbed = "${invcf[0].simpleName}.deletions.bed"

    """
    slivar duo-del \
    --ped $ped \
    --exclude $selfchain \
    ${invcf[0]} > $outbed
    """
}
