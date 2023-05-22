process SLIVAR_DUODEL {
    container "https://depot.galaxyproject.org/singularity/slivar:0.2.7--h2eeb373_0"
    publishDir "$params.outdir", mode: 'copy', overwrite: true

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
