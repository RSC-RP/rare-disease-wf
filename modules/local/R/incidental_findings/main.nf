process INCIDENTAL_FINDINGS {
    container 'incidental_findings.sif'
    publishDir "${params.outdir}/incidental", mode: 'copy', overwrite: true
    label 'process_single'

    input:
    path(invcf)
    val(txdb)
    val(buildver)
    path(sample_list)

    output:
    path("*.incidental.txt")

    script:
    """
    Rscript /opt/incidental_findings/incidental.R \
    --invcf $invcf \
    --txdb $txdb \
    --build $buildver \
    --sample_list $sample_list
    """
}
