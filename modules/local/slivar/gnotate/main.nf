process SLIVAR_GNOTATE {
    container "https://depot.galaxyproject.org/singularity/slivar:0.3.1--h4e814b3_1"
    label 'process_single'

    input:
    path(invcf)
    path(slivar_zip)

    output:
    path(outvcf), emit: vcf

    script:
    outvcf = "${invcf.simpleName}.gnotate.vcf"
    //tmpdir = "$TMPDIR" // This points to /var/temp on Cybertron
    //println tmpdir
    """
    export TMPDIR="temp"
    mkdir temp
    slivar expr -g $slivar_zip -o $outvcf --vcf $invcf
    """
}
