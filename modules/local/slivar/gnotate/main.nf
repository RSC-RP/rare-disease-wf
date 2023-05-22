process SLIVAR_GNOTATE {
    container "https://depot.galaxyproject.org/singularity/slivar:0.2.7--h2eeb373_0"

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
