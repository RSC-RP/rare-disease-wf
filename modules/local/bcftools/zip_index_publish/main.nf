process ZIP_INDEX_PUBLISH {
    container "https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0"
    shell '/bin/bash', '-euo', 'pipefail'
    publishDir "$params.outdir", mode: 'copy', overwrite: true
    label 'process_single'

    input:
        path(invcf)
        val(labl)
        val(sortview) // string 'sort' or 'view'. Only use 'sort' if sorting is needed.
    output:
        tuple path(outvcf), path(outindex)

    script:
    outvcf = "${invcf.simpleName}.${labl}.vcf.gz"
    outindex = "${invcf.simpleName}.${labl}.vcf.gz.csi"
    if(sortview == "view"){
        extargs = "--threads ${task.cpus}"
    }
    else { // 'sort' does not have a multithreading option
        extargs = "-T temp"
    }
    """
mkdir temp
bcftools $sortview $invcf -Oz -o $outvcf $extargs
bcftools index --threads $task.cpus $outvcf
    """
}
