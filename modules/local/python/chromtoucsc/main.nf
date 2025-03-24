process CHROMNAMES_TO_UCSC{
    // chromToUcsc is in the bin folder
    shell '/bin/bash', '-euo', 'pipefail'
    container 'docker://python:latest'

    input: path(invcf)
    output: path(outvcf), emit: vcf

    script:
    outvcf = "${invcf.simpleName}.UCSC.vcf"
    println "Converting to UCSC chromosome names"
    println "Output file: ${outvcf}"
    """
    wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/analysisSet/chromAlias/${params.chromnames}ToUcsc.txt
    chromToUcsc -i $invcf -o $outvcf -a ${params.chromnames}ToUcsc.txt
    """
}
