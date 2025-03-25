process LIFTOVER_HG19_TO_HG38 {
    shell '/bin/bash', '-euo', 'pipefail'
    container 'https://depot.galaxyproject.org/singularity/picard:2.27.4--hdfd78af_0'

    input:
        path(invcf)
        path(fasta)
        path(dict)
    output: path(outvcf), emit: vcf

    script:
    outvcf = "${invcf.simpleName}.liftover.vcf"
    println "Performing liftover"
    println "Output file: ${outvcf}"
    """
    wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz

    picard -Xmx${params.picard_memory_gb}g LiftoverVcf \
    -C hg19ToHg38.over.chain.gz \
    -I $invcf -O $outvcf \
    --REFERENCE_SEQUENCE $fasta \
    --REJECT rejected.vcf \
    --MAX_RECORDS_IN_RAM $params.picard_max_records_in_ram
    """
}
