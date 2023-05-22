process SNPSIFT_ANNOTATE {
    label 'process_medium'

    //conda (params.enable_conda ? "bioconda::snpeff=5.1" : null)
    //container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
    //    'https://depot.galaxyproject.org/singularity/snpeff:5.1--hdfd78af_2' :
    //    'quay.io/biocontainers/snpeff:5.1--hdfd78af_2' }"
    container 'shub://qbicsoftware/qbic-singularity-snpeff:latest'

    input:
    path(vcf)
    path(annotations_vcf)
    val(fields)

    output:
    path("*.snpsift.vcf"), emit: vcf
    path "versions.yml"  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def avail_mem = 6
    if (!task.memory) {
        log.info '[snpEff] Available memory not known - defaulting to 6GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = task.memory.giga
    }
    def prefix = vcf.simpleName
    """
    java -Xmx${avail_mem}g \\
        -jar /usr/local/lib/snpEff/SnpSift.jar annotate \\
        -info $fields \\
        $args \\
        ${annotations_vcf[0]} \\
        $vcf \\
        > ${prefix}.snpsift.vcf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        snpeff: \$(echo \$(snpEff -version 2>&1) | cut -f 2 -d ' ')
    END_VERSIONS
    """
}
