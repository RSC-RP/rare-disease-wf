process CREATE_SEQUENCE_DICTIONARY {
    shell '/bin/bash', '-euo', 'pipefail'
    container 'https://depot.galaxyproject.org/singularity/picard:2.27.4--hdfd78af_0'

    input: path(fasta)
    output: path(outdict), emit: dict

    script:
    basename = fasta.baseName
    if( fasta.extension == "gz" ){
        basename = basename - ~/\.f(a|asta|na)$/
    }
    outdict = "${basename}.dict"
    """
    picard CreateSequenceDictionary -R $fasta -O $outdict
    """
}
