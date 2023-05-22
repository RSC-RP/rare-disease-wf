process CLINVAR_DOWNLOAD {
    input:
    val(build)

    output:
    tuple path('clinvar.vcf.gz'), path('clinvar.vcf.gz.tbi')

    script:
    if(build == 'hg19'){
        assem = 'GRCh37'
    }
    if(build == 'hg38'){
        assem = 'GRCh38'
    }

    """
    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${assem}/clinvar.vcf.gz
    wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_${assem}/clinvar.vcf.gz.tbi
    """
}
