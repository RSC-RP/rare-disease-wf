process TABLE_ANNOVAR {
    container "docker://bioinfochrustrasbourg/annovar"
    label 'process_single'

    input:
    path(invcf)
    path(annovar_db)
    val(annovar_buildver)

    output:
    path(outvcf), emit: vcf

    script:
    prefix = "${invcf.simpleName}"
    outvcf = prefix + "." + annovar_buildver + "_multianno.vcf"
    """
    mkdir temp
    table_annovar.pl -thread ${task.cpus} \
    -vcfinput \
    -buildver $annovar_buildver \
    -outfile $prefix \
    -tempdir temp/ \
    -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a \
    -operation gx,r,f,f,f \
    -polish \
    -nastring . \
    $invcf $annovar_db
    """
}
