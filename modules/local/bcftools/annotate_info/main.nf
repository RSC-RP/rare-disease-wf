process BCFTOOLS_ANNOTATE_INFO {
    container "https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0"
    shell '/bin/bash', '-euo', 'pipefail'

    input:
        path(vcf)      // Input VCFs
        path(anno_vcf) // Bgzipped VCF with annotations
        path(anno_tbi) // Index for the annotations
        val(columns)   // comma-separated list of info columns to transfer
    output: path(outvcf), emit: vcf

    exec:
    prefix = vcf.simpleName
    outvcf = "${prefix}.bcfanno.vcf"

    script:
    """
# see https://github.com/samtools/bcftools/issues/1199#issuecomment-624713745
bcftools query -f'%CHROM\t%POS\n' $vcf > sites.txt

bcftools view -o temp.vcf.gz -Oz $vcf
bcftools index temp.vcf.gz

bcftools annotate \
--threads $task.cpus \
-R sites.txt \
--annotations $anno_vcf \
--columns $columns \
--output $outvcf -Ov \
temp.vcf.gz
    """
}
