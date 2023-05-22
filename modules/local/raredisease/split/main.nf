process SPLIT {
    container = "https://depot.galaxyproject.org/singularity/bcftools:1.9--ha228f0b_4"

    input: tuple(path(vcf), path(tbi))
           path(fai)
    output: path("*.split.*")

    script:
    prefix = vcf.simpleName
    ext = vcf.name - ~/\.gz$/ - ~/.*\./ //vcf.baseName.extension // Should be vcf or gvcf
    """
    # get large chroms and chrM in one file each
    for chrom in \$(awk '\$2 > 40000000 || \$1 ~/(M|MT)\$/' $fai | cut -f 1); do
        bcftools view $vcf --threads $task.cpus -O z -o ${prefix}_\${chrom}.split.${ext}.gz \$chrom
        # Discard files with zero entries
        ENTRIES=\$(zcat ${prefix}_\${chrom}.split.${ext}.gz | grep "^[^#]" | wc -l)
        if [[ "\$ENTRIES" -eq 0 ]]; then
            rm ${prefix}_\${chrom}.split.${ext}.gz
        fi
    done
    # small HLA and gl chroms all to go single file
    awk '!(\$2 > 40000000 || \$1 ~/(M|MT)\$/) { print \$1"\\t0\\t"\$2+1 }' $fai > other_chroms
    # if it's non-empty then create the extras / other split file
    if [ -s other_chroms ]; then
        bcftools view $vcf --threads $task.cpus -R other_chroms -O z -o ${prefix}_OTHER.split.${ext}.gz
        # Discard files with zero entries
        ENTRIES=\$(zcat ${prefix}_OTHER.split.${ext}.gz | grep "^[^#]" | wc -l)
        if [[ "\$ENTRIES" -eq 0 ]]; then
            rm ${prefix}_OTHER.split.${ext}.gz
        fi
    fi
    """

}
