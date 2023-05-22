process BCFTOOLS_NORM_CSQ {
    container = "https://depot.galaxyproject.org/singularity/bcftools:1.17--haef29d1_0"
    shell = ['/bin/bash', '-euo', 'pipefail']

    input:
        path(invcf)
        path(fasta)
        path(fai)
        path(gff)
        path(chrom_convert)
    output: path(outvcf), emit: vcf

    script:
    outvcf = "${invcf.simpleName}.bcftools.vcf"

    if(params.chromnames == 'ucsc')
    """
    bcftools annotate --rename-chrs $chrom_convert --threads $task.cpus -O u $invcf \
| bcftools norm --threads $task.cpus -m - -w 10000 -c wx -f $fasta -O u \
| bcftools csq -s - --ncsq 50 -g $gff -l -f $fasta - -o $outvcf -O v
    """

    else
    """
bcftools norm --threads $task.cpus -m - -w 10000 -f $fasta -O u $invcf \
| bcftools csq -s - --ncsq 50 -g $gff -l -f $fasta - -o $outvcf -O v
    """
}
