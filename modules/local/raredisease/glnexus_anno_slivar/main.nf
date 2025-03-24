process GLNEXUS_ANNO_SLIVAR {
    container 'docker://brentp/rare-disease:v0.2.2'
    shell '/bin/bash', '-euo', 'pipefail'
    publishDir "${params.outdir}/joint-by-chrom/", mode: 'copy'

    input:
        tuple(path(gvcfs), val(chrom))
        path(fasta)
        path(fai)
        path(gff)
        path(slivar_zip)
        val(cohort_name)

    output: tuple(path(output_file), path(output_csi))

    script:
        output_file = "${cohort_name}.${chrom}.glnexus.anno.bcf"
        output_csi = "${output_file}.csi"
        file("$workDir/file.list.${cohort_name}.${chrom}").withWriter { fh ->
            gvcfs.each { gvcf ->
                fh.write(gvcf.toString()); fh.write("\n")
            }
        }
        """
# GRCh38.99
# GRCh37.75

glnexus_cli \
    -t ${task.cpus} \
    --mem-gbytes 128 \
    --config DeepVariant${params.model_type} \
    --list $workDir/file.list.${cohort_name}.${chrom} \
| bcftools norm --threads 3 -m - -w 10000 -f $fasta -O u \
| bcftools csq --threads 3 -s - --ncsq 50 -g $gff -l -f $fasta - -o - -O v \
| snpEff eff -noStats -dataDir $projectDir GRCh38.99 \
| slivar expr -g $slivar_zip -o $output_file --vcf -


bcftools index --threads 6 $output_file
    """

}
