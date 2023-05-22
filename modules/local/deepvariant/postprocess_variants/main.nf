process POSTPROCESS_VARIANTS {
    container = "docker://google/deepvariant:deeptrio-1.5.0-gpu"
    publishDir "${params.outdir}/gvcfs/", mode: 'copy', overwrite: true

    input:
    tuple path(cv_tfrecord), path(gvcf_tfrecord), val(sample_id)
    path(fasta_ensembl)
    path(fai_ensembl)

    output:
    tuple path(out_gvcf), path("${out_gvcf}.tbi")

    script:
    out_gvcf = "${sample_id}.gvcf.gz"
    out_vcf = "${sample_id}.vcf.gz"

    // Build argument for --nonvariant_site_tfrecord_path and --infile
    tfr_prefix = cv_tfrecord.collect{ it.simpleName - "call_variants" - "_${sample_id}"}.sort()
    tfr_string = tfr_prefix.join(' ')

    """
    mkdir tmp
    export TMPDIR=tmp

    for value in $tfr_string
    do
        postprocess_variants \
        --ref $fasta_ensembl \
        --infile "call_variants\${value}_${sample_id}.tfrecord.gz" \
        --outfile temp.vcf.gz \
        --nonvariant_site_tfrecord_path "gvcf\${value}.tfrecord@${params.make_examples_nshards}.gz" \
        --gvcf_outfile temp.gvcf.gz

        if test -f "$out_gvcf"; then
            mv $out_gvcf temp1.gvcf.gz
            bcftools concat -O z -o $out_gvcf --threads ${task.cpus} temp1.gvcf.gz temp.gvcf.gz
        else
            mv temp.gvcf.gz $out_gvcf
        fi
    done

    bcftools index -t --threads ${task.cpus} $out_gvcf
    """

    stub:
    """
    mkdir tmp
    touch "${sample_id}.gvcf.gz"
    touch "${sample_id}.gvcf.gz.tbi"
    touch temp.gvcf.gz
    touch temp.gvcf.gz.tbi
    touch temp1.gvcf.gz
    touch temp.vcf.gz
    touch temp.vcf.gz.tbi
    touch temp.visual_report.html
    """
}
