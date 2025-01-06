process CALL_VARIANTS_TRIO {
    container = "docker://google/deepvariant:deeptrio-1.8.0-gpu"
    maxRetries 2
    maxForks 2 // avoid running out of GPU scratch space
    tag "$meta.id"

    // me_tfrecord is output from make_examples.  role is "child" or "parent".
    input:
    tuple val(meta), path(me_tfrecord), path(gvcf_tfrecord), path(example_info) // meta has id, proband_id, and role

    output:
    tuple val(meta), path("call_variants*.tfrecord.gz"), path(gvcf_tfrecord)

    script:

    sample_id = meta.id
    role = meta.role

    // Build argument for --examples and --outfile (list of call sets to loop through)
    tfr_first = me_tfrecord.findAll{ it.name ==~ /.*tfrecord-00000-of.*/ }
    tfr_prefix = tfr_first.collect{ it.simpleName - "make_examples" }.sort()
    tfr_string = tfr_prefix.join(' ')

    """
    mkdir tmp
    export TMPDIR=tmp
    export TF_GPU_ALLOCATOR=cuda_malloc_async

    for value in $tfr_string
    do
        call_variants \
        --batch_size 384 \
        --outfile "call_variants\${value}_${sample_id}.tfrecord.gz" \
        --examples "make_examples\${value}*.tfrecord@${params.make_examples_nshards}.gz" \
        --checkpoint /opt/models/deeptrio/${params.deepvar_model.toLowerCase()}/${role}
    done
    """

    stub:
    """
    mkdir tmp
    touch call_variants1_${sample_id}.tfrecord.gz
    touch call_variants2_${sample_id}.tfrecord.gz
    """
}
