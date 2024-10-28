process POP_STRUCT {
    container 'vcf_umap_v0.4.sif'
    publishDir "${params.outdir}/qc", mode: 'copy', overwrite: true
    label 'process_single'

    input:
    tuple path(invcf), path(index)
    val(txdb)
    path(sample_metadata)
    val(metadata_columns)
    val(seed)
    val(min_call_rate)
    val(n_exons)
    val(n_pcs)

    output:
    path("*.tsv"), emit: coords
    path("*.tiff"), emit: graphs

    script:
    """
Rscript /run_umap.R --invcf $invcf \
    --txdb $txdb \
    --sample_metadata $sample_metadata \
    --metadata_columns $metadata_columns \
    --seed $seed \
    --min_call_rate $min_call_rate \
    --n_exons $n_exons \
    --n_pcs $n_pcs
    """
}
