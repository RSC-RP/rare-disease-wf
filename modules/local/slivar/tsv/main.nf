process SLIVAR_TSV {
    container 'docker://brentp/rare-disease:v0.2.2'
    publishDir "$params.outdir", mode: 'copy', overwrite: true

    input:
    path(vcf_mendelian)
    path(vcf_comphet)
    path(ped)

    output:
    path(slivar_tsv)

    exec:
    slivar_tsv = "${vcf_mendelian[0].simpleName}.slivar.tsv"

    script:
    // Don't filter by impactful and genic since we already did pathogenicity filter.
    cmd =
    """
gunzip -c ${vcf_mendelian[0]} |
slivar tsv \
  -s denovo \
  -s x_denovo \
  -s recessive \
  -s x_recessive \
  -s dominant \
  -i gnomad_popmax_af -i gnomad_popmax_af_filter -i gnomad_nhomalt \
  -c ANN \
  -g /opt/slivar/pli.lookup \
  -g /opt/slivar/clinvar_gene_desc.txt \
  -p $ped > $slivar_tsv
    """

    if(params.trios_present)
    cmd = cmd +
    """
gunzip -c ${vcf_comphet[0]} |
slivar tsv \
  -s slivar_comphet \
  -i gnomad_popmax_af -i gnomad_popmax_af_filter -i gnomad_nhomalt \
  -c ANN \
  -g /opt/slivar/pli.lookup \
  -g /opt/slivar/clinvar_gene_desc.txt \
  -p $ped | { grep -v ^# || true; } >> $slivar_tsv # || true avoids error if there are no compound hets.
    """

    cmd
}
