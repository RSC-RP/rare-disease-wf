process SLIVAR_MENDEL_FILTER {
    container 'docker://brentp/slivar:v0.3.1' // has slivar-functions.js
    label 'process_single'
    // Don't publish, but run zip_index_publish afterwards

    input:
    path(invcf)
    path(ped)
    path(slivar_zip)
    val(maf_recessive)
    val(maf_dominant)
    val(maf_denovo)
    val(comphet_nhomalt)
    path(slivar_fn_loose)

    output:
    path(outmnd), emit: mendelian_vcf
    path(outcomphet), emit: comphet_vcf

    script:
    outmnd = "${invcf.simpleName}.mendelian.vcf"
    tempch1 = "${invcf.simpleName}.comphetpre.vcf"
    outcomphet = "${invcf.simpleName}.comphetpost.vcf"
    slivarfn = params.check_gq ? '/opt/slivar/slivar-functions.js' : slivar_fn_loose.name
    //tmpdir = "$TMPDIR" // This points to /var/temp on Cybertron
    //println tmpdir
    cmd =
    """
export TMPDIR="temp"
mkdir temp

slivar expr --vcf $invcf \
    --ped $ped \
    -o $outmnd \
    --pass-only \
    -g $slivar_zip \
    --info 'INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"' \
    --js $slivarfn \
    --family-expr 'denovo:variant.CHROM != "X" && variant.CHROM != "chrX" && fam.every(segregating_denovo) && INFO.gnomad_popmax_af < $maf_denovo' \
    --family-expr 'recessive:variant.CHROM != "X" && variant.CHROM != "chrX" && INFO.gnomad_popmax_af < $maf_recessive && fam.every(segregating_recessive)' \
    --family-expr 'x_denovo:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_denovo_x) && INFO.gnomad_popmax_af < $maf_denovo' \
    --family-expr 'x_recessive:(variant.CHROM == "X" || variant.CHROM == "chrX") && fam.every(segregating_recessive_x)' \
    --family-expr 'dominant:variant.CHROM != "X" && variant.CHROM != "chrX" && INFO.gnomad_popmax_af < $maf_dominant && fam.every(segregating_dominant)'
    """

    if(params.trios_present)
    cmd = cmd +
    """

slivar expr --vcf $invcf \
    --ped $ped \
    -o $tempch1 \
    --pass-only \
    -g $slivar_zip \
    --info 'INFO.gnomad_popmax_af < 0.01 && variant.FILTER == "PASS" && variant.ALT[0] != "*"' \
    --js $slivarfn \
    --family-expr 'denovo:fam.every(segregating_denovo) && INFO.gnomad_popmax_af < $maf_denovo' \
    --trio 'comphet_side:comphet_side(kid, mom, dad) && INFO.gnomad_nhomalt < $comphet_nhomalt'

slivar compound-hets -v $tempch1 \
    --sample-field comphet_side --sample-field denovo -p $ped -o $outcomphet
    """

    else
    cmd = cmd +
    """

    touch $outcomphet
    """

    cmd
}
