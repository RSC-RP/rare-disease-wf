nextflow.enable.dsl=2

include { BCFTOOLS_NORM_CSQ } from './modules/local/bcftools/norm_csq/main'
include { SNPEFF } from './modules/nf-core/snpeff/main'
include { SLIVAR_GNOTATE } from './modules/local/slivar/gnotate/main'
include { TABLE_ANNOVAR } from './modules/local/annovar/main'
include { MVP_ANNO } from './modules/local/python/mvp/main'
include { ZIP_INDEX_PUBLISH } from './modules/local/bcftools/zip_index_publish/main'
include { SPLIT } from './modules/local/raredisease/split/main'
include { BCFTOOLS_CONCAT } from './modules/local/bcftools/concat/main'
include { CHROMNAMES } from './modules/local/python/chromnames/main'
include { FILEORDER } from './modules/local/python/fileorder/main'
include { SLIVAR_MENDEL_FILTER } from './modules/local/slivar/mendel_filter/main'
include { CLINVAR_DOWNLOAD } from './modules/local/clinvar_download/main'
include { SNPSIFT_ANNOTATE } from './modules/local/snpsift/annotate/main'
include { FILTER_PATHOGENIC } from './modules/local/python/filter_pathogenic/main'
include { POP_STRUCT } from './modules/local/R/pop_struct/main'
include { SLIVAR_TSV } from './modules/local/slivar/tsv/main'
include { INCIDENTAL_FINDINGS } from './modules/local/R/incidental_findings/main'
include { MAKE_EXAMPLES_TRIO } from './modules/local/deepvariant/make_examples_trio/main'
include { CALL_VARIANTS_TRIO } from './modules/local/deepvariant/call_variants_trio/main'
include { MAKE_EXAMPLES_SINGLE } from './modules/local/deepvariant/make_examples_single/main'
include { CALL_VARIANTS_SINGLE } from './modules/local/deepvariant/call_variants_single/main'
include { POSTPROCESS_VARIANTS } from './modules/local/deepvariant/postprocess_variants/main'
include { GLNEXUS } from './modules/local/glnexus/main'
include { SLIVAR_DUODEL } from './modules/local/slivar/duodel/main'
include { BWA_MEM } from './modules/nf-core/bwa/mem/main'
include { SAMTOOLS_INDEX } from './modules/nf-core/samtools/index/main'
//include { preprocessvcf } from './modules/local/preprocessvcf/main'

// Channel setup for Ensembl references for BcfTools and Slivar
Channel.fromPath(file(params.fasta_ensembl, checkIfExists: true))
    .collect() // allows the reference to be used with multiple input VCFs
    .set{ fasta_ensembl }
Channel.fromPath(file(params.gff_ensembl, checkIfExists: true))
    .collect()
    .set{ gff_ensembl }
// FAIDX for reference genome
fai_file = file("${params.fasta_ensembl}.fai")
if( fai_file.exists() ){
    Channel.fromPath(fai_file)
        .collect()
        .set{ fai_ensembl }
}
else {
    // Run Samtools Faidx
    SAMTOOLS_FAIDX(fasta_ensembl)
    SAMTOOLS_FAIDX.out.fai
        .collect()
        .set{ fai_ensembl }
}

// Reference sequence for alignment and genotype calling, which may be different from Ensembl.
Channel.fromPath(file(params.fasta_bams, checkIfExists: true))
    .collect() // allows the reference to be used with multiple input VCFs
    .set{ fasta_bams }
fai_file2 = file("${params.fasta_bams}.fai")
if( fai_file2.exists() ){
    Channel.fromPath(fai_file2)
        .collect()
        .set{ fai_bams }
}
else {
    // Run Samtools Faidx
    SAMTOOLS_FAIDX(fasta_bams)
    SAMTOOLS_FAIDX.out.fai
        .collect()
        .set{ fai_bams }
}

// Slivar file with gnomAD data
Channel.fromPath(file(params.slivar_zip, checkIfExists: true))
    .collect()
    .set{ slivar_zip }
// Slivar selfchain file for duo-del
Channel.fromPath(file(params.slivar_selfchain, checkIfExists: true))
    .collect()
    .set{ slivar_selfchain }

// Chromsome conversion file
Channel.fromPath(file(params.chrom_convert, checkIfExists: true))
    .collect()
    .set{ chrom_convert }

// VCF
Channel.fromPath(file(params.vcf, checkIfExists: true))
    .set{ vcf_ch }
// Pedigree
Channel.fromPath(file(params.ped, checkIfExists: true))
    .collect()
    .set{ ped_ch }
// Metadata
Channel.fromPath(file(params.metadata_tab, checkIfExists: true))
    .set{metadata_ch}

// workflow to sort and index vcf and split it up by chromosome
workflow index_split{
    take:
        vcf_ch
    emit:
        vcf_ch
        vcf_ch_sp
    main:
        if(params.sort){
            sortview = 'sort'
        }
        else{
            sortview = 'view'
        }
        ZIP_INDEX_PUBLISH(vcf_ch, 'original', sortview)
        ZIP_INDEX_PUBLISH.out
            .set{ vcf_ch }
        SPLIT(ZIP_INDEX_PUBLISH.out, fai_bams)
        SPLIT.out
            .set{ vcf_ch_sp }
}

// workflow for QC
workflow qc {
    take:
        vcf_ch // tuple of VCF and index
        metadata_ch
    main:
        POP_STRUCT(vcf_ch,
            "TxDb.Hsapiens.UCSC.${params.annovar_buildver}.knownGene",
            metadata_ch,
            params.metadata_columns,
            params.umap_seed,
            params.umap_min_call_rate,
            params.umap_n_exons,
            params.umap_n_pcs)
}

// workflow to merge VCFs into one file
workflow merge_vcfs {
    take:
        vcf_ch_sp
    emit:
        vcf_ch
    main:
        CHROMNAMES(fai_ensembl)
        FILEORDER(vcf_ch_sp, CHROMNAMES.out)
        BCFTOOLS_CONCAT(vcf_ch_sp, FILEORDER.out, file(params.vcf).simpleName)
            .set{ vcf_ch }
}

// workflow for annotation
workflow annotate {
    take:
        vcf_ch_sp
        slivar_zip
    emit:
        vcf_ch_sp
    main:
        // Convert from list to channel emitting individual items
        vcf_ch_sp
            .flatten()
            .set{ vcf_ch_sp }

        // Annotation steps
        BCFTOOLS_NORM_CSQ(vcf_ch_sp, fasta_ensembl, fai_ensembl, gff_ensembl, chrom_convert)
        genome = params.gff_ensembl - ~/.*Homo_sapiens\./ - ~/\.gff.*/
        if(params.annovar_buildver == "hg38" && ! file("${params.snpEff_dir}/${genome}").exists()){
            genome = "GRCh38.105" // most recent snpEff build
        }
        SNPEFF(BCFTOOLS_NORM_CSQ.out.vcf, genome, params.snpEff_dir)
        SLIVAR_GNOTATE(SNPEFF.out.vcf, slivar_zip)
        Channel.fromPath(file(params.annovar_db, checkIfExists: true))
            .collect()
            .set{ annovar_db }
        TABLE_ANNOVAR(SLIVAR_GNOTATE.out.vcf, annovar_db, params.annovar_buildver)
        TABLE_ANNOVAR.out.vcf
            .set{ vcf_ch_sp }
        if(params.annovar_buildver == "hg19"){
            Channel.fromPath(file(params.mvp_tab, checkIfExists: true))
                .collect()
                .set{ mvp_tab }
            MVP_ANNO(vcf_ch_sp, mvp_tab)
            MVP_ANNO.out.vcf
                .set{ vcf_ch_sp }
        }
        CLINVAR_DOWNLOAD(params.annovar_buildver)
        CLINVAR_DOWNLOAD.out
            .collect()
            .set{ clinvar_vcf }
        SNPSIFT_ANNOTATE(vcf_ch_sp, clinvar_vcf, 'ALLELEID,CLNSIG,CLNDN')
        SNPSIFT_ANNOTATE.out.vcf
            .set{ vcf_ch_sp }

        // Convert back from individual items to list
        vcf_ch_sp
            .collect()
            .set{ vcf_ch_sp }

        if(params.export_prefilt){
            // Concatenate individual VCFs into one
            merge_vcfs(vcf_ch_sp)
            // Export annotated VCF before filtering
            ZIP_INDEX_PUBLISH(merge_vcfs.out, 'annotated', 'view')
        }
}

// Short workflows for publishing after filtering
workflow publish_mendelian {
    take: vcf_ch_sp
    emit: vcf_ch

    main:
    // Concatenate individual VCFs into one
    merge_vcfs(vcf_ch_sp)
    // Export annotated VCF after filtering
    ZIP_INDEX_PUBLISH(merge_vcfs.out, 'mendelian', 'view')
    ZIP_INDEX_PUBLISH.out
        .set{ vcf_ch }
}

workflow publish_comphet {
    take: vcf_ch_sp
    emit: vcf_ch

    main:
    if(params.trios_present){
        // Concatenate individual VCFs into one
        merge_vcfs(vcf_ch_sp)
        // Export annotated VCF after filtering
        ZIP_INDEX_PUBLISH(merge_vcfs.out, 'comphet', 'view')
        ZIP_INDEX_PUBLISH.out
            .set{ vcf_ch }
    } else {
        Channel.fromPath('dummy.vcf.gz')
            .set{ vcf_ch }
    }
}

// subworkflow for incidental findings
workflow incidental {
    take:
        vcf_ch_sp
        samples_ch
    main:
        merge_vcfs(vcf_ch_sp)
        INCIDENTAL_FINDINGS(merge_vcfs.out,
            "TxDb.Hsapiens.UCSC.${params.annovar_buildver}.knownGene",
            params.annovar_buildver,
            samples_ch)
}

// subworkflow for filtering
workflow filtervcf {
    take:
        vcf_ch_sp
        ped_ch
        slivar_zip
    emit:
        vcf_ch_mendelian
        vcf_ch_comphet
    main:
        // Convert from list to channel emitting individual items
        vcf_ch_sp
            .flatten()
            .set{ vcf_ch_sp }
        
        // Filtering steps and incidental findings
        FILTER_PATHOGENIC(vcf_ch_sp)
        if(params.do_incidental){
            Channel.fromPath(file(params.incidental_samples), checkIfExists: true)
                .set{ samples_ch }
            incidental(FILTER_PATHOGENIC.out.collect(), samples_ch)
        }
        SLIVAR_MENDEL_FILTER(FILTER_PATHOGENIC.out, ped_ch, slivar_zip,
            params.maf_recessive, params.maf_dominant, params.maf_denovo, params.comphet_nhomalt)
            .set{ vcf_ch_sp }
        
        // Convert back from individual items to list
        vcf_ch_sp.mendelian_vcf
            .collect()
            .set{ vcf_ch_sp_mendelian }
        vcf_ch_sp.comphet_vcf
            .collect()
            .set{ vcf_ch_sp_comphet }

        // Export
        publish_mendelian(vcf_ch_sp_mendelian)
        publish_mendelian.out
            .set{ vcf_ch_mendelian }
        publish_comphet(vcf_ch_sp_comphet)
        publish_comphet.out
            .set{ vcf_ch_comphet }
}

// main workflow
workflow {
    //preprocessvcf(vcf_ch)
    //preprocessvcf.out.view()
    index_split(vcf_ch)
    //index_split.out.view()
    qc(index_split.out[0], metadata_ch)
    SLIVAR_DUODEL(index_split.out[0], ped_ch, slivar_selfchain)
    annotate(index_split.out[1], slivar_zip)
    //annotate.out.view()
    filtervcf(annotate.out, ped_ch, slivar_zip)
    SLIVAR_TSV(filtervcf.out[0], filtervcf.out[1], ped_ch)
}

// optional workflow for variant calling on trios, duos, or singletons
workflow calltrios {
    // Read in sample list
    Channel.fromPath(file(params.sample_bams, checkIfExists: true))
        .splitCsv(header: true, sep: ',')
        .map{
            row ->
                [ [proband_sex: row.proband_sex, proband_id: row.proband_id, father_id: row.father_id, mother_id: row.mother_id], // meta for family
                  [[id: row.proband_id, proband_id: row.proband_id, ord: 0],
                   [id: row.father_id, proband_id: row.proband_id, ord: 1],
                   [id: row.mother_id, proband_id: row.proband_id, ord: 2]], // meta for individuals
                  [row.proband_bam, row.father_bam, row.mother_bam],
                  [row.proband_index, row.father_index, row.mother_index]
                ]
            }
        .transpose()
        .filter{ it[1].id != "" & it[2] != "" }
        .map{ [it[0], it[1], file(it[2], checkIfExists: true), it[3]] }
        .set{ input_ch }

    // Get a channel just for looking up metadata again
    input_ch
        .map{ [it[1], it[0]] }
        .set{ meta_lookup }

    // Split off into samples input as bams vs fastqs
    input_ch
        .branch{ it ->
            fastq: it[2].name.endsWith(".fq.gz") | it[2].name.endsWith(".fastq.gz") | it[2].name.endsWith(".fq") | it[2].name.endsWith(".fastq")
            bam: true // get everything else
        }
        .set{ input_ch }

    // Align FASTQ files
    input_ch.fastq
        .map{ [it[1], it[2]] } // get metadata and fastq
        .set{ fastq_ch }
    fasta_bams
        .map{ [[id: it.simpleName], it]}
        .collect()
        .set{ fasta_for_bwa }
    Channel.fromPath(params.bwa_index)
        .map{ [[id: it.simpleName], it]}
        .collect()
        .set{ bwa_index }
    BWA_MEM(fastq_ch, bwa_index, fasta_for_bwa, true)

    // Combine new aligned bams with any existing bams
    BWA_MEM.out.bam
        .join(meta_lookup)
        .map{ [it[2], it[0], it[1], ""]}
        .concat(input_ch.bam)
        .set{ bam_ch }
    
    // Index any BAMS with missing index
    bam_ch
        .branch{ it ->
            needs_index: it[3] == ""
            indexed: true
        }
        .set{ bam_ch }
    bam_ch.needs_index
        .map{ [it[1], it[2]] }
        .set{ bams_to_index }
    SAMTOOLS_INDEX(bams_to_index)

    // Join back in with indexed BAMs
    bam_ch.indexed
        .map{ [it[0], it[1], it[2], file(it[3], checkIfExists: true)] }
        .set{ bams_indexed }
    bams_to_index
       .join(SAMTOOLS_INDEX.out.bai)
       .join(meta_lookup)
       .map{ [it[3], it[0], it[1], it[2]] }
       .concat(bams_indexed)
       .set{ bam_ch }

    // Get BAMS back into proband-father-mother order and split families from singletons
    bam_ch
        .map{ [ it[0], it[1].plus([bam: it[2], index: it[3]])] }
        .groupTuple(sort: { it.ord } )
        .map{ [ it[0], it[1].bam, it[1].index ] }
        .branch{ it ->
            single: it[0].father_id == "" & it[0].mother_id == ""
            family: true
        }
        .set{bam_ch}
    
    // Variant calling on families
    MAKE_EXAMPLES_TRIO(bam_ch.family, fasta_bams, fai_bams)
    MAKE_EXAMPLES_TRIO.out.proband_tfrecord
        .map{ [[proband_id: it[0].proband_id], it[0], it[1], it[2] ] }
        .join(MAKE_EXAMPLES_TRIO.out.example_info)
        .map{ [it[1], it[2], it[3], it[4]] }
        .set{ all_proband_tfrecords }
    MAKE_EXAMPLES_TRIO.out.father_tfrecord
        .map{ [[proband_id: it[0].proband_id], it[0], it[1], it[2] ] }
        .join(MAKE_EXAMPLES_TRIO.out.example_info)
        .map{ [it[1], it[2], it[3], it[4]] }
        .set{ all_father_tfrecords }
    MAKE_EXAMPLES_TRIO.out.mother_tfrecord
        .map{ [[proband_id: it[0].proband_id], it[0], it[1], it[2] ] }
        .join(MAKE_EXAMPLES_TRIO.out.example_info)
        .map{ [it[1], it[2], it[3], it[4]] }
        .set{ all_mother_tfrecords }
    all_proband_tfrecords
        .concat(all_father_tfrecords, all_mother_tfrecords)
        .filter{ it[0].id != "" }
        .set{ all_me_tfrecords }
    CALL_VARIANTS_TRIO(all_me_tfrecords)

    // Variant calling on singletons
    Channel.fromPath(file(params.par_bed, checkIfExists: true))
        .collect()
        .set{ par_bed }
    MAKE_EXAMPLES_SINGLE(bam_ch.single, fasta_bams, fai_bams, par_bed)
    MAKE_EXAMPLES_SINGLE.out.proband_tfrecord
        .join(MAKE_EXAMPLES_SINGLE.out.example_info)
        .set{ single_tfrecords }
    CALL_VARIANTS_SINGLE(single_tfrecords)

    // Join together families and singletons
    CALL_VARIANTS_TRIO.out
        .concat(CALL_VARIANTS_SINGLE.out)
        .set{ call_variants_out }

    // Postprocess both together
    POSTPROCESS_VARIANTS(call_variants_out, fasta_bams, fai_bams)
    POSTPROCESS_VARIANTS.out
        .map{ [it[1], it[2]] }
        .collect()
        .set{ all_ind_vcfs }

    // gl_nexus for joint VCF
    GLNEXUS(all_ind_vcfs, params.cohort_name) 
}
