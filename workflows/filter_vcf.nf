include { BCFTOOLS_NORM_CSQ } from '../modules/local/bcftools/norm_csq/main'
include { SNPEFF } from '../modules/nf-core/snpeff/main'
include { SLIVAR_GNOTATE } from '../modules/local/slivar/gnotate/main'
include { TABLE_ANNOVAR } from '../modules/local/annovar/main'
include { MVP_ANNO } from '../modules/local/python/mvp/main'
include { ZIP_INDEX_PUBLISH } from '../modules/local/bcftools/zip_index_publish/main'
include { SPLIT } from '../modules/local/raredisease/split/main'
include { BCFTOOLS_CONCAT } from '../modules/local/bcftools/concat/main'
include { CHROMNAMES } from '../modules/local/python/chromnames/main'
include { FILEORDER } from '../modules/local/python/fileorder/main'
include { SLIVAR_MENDEL_FILTER } from '../modules/local/slivar/mendel_filter/main'
include { CLINVAR_DOWNLOAD } from '../modules/local/clinvar_download/main'
include { SNPSIFT_ANNOTATE } from '../modules/local/snpsift/annotate/main'
include { FILTER_PATHOGENIC } from '../modules/local/python/filter_pathogenic/main'
include { POP_STRUCT } from '../modules/local/R/pop_struct/main'
include { SLIVAR_TSV } from '../modules/local/slivar/tsv/main'
include { SLIVAR_DUODEL } from '../modules/local/slivar/duodel/main'
include { SAMTOOLS_FAIDX } from '../modules/nf-core/samtools/faidx/main'
include { INCIDENTAL_FINDINGS } from '../modules/local/R/incidental_findings/main'

// workflow to sort and index vcf and split it up by chromosome
workflow index_split{
    take:
        vcf_ch
    main:
        // Reference sequence for alignment and genotype calling, which may be different from Ensembl.
        fai_file2 = file("${params.fasta_bams}.fai")
        if( fai_file2.exists() ){
            Channel.fromPath(fai_file2)
                .collect()
                .set{ fai_bams }
        }
        else {
            Channel.fromPath(file(params.fasta_bams, checkIfExists: true))
                .collect()
                .set{ fasta_bams }
            // Run Samtools Faidx
            SAMTOOLS_FAIDX(fasta_bams)
            SAMTOOLS_FAIDX.out.fai
                .collect()
                .set{ fai_bams }
        }

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
    emit:
        vcf_ch
        vcf_ch_sp
}

// workflow for QC
workflow qc {
    take:
        vcf_ch // tuple of VCF and index
        metadata_ch
    main:
        if(params.run_popstruct){
            POP_STRUCT(vcf_ch,
                "TxDb.Hsapiens.UCSC.${params.annovar_buildver}.knownGene",
                metadata_ch,
                params.metadata_columns,
                params.umap_seed,
                params.umap_min_call_rate,
                params.umap_n_exons,
                params.umap_n_pcs)
        }
}

// workflow to merge VCFs into one file
workflow merge_vcfs {
    take:
        vcf_ch_sp
        fai
    main:
        CHROMNAMES(fai)
        FILEORDER(vcf_ch_sp, CHROMNAMES.out)
        BCFTOOLS_CONCAT(vcf_ch_sp, FILEORDER.out, file(params.vcf).simpleName)
            .set{ vcf_ch }
    emit:
        vcf_ch
}

// workflow for annotation
workflow annotate {
    take:
        vcf_ch_sp
        slivar_zip
        fai
    main:
        // Convert from list to channel emitting individual items
        vcf_ch_sp
            .flatten()
            .set{ vcf_ch_sp }
        
        // Channel setup for Ensembl references for BcfTools and Slivar
        Channel.fromPath(file(params.fasta_ensembl, checkIfExists: true))
            .collect() // allows the reference to be used with multiple input VCFs
            .set{ fasta_ensembl }
        Channel.fromPath(file(params.gff_ensembl, checkIfExists: true))
            .collect()
            .set{ gff_ensembl }
        
        // Chromsome conversion file
        Channel.fromPath(file(params.chrom_convert, checkIfExists: true))
            .collect()
            .set{ chrom_convert }

        // Annotation steps
        BCFTOOLS_NORM_CSQ(vcf_ch_sp, fasta_ensembl, fai, gff_ensembl, chrom_convert)
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
            merge_vcfs(vcf_ch_sp, fai)
            // Export annotated VCF before filtering
            ZIP_INDEX_PUBLISH(merge_vcfs.out, 'annotated', 'view')
        }
    emit:
        vcf_ch_sp
}

// Short workflows for publishing after filtering
workflow publish_mendelian {
    take:
        vcf_ch_sp
        fai

    main:
    // Concatenate individual VCFs into one
    merge_vcfs(vcf_ch_sp, fai)
    // Export annotated VCF after filtering
    ZIP_INDEX_PUBLISH(merge_vcfs.out, 'mendelian', 'view')
    ZIP_INDEX_PUBLISH.out
        .set{ vcf_ch }

    emit: vcf_ch
}

workflow publish_comphet {
    take:
        vcf_ch_sp
        fai

    main:
    if(params.trios_present){
        // Concatenate individual VCFs into one
        merge_vcfs(vcf_ch_sp, fai)
        // Export annotated VCF after filtering
        ZIP_INDEX_PUBLISH(merge_vcfs.out, 'comphet', 'view')
        ZIP_INDEX_PUBLISH.out
            .set{ vcf_ch }
    } else {
        Channel.fromPath('dummy.vcf.gz')
            .set{ vcf_ch }
    }

    emit: vcf_ch
}

// subworkflow for incidental findings
workflow incidental {
    take:
        vcf_ch_sp
        samples_ch
        fai
    main:
        merge_vcfs(vcf_ch_sp, fai)
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
        fai
    main:
        // Convert from list to channel emitting individual items
        vcf_ch_sp
            .flatten()
            .set{ vcf_ch_sp }
        
        // Filtering steps and incidental findings
        FILTER_PATHOGENIC(vcf_ch_sp)
        if(params.do_incidental){
            Channel.fromPath(file(params.incidental_samples, checkIfExists: true))
                .set{ samples_ch }
            incidental(FILTER_PATHOGENIC.out.collect(), samples_ch, fai)
        }
        Channel.fromPath(file(params.slivar_fn_loose, checkIfExists: true))
            .collect()
            .set{ slivar_fn_loose }
        SLIVAR_MENDEL_FILTER(FILTER_PATHOGENIC.out, ped_ch, slivar_zip,
            params.maf_recessive, params.maf_dominant, params.maf_denovo, params.comphet_nhomalt,
            slivar_fn_loose)
            .set{ vcf_ch_sp }
        
        // Convert back from individual items to list
        vcf_ch_sp.mendelian_vcf
            .collect()
            .set{ vcf_ch_sp_mendelian }
        vcf_ch_sp.comphet_vcf
            .collect()
            .set{ vcf_ch_sp_comphet }

        // Export
        publish_mendelian(vcf_ch_sp_mendelian, fai)
        publish_mendelian.out
            .set{ vcf_ch_mendelian }
        publish_comphet(vcf_ch_sp_comphet, fai)
        publish_comphet.out
            .set{ vcf_ch_comphet }
    emit:
        vcf_ch_mendelian
        vcf_ch_comphet
}

// main workflow
workflow FILTER_VCF {
    // Slivar file with gnomAD data
    Channel.fromPath(file(params.slivar_zip, checkIfExists: true))
        .collect()
        .set{ slivar_zip }
    // Slivar selfchain file for duo-del
    Channel.fromPath(file(params.slivar_selfchain, checkIfExists: true))
        .collect()
        .set{ slivar_selfchain }

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
    
    // FAIDX for reference genome
    fai_file = file("${params.fasta_ensembl}.fai")
    if( fai_file.exists() ){
        Channel.fromPath(fai_file)
            .collect()
            .set{ fai_ensembl }
    }
    else {
        Channel.fromPath(file(params.fasta_ensembl, checkIfExists: true))
            .collect()
            .set{ fasta_ensembl }
        // Run Samtools Faidx
        SAMTOOLS_FAIDX(fasta_ensembl)
        SAMTOOLS_FAIDX.out.fai
            .collect()
            .set{ fai_ensembl }
    }

    //preprocessvcf(vcf_ch)
    //preprocessvcf.out.view()
    index_split(vcf_ch)
    //index_split.out.view()
    qc(index_split.out[0], metadata_ch)
    SLIVAR_DUODEL(index_split.out[0], ped_ch, slivar_selfchain)
    annotate(index_split.out[1], slivar_zip, fai_ensembl)
    //annotate.out.view()
    filtervcf(annotate.out, ped_ch, slivar_zip, fai_ensembl)
    SLIVAR_TSV(filtervcf.out[0], filtervcf.out[1], ped_ch)
}
