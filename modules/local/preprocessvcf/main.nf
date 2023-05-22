nextflow.enable.dsl=2

include { LIFTOVER_HG19_TO_HG38 } from '../picard/liftovervcf/main'
include { CHROMNAMES_TO_UCSC } from '../python/chromtoucsc/main'
include { CREATE_SEQUENCE_DICTIONARY } from '../picard/createsequencedictionary/main'
include { PAD_VCF_INDELS } from '../python/pad_vcf_indels/main'
include { SAMTOOLS_FAIDX } from '../../nf-core/samtools/faidx/main'

// Channel setup for UCSC references for preprocessvcf
// UCSC Reference genome for Hg38
Channel.fromPath(file(params.fasta38ucsc, checkIfExists: true))
    .set{ fasta38ucsc }
// UCSC Reference genome for Hg19
Channel.fromPath(file(params.fasta19ucsc, checkIfExists: true))
    .set{ fasta19ucsc }
// Dictionary
fasta_basename = file(params.fasta38ucsc).baseName
if( file(params.fasta38ucsc).extension == "gz" ){
    fasta_basename = fasta_basename - ~/\.f(a|asta|na)$/
}
dict_file = file("${file(params.fasta38ucsc).parent}/${fasta_basename}.dict")
if( dict_file.exists() ){
    Channel.fromPath(dict_file)
        .set{ dict }
}
else {
    // If the dictionary doesn't exist already, make it but keep it in a
    // channel rather than writing to the reference genome directory.
    CREATE_SEQUENCE_DICTIONARY(fasta38ucsc)
    CREATE_SEQUENCE_DICTIONARY.out.dict
        .set{ dict }
}
// FAIDX for Hg19 for deletion padding
fai_file = file("${params.fasta19ucsc}.fai")
if( fai_file.exists() ){
    Channel.fromPath(fai_file)
        .set{ fai }
}
else {
    // Run Samtools Faidx
    SAMTOOLS_FAIDX(fasta19ucsc)
    SAMTOOLS_FAIDX.out.fai
        .set{ fai }
}

workflow preprocessvcf {
    take:
        vcf_ch
    //    fasta19ucsc
    //    fasta38ucsc
    //    dict
    //    fai
    emit:
        vcf = vcf_ch
    main:
        if( params.do_liftover ){
            // Get chromosome names into UCSC format if needed for liftover
            // Possibly insert something here to run referencedict workflow if needed
            if( params.chromnames != "ucsc" ){
                CHROMNAMES_TO_UCSC(vcf_ch)
                CHROMNAMES_TO_UCSC.out.vcf
                    .set{ vcf_ch }
            }
            // Make sure there are no malformed deletions (blank ALT alleles)
            PAD_VCF_INDELS(vcf_ch, fasta19ucsc, fai)
            PAD_VCF_INDELS.out.vcf
                .set{ vcf_ch }
            // Convert from Hg19 to Hg38
            LIFTOVER_HG19_TO_HG38(vcf_ch, fasta38ucsc, dict)
            LIFTOVER_HG19_TO_HG38.out.vcf
                .set{ vcf_ch }
        }
}

workflow referencedict {
    // Use Picard to create a sequence dictionary needed for liftover.
    //take:
    //    fasta38
    emit:
        dict = CREATE_SEQUENCE_DICTIONARY.out.dict
    main:
        CREATE_SEQUENCE_DICTIONARY(fasta38ucsc)
        // Write dictionary to same directory as reference genome.
        CREATE_SEQUENCE_DICTIONARY.out.dict
            .subscribe{file("$it").copyTo("${file(params.fasta38ucsc).parent}/")}
}
