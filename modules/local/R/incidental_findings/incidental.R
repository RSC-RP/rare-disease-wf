# Script to process a VCF and output a short report for each subject with one
# or more incidental findings relating to known genetic disorders.

args <- R.utils::commandArgs(asValues = TRUE,
                             defaults = list(clinvar_only = TRUE))
# args <- list(invcf = "~/heike_c_2022.05_cause_wgs/results/luquetti_grc_combined_2.filtered_pathogenic.vcf",
#              txdb = "TxDb.Hsapiens.UCSC.hg19.knownGene",
#              build = "hg19",
#              sample_list = "/active/taylor_s/people/lclar5/CAUSE/heike_c_2022.05_cause_wgs/samples/CAUSE_samples.txt",
#              clinvar_only = TRUE)

library(VariantAnnotation)
library(args$txdb, character.only = TRUE)

# Read in list of genes ####
# From https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/

af <- "/opt/incidental_findings/acmg_genes.txt" # location in Docker container
# af <- "acmg_genes.txt"
acmg_genes <- read.table(af, header = TRUE, sep = "\t",
                         colClasses = "character")

txdb_genes <- genes(get(args$txdb))

acmg_genes_gr <- txdb_genes[mcols(txdb_genes)$gene_id %in% acmg_genes$ENTREZID]

temp <- match(acmg_genes_gr$gene_id, acmg_genes$ENTREZID)

acmg_genes_gr$SYMBOL <- acmg_genes$SYMBOL[temp]
acmg_genes_gr$Disease <- acmg_genes$Disease.name.and.MIM.number[temp]
acmg_genes_gr$MIM_gene <- acmg_genes$MIM_gene[temp]

seqlevels(acmg_genes_gr) <- seqlevelsInUse(acmg_genes_gr)

# Read in VCF ####
#bg <- paste0(args$invcf, ".bgz")
bg <- bgzip(args$invcf)
vf <- indexVcf(bg)

# switch chromosome name format if necessary
if(!any(startsWith(seqnamesTabix(vf), "chr"))){
    seqlevels(acmg_genes_gr) <- sub("^chr", "", seqlevels(acmg_genes_gr))
}

# subset genes based on chromosome names present in VCF
acmg_genes_gr <- acmg_genes_gr[seqnames(acmg_genes_gr) %in% seqnamesTabix(vf)]

my_variants <- readVcf(vf, param = ScanVcfParam(which = acmg_genes_gr))

# Filter to rare variants
if("gnomad_popmax_af" %in% colnames(info(my_variants))){
  my_variants <- my_variants[info(my_variants)$gnomad_popmax_af < 0.03,]
}

# Filter to Pathogenic and Likely pathogenic in CLNSIG
if(args$clinvar_only){
  clnsig_strings <-
    sapply(info(my_variants)$CLNSIG,
           function(x) paste(x, collapse = " "))
  my_variants <-
    my_variants[grepl("(P|Likely_p)athogenic", clnsig_strings),]
}

# Read in list of samples to process ####
my_samples <- readLines(args$sample_list)

# Identify variants detected in samples of interest ####
GT <- geno(my_variants)$GT[,my_samples, drop = FALSE]
al_present <- matrix(grepl("1", GT), nrow = nrow(GT), ncol = ncol(GT),
                     dimnames = dimnames(GT))
my_variants_present <- my_variants[rowSums(al_present) > 0,]

if(nrow(my_variants_present) == 0){
  writeLines("None found",
             con = sub("\\..*", ".incidental.txt", args$invcf))
  q(save = "no", status = 0)
}

# Table of variant info ####
vardf <- data.frame(chromosome = seqnames(my_variants_present),
                    position = start(my_variants_present),
                    ref = as.character(rowRanges(my_variants_present)$REF),
                    alt = as.character(unlist(rowRanges(my_variants_present)$ALT)))
vardf$build <- args$build
acmg_rows <- sapply(seq_len(nrow(my_variants_present)),
                    function(i){
                        which(acmg_genes$SYMBOL %in%
                                  info(my_variants_present)$Gene.refGene[[i]])[1]
                    })
vardf <- cbind(vardf, acmg_genes[acmg_rows,])
vardf$ClinVar_alleleID <- unlist(info(my_variants_present)$ALLELEID)
vardf$ClinVar_significance <- unlist(info(my_variants_present)$CLNSIG)

# Build table of information about each sample*variant combination ####
sv_df <- data.frame(var_index = integer(0),
                    Participant = character(0),
                    GT = character(0),
                    AD = character(0),
                    DP = integer(0),
                    GQ = integer(0))

# Loop through variants and add samples
al_present2 <- al_present[rownames(my_variants_present),, drop = FALSE]
for(i in seq_len(nrow(vardf))){
    these_samples <- colnames(al_present2)[al_present2[i,]]
    for(sam in these_samples){
        newrow <- list(var_index = i, Participant = sam,
                       GT = geno(my_variants_present)$GT[i,sam],
                       AD = paste(geno(my_variants_present)$AD[[i,sam]], collapse = ","),
                       DP = geno(my_variants_present)$DP[i,sam],
                       GQ = geno(my_variants_present)$GQ[i,sam])
        sv_df <- rbind(sv_df, newrow)
    }
}

out <- cbind(vardf[sv_df$var_index,], sv_df[,-1])

write.table(out, sep = "\t", row.names = FALSE, col.names = TRUE,
            file = sub("\\..*", ".incidental.txt", args$invcf))
