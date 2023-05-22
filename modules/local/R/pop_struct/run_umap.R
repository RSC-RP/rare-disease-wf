#!/usr/bin/env Rscript

# Script to extract some random variants from a VCF and run UMAP to get a
# visualization of population structure.

args <- R.utils::commandArgs(asValues = TRUE,
                             defaults = c(min_call_rate = "0.95", n_exons = "5000", n_pcs = "50"))


if(!require(args$txdb, quietly = TRUE, character.only = TRUE)){
    BiocManager::install(args$txdb)
}

library(VariantAnnotation)
library(args$txdb, character.only = TRUE)
library(pcaMethods)
library(umap)
library(ggplot2)

set.seed(as.integer(args$seed))

# Error check VCF against sample metadata
hdr <- scanVcfHeader(args$invcf)
metadata <- read.table(args$sample_metadata, header = TRUE, sep = "\t")
if(!all(samples(hdr) %in% metadata[[1]])){
    stop("Sample names must match between VCF and sample metadata")
}

# Setup to read VCF
ind1 <- paste0(args$invcf, ".tbi")
ind2 <- paste0(args$invcf, ".csi")
vcf_file <- VcfFile(args$invcf, index = ifelse(file.exists(ind1), ind1, ind2))
svp <- ScanVcfParam(geno = "GT")

# If file is big, subset by random exons, otherwise read whole thing.
if(file.info(args$invcf)$size > 1e9){
    # Get some random exons and match the chromosome names with the VCF
    ex <- exons(get(args$txdb))
    
    random_exons <- ex[sort(sample(length(ex), as.integer(args$n_exons), replace = FALSE))]
    
    chrnames_vcf <- seqnames(seqinfo(hdr))
    chrnames_ex <- seqnames(seqinfo(random_exons))
    
    chrnames_ex_new <- chrnames_ex
    not_found <- !chrnames_ex %in% chrnames_vcf
    chrnames_ex_new[not_found] <- sub("chr", "", chrnames_ex_new[not_found])
    chrnames_ex_new[!chrnames_ex_new %in% chrnames_vcf] <- NA
    
    random_exons <- random_exons[seqnames(random_exons) %in% chrnames_ex[!is.na(chrnames_ex_new)]]
    seqlevels(random_exons) <- intersect(seqlevels(random_exons), chrnames_ex[!is.na(chrnames_ex_new)])
    seqlevels(random_exons) <- chrnames_ex_new[!is.na(chrnames_ex_new)]
    
    vcfWhich(svp) <- random_exons
}

# Import data from VCF
random_vcf <- readVcf(vcf_file, param = svp)

# Convert to numeric
biallelic <- apply(geno(random_vcf)$GT, 1,
                   function(x){
                       all(x %in% c("./.", ".|.", "0/0", "0/1", "1/1", "0|0", "0|1", "1|1"))
                   })

temp <- geno(random_vcf)$GT[biallelic,]
numgeno <- matrix(0L, nrow = nrow(temp), ncol = ncol(temp),
                  dimnames = dimnames(temp))
numgeno[temp %in% c("./.", ".|.")] <- NA_integer_
numgeno[temp %in% c("0/1", "0|1")] <- 1L
numgeno[temp %in% c("1/1", "1|1")] <- 2L

# Filter by call rate
numgeno2 <- numgeno[rowMeans(!is.na(numgeno)) >= as.numeric(args$min_call_rate),]

# PCA and UMAP
mypca <- pca(t(numgeno2), method = "ppca", scale = "uv",
             nPcs = min(c(as.integer(args$n_pcs), ncol(numgeno2) %/% 2L)))
umap.custom <- umap.defaults
umap.custom$n_neighbors <- min(c(umap.defaults$n_neighbors, ncol(numgeno2) - 1L))
myumap <- umap(scores(mypca), method = "naive", config = umap.custom)

# Put coordinates in table with sample metadata
temp <- match(metadata[[1]], rownames(myumap$layout))
metadata$UMAP_1 <- myumap$layout[temp,1]
metadata$UMAP_2 <- myumap$layout[temp,2]
metadata <- cbind(metadata, scores(mypca)[temp,])

# Export table of coordinates
write.table(metadata, paste0("pop_struct_ordination_coordinates_", Sys.Date(), ".tsv"),
            col.names = TRUE, sep = "\t", row.names = FALSE)

# Make plots
meta_cols <- strsplit(args$metadata_columns, ",")[[1]]

for(var in meta_cols){
    if(is.numeric(metadata[[var]]) & length(unique(metadata[[var]])) < 11){
        metadata[[var]] <- as.factor(metadata[[var]])
    }

    p <- ggplot(metadata, aes(x = UMAP_1, y = UMAP_2)) +
        geom_point(aes(color = .data[[var]]))
    if(is.numeric(metadata[[var]])){
        p <- p + scale_color_viridis_c()
    } else {
        p <- p + scale_color_manual(values = dittoSeq::dittoColors())
    }
    tiff(paste0("UMAP_", var, "_", Sys.Date(), ".tiff"),
         res = 200, width = 1800, height = 1400, compression = "lzw")
    print(p)
    dev.off()
}
