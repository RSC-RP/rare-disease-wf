FROM bioconductor/bioconductor_docker:RELEASE_3_16

LABEL \
    author="Lindsay Clark" \
    description="Bioconductor image for visualizing population structure from a VCF" \
    maintainer="Lindsay.Clark@seattlechildrens.org"

ADD ./lentschroot.cer.crt /usr/local/share/ca-certificates/lentschroot.crt
ADD ./lpki30v2.cer.crt /usr/local/share/ca-certificates/lpki30.crt
ADD ./lschintrootca2.cer.crt /usr/local/share/ca-certificates/lschintrootca2.crt
ADD ./lpki10.cer.crt /usr/local/share/ca-certificates/lpki10.crt
ADD ./lzscalerInt.cer.crt /usr/local/share/ca-certificates/lzscalerInt.crt
ADD ./lentschissuing.cer.crt /usr/local/share/ca-certificates/lentschissuing.crt
ADD ./lpki20.cer.crt /usr/local/share/ca-certificates/lpki20.crt
ADD ./lschintrootca.cer.crt /usr/local/share/ca-certificates/lschintrootca.crt
ADD ./lzscalerroot.cer.crt /usr/local/share/ca-certificates/lzscalerroot.crt

RUN update-ca-certificates

RUN R -e 'BiocManager::install(c("VariantAnnotation", "pcaMethods", "umap", "ggplot2", "dittoSeq"))'
RUN R -e 'BiocManager::install(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "TxDb.Hsapiens.UCSC.hg38.knownGene"))'
RUN R -e 'BiocManager::install(c("R.utils"))'

COPY run_umap.R /
