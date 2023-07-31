## Filtering a VCF to find rare disease causing variants

This workflow is heavily modified from https://github.com/brentp/rare-disease-wf.
Files from the original workflow can be found in the `OLD` folder.

Special thanks also to https://github.com/dpryan79/ChromosomeMappings for
chromsome name conversion keys.

For the main workflow,
a VCF from WGS and/or WES is required, along with a pedigree file indicating
relationships between probands and parents or other relatives.  Functional
annotation is performed with a variety of tools.  The VCF is then filtered
just to variants that are annotated as potentially pathogenic, and then filtered
again based on segregation patterns in families.

Within `main.nf` there is also a workflow called `calltrios ` that can be used
to perform variant calling, starting from BAMs from trios and duos. DeepTrio
is used for variant calling and outputs one GVCF per sample.  GLNexus is then
used to generate the joint called VCF.

The `modules/local` folder contains documented custom modules for this workflow,
which I hope can be reused by our group.

Although I have not integrated it into this workflow, I recommend running
Somalier (https://github.com/brentp/somalier) on the VCF (or, if available,
BAMs) and pedigree file before running the main workflow.  I have found it to be
very helpful for identifying sample mixups, incorrect paternity, incorrect sex,
and Klinefelter syndrome. See `workbooks/somalier_template.Rmd` for example
code.

Older commits of this repository are available to RSC team members via
https://childrens-atlassian/bitbucket/projects/RPDEV/repos/rare-disease-wf/browse.
The full history is not available in the public version since it contains data
that I do not have permission to distribute.

## How to use the workflow

Make your own fork on Bitbucket (recommended), or simply clone this repo using

``` bash
git clone ssh://git@childrens-atlassian:7999/rp/rare-disease-wf.git
```

**To run your own data:** You should open `nextflow.config` and edit it to
point to your own files. The params that you are likely to need to change are
closest to the top, starting from `project` and going down to `outdir`.

* You can use `//` to indicate a comment line in `nextflow.config` if you don't
want to delete the original parameters.
* If you are running the variant calling workflow on your own data (i.e. you are
starting from BAMs rather than from a VCF), be sure to set `test_bams = false`
and `make_examples_nshards = 32`.  The original settings (`test_bams = true` and
`make_examples_nshards = 4`) are intended for the tiny example dataset.
* The `nextflow.config` file is set up for hg19/GRCh37. If you have aligned to
hg38/GRCh38, open the file `hg38.config` and copy the six parameters from it over
to `nextflow.config`, overwriting or commenting out the originals.
* If your chromosome names start with `chr`, you should set `chromnames = 'ucsc'`.
You should also set `fasta_bams` to the reference genome that was used for alignment,
for example
`/gpfs/shared_data/references/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/hg19.fa` or
`/gpfs/shared_data/references/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/hg38.p13.fa`.
You can check the file header of your BAMs or VCF to be sure.
This needs to be done regardless of whether you are running the variant calling
pipeline or just the trio analysis (variant filtering) pipeline.
* You can look at the datasets in `testdata/` to see how to format your files.
Note that although the example dataset only contains one family, yours can
contain many. The example dataset is a VCF, but it is fine to use a gzipped
VCF or a BCF.

**To run an example dataset:** This workflow comes with tiny example BAMs and a
VCF from the Genome in a Bottle (GIAB) Ashkenazim trio. To run the pipeline on
this dataset, the only thing you have to change in `nextflow.config` is the
`project` parameter.

When you are ready to run the workflow, you may want to start a `screen` or `tmux`
session so that you can close the window if needed.  Use `qsub` to start an interactive
session, and give it at least 8 hours or so of walltime.

Set up the conda environment for Nextflow:

``` bash
cd rare-disease-wf

mamba env create -f env/nextflow.yaml
```

Load the conda environment.

``` bash
mamba activate nextflow
```

The first time you run the workflow, you will need to load apptainer (previously called "singularity") in order to download the container images.

``` bash
module load apptainer
```

You should also copy over the custom apptainer images for this workflow before running the workflow for the first time:

``` bash
mkdir -p ~/apptainer
cp /gpfs/shared_data/singularity/rare-disease-wf/* ~/apptainer/
```

This is how you would run variant calling on trios and duos:

``` bash
nextflow \
    -c nextflow.config \
    run main.nf \
    -entry calltrios \
    -profile PBS_apptainer \
    -resume
```

And this is how you would run trio analysis:

``` bash
nextflow \
    -c nextflow.config \
    run main.nf \
    -profile PBS_apptainer \
    -resume
```

## Outputs

In the directory specified by `outdir` in your `nextflow.config` (the `results`
directory of this repo by default), output files will be copied over from your
Nextflow work directory (`~/temp` by default).

For the `calltrios` workflow, you will get a subfolder called `gvcfs` containing
one gVCF per sample, generated by Deep Trio. gVCFs will be named according to
sample ID from the CSV passed to the `sample_bams` parameter. In the main output
directory you will get a file ending in `.glnexus.bcf`, named according to
`cohort_name`, containing the joint genotype calls from GLnexus. This BCF is
suitable to be passed directly to the main workflow (i.e. you can edit the `vcf`
parameter to point to this file after running variant calling).

From the main workflow, you will get several indexed, gzipped VCFs:

* One ending in `.original.vcf.gz`, containing all variants from the original VCF.
* If you set `export_prefilt = true`, you will get a VCF ending in `.annotated.vcf.gz`
containing all of the original variants, with `INFO` columns added containing the
annotations from various tools.
* One ending in `.mendelian.vcf.gz` containing filtered variants tagged by Slivar
as de novo, dominant, recessive, X de novo, and/or X recessive.
* One ending in `.comphet.vcf.gz` containing filtered variants tagged by Slivar as
complementary heterozygotes.

You will also get a file named according to your VCF and ending in `.slivar.tsv`
containing summarized results on variants of interest, output by Slivar. If you
set `do_incidental = true`, there will be a folder called `incidental` containing
a tab-delimited text file named after your VCF, listing pathogenic/likely pathogenic
variants that are recommended by the American College of Medical Genetics to be
reported to participants if they consented to receive such information. There will
also be a folder called `qc` containing population structure results on the dataset,
including TIFF files of UMAP plots colored by demographic data, and a tab-delimited
text file listing coordinates for all samples.

## Methods

Below is some suggested text for methods writeups. Edit as necessary.
If line breaks are a problem for copying and pasting into Word, just copy from
the rendered README on Bitbucket.

### Genotype calling

Variant calling was performed on BAM files using DeepTrio v1.5.0
(Kolesnikov et al. 2021, doi:10.1101/2021.04.05.438434), which uses
parent-offspring relationships to improve calling accuracy.  GVCFs generated
by DeepTrio were then passed to GLNexus v1.4.1
(Lin et al. 2018, doi:10.1101/343970) for joint genotype calling.

### Quality control

To visualize population structure, variants within 5000 random exons were
imported into R, then filtered down to biallelic markers with a minumum call
rate of 95%. Principal components analysis was performed using the probabalistic
PCA (PPCA) method in the pcaMethods Bioconductor package
(Stacklies et al. 2007, doi:10.1093/bioinformatics/btm069).  The first 50 axes
were then subjected to UMAP (McInnes et al. 2020, doi:10.48550/arXiv.1802.03426)
to generate a two-dimensional visualization of participants. Points were colored
by race, ethnicity, study site, [other variables] in order to confirm that
clustering was primarily by ancestry and that no batch effects were present.

### Annotation and filtering

The Nextflow workflow available at https://github.com/brentp/rare-disease-wf was
modified to take a joint-called VCF as input and to perform additional annotation
and filtering. BCFtools (Li et al. 2009, doi:10.1093/bioinformatics/btp352),
SnpEff (Cingolani et al. 2012, doi:10.4161/fly.19695),
Slivar (Pedersen et al. 2021, doi:10.1038/s41525-021-00227-3),
ANNOVAR (Wang et al. 2010, doi:10.1093/nar/gkq603),
and SnpSift (Cingolani et al. 2012, doi:10.3389/fgene.2012.00035) were used to add
functional annotation, predictions of pathogenicity,
gnomAD allele frequencies (Karczewski et al. 2020, doi:10.1038/s41586-020-2308-7),
and data from ClinVar (Landrum et al. 2014, doi:10.1093/nar/gkt1113) to the VCF.
The VCF was then filtered to retain variants that met any of the following
criteria: (1) predicted as deleterious by
MetaSVM (Dong et al. 2015, doi:10.1093/hmg/ddu733),
FATHMM (Shihab et al. 2013, doi:10.1002/humu.22225),
Polyphen2 HDIV or HVAR (Adzhubei et al. 2013, doi:10.1002/0471142905.hg0720s76),
or fathmm-MKL_coding (Shihab et al. 2015, doi:10.1093/bioinformatics/btv009);
(2) predicted splicing or frameshift variant by ANNOVAR;
(3) MVP (Qi et al 2021, doi:10.1038/s41467-020-20847-0) score greater than 0.9,
CADD (Kircher et al. 2014, doi:10.1038/ng.2892) score greater than 7,
DANN (Quant et al. 2015, doi:10.1093/bioinformatics/btu703) score of at least 0.999,
or GERP++ (Davydov et al. 2010, doi:10.1371/journal.pcbi.1001025) score greater than 5.5; or
(4) listed as pathogenic or likely pathogenic in ClinVar.
Slivar was then used to identify potential causative variants based on inheritance
in families, with a gnomAD minor allele frequency cutoff of 0.001 for recessive
or de novo variants and 0.00002 for dominant variants, or fewer than 10
homozygous individuals in gnomAD for complementary heterozygous variants.

### Incidental findings

Across XXX participants who had consented to receive incidental findings, the
VCF containing genotype calls was filtered to variants that were within any of
73 genes listed by the American College of Medical Genetics and Genomics’
recommendations for reporting incidental findings
(ACMB Board of Directors 2015, doi: 10.1038/gim.2014.151; https://www.ncbi.nlm.nih.gov/clinvar/docs/acmg/).
These variants were further filtered to those that were listed as Pathogenic or
Likely Pathogenic by ClinVar, yielding a total of XX variants across XX
participants, [all of which were found in a heterozygous state]
[with X being found in a homozygous state and the rest in a heterozygous state].
Variants were then reviewed individually for genotyping quality, as well as to
determine whether they were likely to be dominant or recessive.
