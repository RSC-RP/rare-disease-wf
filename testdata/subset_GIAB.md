## Subset GIAB Ashkenazim trio BAMs for example dataset

Starting from the large, public dataset, I want to make much smaller BAMs to
enable a quick test of the pipeline.

Regions to capture:

* chr2:179300000-179600000 : Comphet mutations in TTN found
* chrX:122530000-122550000 : X-linked recessive for mental retardation, inherited from mother
* chr16:450000-470000 : frameshift in DECR2, de-novo

I put these regions into the file subset_GIAB.bed.

Subset using SAMtools

``` bash
qsub -I -q sceaq -l select=1:ncpus=1:mem=16g -l walltime=2:00:00 -P 207f23bf-acb6-4835-8bfe-142436acb58c

cd /active/taylor_s/people/lclar5/RPDEV/rare-disease-wf

module load SAMtools

samtools view -H /gpfs/shared_data/GIAB/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam

samtools view -b -M -L testdata/subset_GIAB.bed \
/gpfs/shared_data/GIAB/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam  > \
testdata/bams/HG002_subset.bam

samtools view -b -M -L testdata/subset_GIAB.bed \
/gpfs/shared_data/GIAB/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bam  > \
testdata/bams/HG003_subset.bam

samtools view -b -M -L testdata/subset_GIAB.bed \
/gpfs/shared_data/GIAB/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup.bam  > \
testdata/bams/HG004_subset.bam

samtools index testdata/bams/HG002_subset.bam
samtools index testdata/bams/HG003_subset.bam
samtools index testdata/bams/HG004_subset.bam
```

Copied these to /gpfs/shared_data/GIAB/OsloUniversityHospital_Exome/ using the service account.
