Wrap-up R package to coordinate [SimusCop](https://github.com/qasimyu/simuscop) simulations. 

Also, here we present a easy conversion of VCF files in variations and snp files for SimuSCoP simuReads. This tools is useful if you want to simulate Illumina reads based on empirical variations already stored in a VCF file.

# Install simuscopR

Install from GitHub:

```{r}
library(devtools)

install_github("Cristianetaniguti/simuscopR")
```

*It also requires samtools installed.*

Or use the Docker image:

```{bash}
docker pull cristaniguti/simuscopr:latest
```

# Input files

All inputs must have information for just one chromosome and one sample.

* Reference genome and indexes from picard, bwa and samtools
* Sorted BAM file from real data
* VCF file
* BED file

If you need to filter the files to keep only one chromosome, we suggest:

```{bash, eval=FALSE}
docker run -v $(pwd):/opt cristaniguti/r-samtools samtools view -b /opt/your_bam_file Chr10 > filtered_bam

docker run -v $(pwd):/opt cristaniguti/vcftools vcftools --gzvcf /opt/vcf.file.gz --indv ID1 --chr Chr01 --recode --out vcf_filt
```

To obtain the BED file from the BAM file

```{bash}
docker run -it -v $(pwd):/data  biocontainers/bedtools:v2.28.0_cv1 bamtobed -i /data/PT_F.chr.sorted.bam > /data/PT_F.sorted_filt.bed
```

# Create profile

```{bash}
docker run -it -v $(pwd):/opt cristaniguti/simuscopr 

```

You can find these test files here `system.file("extdata", package="simuscopR")`.

```{r}
setwd("/opt")
seqToProfile("PT_F.chr.sorted.bam", "PT_F.sorted_filt.bed", "gatk_chr10_PT_F.recode.vcf",
             "reference/Chr10.populus.fa", "PT_F.profile")

```

# Variants files

```{r}
library(vcfR)

vcfR.object <- read.vcfR("gatk_chr10_PT_F.recode.vcf")

variants <- vcf2variants(vcfR.object, sample = "PT_F", chrom = "Chr10")

write.table(variants$SNVs, file = "SNVs.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(variants$indels, file = "indels.txt", sep = "\t", quote = F, col.names = F, row.names = F)
write.table(variants$insertions, file = "insertions.txt", sep = "\t", quote = F, col.names = F, row.names = F)

system("cat SNVs.txt indels.txt insertions.txt > variants.txt")
```

# Run simulation

```{r}
simuReads(ref = "reference/Chr10.populus.fa",
          profile = "PT_F.profile",
          variation = "variants.txt",
          target = "PT_F.sorted_filt.bed",
          name = "PT_F",
          output = ".",
          layout = "SE",
          threads = 6,
          verbose = 1,
          coverage = 15)
```


# Test to see if worked

```{bash}
docker run -v $(pwd):/data  kfdrc/bwa-picard:latest-dev  bin/bwa mem -t 2  \
/data/reference/Chr10.populus.fa /data/PT_F.fq  > file.bam

docker run -v $(pwd):/data  kfdrc/bwa-picard:latest-dev java -jar /picard.jar SortSam \
I="/data/file.bam" \
O="/data/PT_F.simulated.bam" \
TMP_DIR=./tmp \
SORT_ORDER=coordinate \
CREATE_INDEX=true

# Check images on IGV
```