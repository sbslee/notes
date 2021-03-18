# Miscellaneous

## Commonly used commands

To zip a VCF file:

```
bgzip -c sample.vcf > sample.vcf.gz
```

To index a VCF file:

```
tabix -p vcf sample.vcf.gz
```

To slice a VCF file:

```
tabix -h sample.vcf.gz chr1:1000-2000 > sliced.vcf.gz
```

To count the number of sequence reads in a FASTQ file:

```
echo $(cat sample.fastq | wc -l) / 4 | bc
```

To count the number of sequence reads in a zipped FASTQ file:

```
echo $(zcat sample.fastq | wc -l) / 4 | bc
```

To count the number of sequence reads in a zipped FASTQ file (macOS):

```
echo $(zcat < sample.fastq | wc -l) / 4 | bc
```
