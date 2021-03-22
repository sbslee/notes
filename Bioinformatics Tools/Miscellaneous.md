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




```
Link to GTCtoVCF (Illumina’s):
	https://github.com/Illumina/GTCtoVCF

Link to gtc2vcf:
  	https://github.com/freeseek/gtc2vcf
```


```
Getting exon coordinates for a gene
grep -w "CYP2A6" Homo_sapiens.GRCh37.75.gtf | grep "CYP2A6-001" | grep -w "exon" | cut -f 1,4,5,9 -d$'\t' | cut -f 1,3 -d';' | sed 's/gene_id "ENSG00000255974"; //g'

# How to get rpkm values from a .gct file:
for sample in `cat liverbank_rnaseq_id.list`; do
  printf "`echo $sample`\t`grep -w "CYP2A7" /net/grc/vol6/data/processed/samples/$sample/RNA_SEQ/qc/genes.rpkm.gct`\n";
done > liverbank_cyp2a7_rpkm.txt
```
