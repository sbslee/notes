# BCFtools

# Variant calling pipeline

1. Calculate genotype likelihoods at each genomic position with coverage.

```
bcftools mpileup -Ou -q 1 -a AD --max-depth 1000 -f ref.fa -r chr1:1000-2000 -o sample.bcf sample.bam
```

2. Make variant calls.

```
bcftools call -Oz -mv -o sample.vcf.gz sample.bcf
```

3. Index the VCF file.

```
bcftools index sample.vcf.gz
```

4. Left-align and normalize indels.

```
bcftools norm -Ob -f ref.fa -o sample.normalized.bcf sample.vcf.gz
```

5. Filter variants.

```
bcftools filter -Ov --IndelGap 5 -o sample.normalized.filtered.vcf sample.normalized.bcf
```
