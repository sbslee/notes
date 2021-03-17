# BCFtools

# Variant calling pipeline

```
bcftools mpileup -Ou -q 1 -a AD --max-depth 1000 -f ref.fa -r chr1:1000-2000 -o sample.bcf sample.bam
```

```
bcftools call -Oz -mv -o sample.vcf.gz sample.bcf
```

```
bcftools index sample.vcf.gz
```

```
bcftools norm -Ob -f ref.fa -o sample.normalized.bcf sample.vcf.gz
```

```
bcftools filter -Ov --IndelGap 5 -o sample.normalized.filtered.vcf sample.normalized.bcf
```
