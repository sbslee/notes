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
