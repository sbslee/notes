# SnpEff and SnpSift

To download the pre-built human database (GRCh37.75):

```
java -jar snpEff.jar download -v GRCh37.75
```

To run annotation:

```
java -jar snpEff.jar eff hg19 in.vcf > ann.vcf
```
