# Genome Analysis Toolkit

## Table of Contents

* [Pipeline for germline short variant discovery](#Pipeline-for-germline-short-variant-discovery)
* [GATK resource bundle](#GATK-resource-bundle)
* [Process the reference genome](#Process-the-reference-genome)

## Pipeline for germline short variant discovery <a name="Pipeline-for-germline-short-variant-discovery"></a>

This pipeline is based on GATK Team's Best Practices Workflows for [Germline short variant discovery (SNPs + Indels)](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)

```
HaplotypeCaller
```

```
gatk GenomicsDBImport
```

## GATK resource bundle <a name="GATK-resource-bundle"></a>

The GATK resource bundle is a collection of standard files for working with human resequencing data with the GATK. For example, it can be used for Base Quality Score Recalibration (BQSR). See this [page](https://gatk.broadinstitute.org/hc/en-us/articles/360035890811-Resource-bundle) for more details.

**FTP server access was disabled on June 1, 2020.**

```
ftp ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/
ftp> cd /bundle/b37
ftp> mget 1000G_phase1.indels.b37.*
ftp> ls Mills_and_1000G_gold_standard.indels.b37.vcf*
```

## Process the reference genome <a name="Process-the-reference-genome"></a>

According to this [post](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format) by GATK Team: "Most GATK tools additionally require that the main FASTA file be accompanied by a dictionary file ending in `.dict` and an index file ending in `.fai`, because it allows efficient random access to the reference bases. GATK will look for these index files based on their name, so it is important that they have the same basename as the FASTA file. If you do not have these files available for your organism's reference file, you can generate them very easily; instructions are included below."

To create to create a `.dict` file:

```
gatk CreateSequenceDictionary -R ref.fasta
```

To create a `.fai` file:

```
samtools faidx ref.fasta
```
