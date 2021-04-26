# Genome Analysis Toolkit (GATK)

## Table of contents

* [Pipeline for germline short variant discovery](#Pipeline-for-germline-short-variant-discovery)
* [Pipeline for somatic short variant discovery](#Pipeline-for-somatic-short-variant-discovery)
* [GATK resource bundle](#GATK-resource-bundle)
* [Process the reference genome](#Process-the-reference-genome)
* [VCF filters](#VCF-filters)

## Pipeline for germline short variant discovery <a name="Pipeline-for-germline-short-variant-discovery"></a>

This pipeline is based on GATK Team's Best Practices Workflows for [Germline short variant discovery (SNPs + Indels)](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-).

### Call variants per-sample

```
gatk HaplotypeCaller \
-R ref.fa \
--emit-ref-confidence GVCF \
-I sample.bam \
-O sample.g.vcf
-L chr5:500-1000 \
--QUIET \
--java-options "-Xmx4G"
```

### Consolidate GVCFs

```
gatk GenomicsDBImport \
--intervals chr5:500-1000 \
--genomicsdb-workspace-path output_dir/temp/datastore \
--merge-input-intervals \
--QUIET \
--java-options "-Xmx4G" \
-V sample1.g.vcf \
-V sample2.g.vcf
```

### Joint-Call Cohort

```
gatk GenotypeGVCFs \
-R ref.fa \
-V gendb://output_dir/temp/datastore \
-O output_dir/temp/germline.joint.vcf \
--QUIET \
--java-options "-Xmx4G" \
-D dbsnp.vcf
```

### Filter variants

```
gatk VariantFiltration \
-R ref.fa \
-L chr5:500-1000 \
-O germline.joint.filtered.vcf \
--variant $output_dir/temp/germline.joint.vcf \
--filter-expression 'QUAL <= 50.0' \
--filter-name QUALFilter \
--QUIET \
--java-options "-Xmx4G"
```

## Pipeline for somatic short variant discovery <a name="Pipeline-for-somatic-short-variant-discovery"></a>

This pipeline is based on GATK Team's Best Practices Workflows for [Somatic short variant discovery (SNVs + Indels)](https://gatk.broadinstitute.org/hc/en-us/articles/360035894731).

### Tumor with matched normal

```
gatk Mutect2 \
-R reference.fa \
-I tumor.bam \
-I normal.bam \
-normal normal_sample_name \
--germline-resource af-only-gnomad.vcf.gz \
--panel-of-normals pon.vcf.gz \
-O somatic.vcf.gz
```

### Filter variants in a Mutect2 VCF callset

```
gatk FilterMutectCalls \
-R reference.fasta \
-V somatic.vcf.gz \
--contamination-table contamination.table \
--tumor-segmentation segments.tsv \
-O filtered.vcf.gz
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

## VCF filters <a name="VCF-filters"></a>

| Tool                    | ID               | Description                                                                                           |
| ------------------------| ---------------- | ----------------------------------------------------------------------------------------------------- |
| N/A                     | PASS             | All filters passed                                                                                    |
| N/A                     | FAIL             | Fail the site if all alleles fail but for different reasons.                                          |
| Mutect2                 | base_qual        | alt median base quality                                                                               |
| Mutect2                 | clustered_events | Clustered events observed in the tumor                                                                |
| Mutect2                 | contamination    | contamination                                                                                         |
| Mutect2                 | duplicate        | evidence for alt allele is overrepresented by apparent duplicates                                     |
| Mutect2                 | fragment         | abs(ref - alt) median fragment length                                                                 |
| Mutect2                 | germline         | Evidence indicates this site is germline, not somatic                                                 |
| Mutect2                 | haplotype        | Variant near filtered variant on same haplotype.                                                      |
| Mutect2                 | low_allele_frac  | Allele fraction is below specified threshold                                                          |
| Mutect2                 | map_qual         | ref - alt median mapping quality                                                                      |
| Mutect2                 | multiallelic     | Site filtered because too many alt alleles pass tumor LOD                                             |
| Mutect2                 | n_ratio          | Ratio of N to alt exceeds specified ratio                                                             |
| Mutect2                 | normal_artifact  | artifact_in_normal                                                                                    |
| Mutect2                 | orientation      | orientation bias detected by the orientation bias mixture model                                       |
| Mutect2                 | panel_of_normals | Blacklisted site in panel of normals                                                                  |
| Mutect2                 | position         | median distance of alt variants from end of reads                                                     |
| Mutect2                 | possible_numt    | Allele depth is below expected coverage of NuMT in autosome                                           |
| Mutect2                 | slippage         | Site filtered due to contraction of short tandem repeat region                                        |
| Mutect2                 | strand_bias      | Evidence for alt allele comes from one read direction only                                            |
| Mutect2                 | strict_strand    | Evidence for alt allele is not represented in both directions                                         |
| Mutect2                 | weak_evidence    | Mutation does not meet likelihood threshold                                                           |
| FilterByOrientationBias | orientation_bias | Orientation bias (in one of the specified artifact mode(s) or complement) seen in one or more samples |
