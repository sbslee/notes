# Genome Analysis Toolkit

## Table of Contents

* [Pipeline for germline short variant discovery](#Pipeline-for-germline-short-variant-discovery)
* [Pipeline for somatic short variant discovery](#Pipeline-for-somatic-short-variant-discovery)
* [GATK resource bundle](#GATK-resource-bundle)
* [Process the reference genome](#Process-the-reference-genome)
* [Mutect2 filters](#Mutect2-filters)

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

## Mutect2 filters <a name="Mutect2-filters"></a>

| ID     | Description        |
| ------ | ------------------ |
| PASS   | All filters passed |





##FILTER=<ID=PASS,Description="All filters passed">
##FILTER=<ID=FAIL,Description="Fail the site if all alleles fail but for different reasons.">
##FILTER=<ID=base_qual,Description="alt median base quality">
##FILTER=<ID=clustered_events,Description="Clustered events observed in the tumor">
##FILTER=<ID=contamination,Description="contamination">
##FILTER=<ID=duplicate,Description="evidence for alt allele is overrepresented by apparent duplicates">
##FILTER=<ID=fragment,Description="abs(ref - alt) median fragment length">
##FILTER=<ID=germline,Description="Evidence indicates this site is germline, not somatic">
##FILTER=<ID=haplotype,Description="Variant near filtered variant on same haplotype.">
##FILTER=<ID=low_allele_frac,Description="Allele fraction is below specified threshold">
##FILTER=<ID=map_qual,Description="ref - alt median mapping quality">
##FILTER=<ID=multiallelic,Description="Site filtered because too many alt alleles pass tumor LOD">
##FILTER=<ID=n_ratio,Description="Ratio of N to alt exceeds specified ratio">
##FILTER=<ID=normal_artifact,Description="artifact_in_normal">
##FILTER=<ID=orientation,Description="orientation bias detected by the orientation bias mixture model">
##FILTER=<ID=panel_of_normals,Description="Blacklisted site in panel of normals">
##FILTER=<ID=position,Description="median distance of alt variants from end of reads">
##FILTER=<ID=possible_numt,Description="Allele depth is below expected coverage of NuMT in autosome">
##FILTER=<ID=slippage,Description="Site filtered due to contraction of short tandem repeat region">
##FILTER=<ID=strand_bias,Description="Evidence for alt allele comes from one read direction only">
##FILTER=<ID=strict_strand,Description="Evidence for alt allele is not represented in both directions">
##FILTER=<ID=weak_evidence,Description="Mutation does not meet likelihood threshold">
