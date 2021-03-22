# Genome Analysis Toolkit <a name="Genome-Analysis-Toolkit"></a>

```
I accessed to GATK Resource Bundle to download some files necessary for GATK's Base Quality Score Recalibration (BQSR).

				SEUNGs-MacBook-Pro-2:~ seungbeenlee$ ftp ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/
				ftp> cd /bundle/b37
				ftp> mget 1000G_phase1.indels.b37.*
				ftp> ls Mills_and_1000G_gold_standard.indels.b37.vcf*
```

## Process reference genome <a name="Process-reference-genome"></a>

According to this [post](https://gatk.broadinstitute.org/hc/en-us/articles/360035531652-FASTA-Reference-genome-format) by GATK Team:

> Most GATK tools additionally require that the main FASTA file be accompanied by a dictionary file ending in `.dict` and an index file ending in `.fai`, because it allows efficient random access to the reference bases. GATK will look for these index files based on their name, so it is important that they have the same basename as the FASTA file. If you do not have these files available for your organism's reference file, you can generate them very easily; instructions are included below.

To create to create a `.dict` file:

```
gatk CreateSequenceDictionary -R ref.fasta
```

To create a `.fai` file:

```
samtools faidx ref.fasta
```
