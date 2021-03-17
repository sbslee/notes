# SAMtools

To extract the sequence reads of a BAM file:

```
samtools view sample.bam
```

To extract the header of a BAM file:

```
samtools view -H sample.bam
```

To index a BAM file:

```
samtools index sample.bam
```

To index a FASTA file:

```
samtools faidx ref.fa -o ref.fa.fai
```

To slice a BAM file:

```
samtools view -b sample.bam "chr1:1000-2000" > sliced.bam
```

To merge two BAM files:

```
samtools merge merged.bam run1.bam run2.bam
```

To extract the sample ID from a BAM file:

```
samtools view -H sample.bam | grep "^@RG" | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq
```

To estimate the read length of a BAM file:

```
samtools view sample.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c
```
