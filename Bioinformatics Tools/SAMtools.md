# SAMtools

To extract the sequence reads of a BAM file:

```samtools view sample.bam```

To extract the header of a BAM file:

```samtools view -H sample.bam```



```
View a BAM file.
$ samtools view IA_1.bam | less

Slice a BAM file.
$ samtools view -b CYP2D6.NWD100436.bam "chr22:42116498-42135906" > CYP2D6.sliced.NWD100436.bam

Merge BAM files.
$ samtools merge NWD100436.merged.bam CYP2D6.sliced.NWD100436.bam CYP2D7.sliced.NWD100436.bam

Get sample ID from BAM file.
$ samtools view -H $bam | grep "@RG" | head -1 | awk '{print $8}' | sed 's/SM://

Index a BAM file.
$ samtools index IA_1.bam

Index a FASTA file.
$ samtools faidx hs37d5.fa -o hs37d5.fa.fai

Index a VCF file.
$ bgzip -c file.vcf > file.vcf.gz # keeps the original file
$ tabix -p vcf file.vcf.gz

# How to index all BAM files in a directory:
module load samtools/latest
cd /net/grc/vol6/nobackup/nocleanup/sbslee/pgrn_seq_bams
for bam in `ls`; do samtools index $bam; done

# How to get the sample names from multiple BAM files:
for f in *.bam; do echo -ne "$f\t" ; samtools view -H $f | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq; done

# How to get the sample name from a BAM file:
samtools view -H test.bam | grep '^@RG' | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq

Find read lengths
samtools view /net/grc/vol5/data/processed/samples/119707/WHOLE_GENOME/119707.merged.sorted.markeddups.realigned.recal.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c
```
