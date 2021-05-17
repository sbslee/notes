Bioinformatics
**************

Frequently used commands for Bioinformatics
===========================================

* To extract regions from a BED file:

    .. code-block:: console

        $ awk '{print $1":"$2"-"$3}' example.bed | sed 's/chr//g' > regions.list

* To zip a VCF file:

    .. code-block:: console

        $ bgzip -c sample.vcf > sample.vcf.gz

* To index a VCF file:

    .. code-block:: console

        $ tabix -p vcf sample.vcf.gz

* To slice a VCF file:

    .. code-block:: console

        $ tabix -h sample.vcf.gz chr1:1000-2000 > sliced.vcf.gz

* To slice a VCF file without using tabix:

    .. code-block:: console

        $ zcat sample.vcf.gz | awk '{OFS="\t"; if ($2 > 1000 && $2 < 2000){ print }}'

* To count the number of sequence reads in a FASTQ file:

    .. code-block:: console

        $ echo $(cat sample.fastq | wc -l) / 4 | bc

* To count the number of sequence reads in a zipped FASTQ file:

    .. code-block:: console

        $ echo $(zcat sample.fastq | wc -l) / 4 | bc

* To count the number of sequence reads in a zipped FASTQ file (macOS):

    .. code-block:: console

        $ echo $(zcat < sample.fastq | wc -l) / 4 | bc

* To extract only sequence reads from a zipped FASTQ file:

    .. code-block:: console

        $ zcat sample.fastq.gz | awk '{if (NR% 4 == 2) print $0}'

* To extract exon coordinates for a gene:

    .. code-block:: console

        $ grep -w "CYP2A6" Homo_sapiens.GRCh37.75.gtf | grep "CYP2A6-001" | grep -w "exon" | cut -f 1,4,5,9 -d$'\t' | cut -f 1,3 -d';' | sed 's/gene_id "ENSG00000255974"; //g'

* To extract rpkm values from a .gct file:

    .. code-block:: console

        $ printf "`echo $sample`\t`grep -w "CYP2A7" /net/grc/vol6/data/processed/samples/$sample/RNA_SEQ/qc/genes.rpkm.gct`\n"

GTCtoVCF
========

Here's the `link <https://github.com/Illumina/GTCtoVCF>`__ to the GTCtoVCF program.

gtc2vcf
=======

Here's the `link <https://github.com/freeseek/gtc2vcf>`__ to the gtc2vcf program.

bcl2fastq
=========

Introduction
------------

bcl2fastq is a software tool developed by Illumina Inc. for demultiplexing sequence read data. The official documentation is available `here <https://sapac.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf>`__.

+--------------------------+-------+----------------------+
| Platform                 | Lanes | Fluidically distinct |
+==========================+=======+======================+
| MiSeq                    | 1     | N/A                  |
+--------------------------+-------+----------------------+
| HiSeq - Rapid Mode       | 2     | Yes                  |
+--------------------------+-------+----------------------+
| HiSeq - High Output Mode | 8     | Yes                  |
+--------------------------+-------+----------------------+
| HiSeq X Ten              | 8     | Yes                  |
+--------------------------+-------+----------------------+
| NextSeq                  | 4     | No                   |
+--------------------------+-------+----------------------+
| NovaSeq S1, S2           | 2     | Yes with XP protocol |
+--------------------------+-------+----------------------+
| NovaSeq S3, S4           | 4     | Yes with XP protocol |
+--------------------------+-------+----------------------+

Commonly used options
---------------------

* ``--no-lane-splitting``

    Do not split FASTQ files by lane.

* ``--barcode-mismatches``

    | Specifies how to process each cycle:
    | * ``n`` - Ignore the cycle.
    | * ``Y`` (or ``y``) - Use the cycle.
    | * ``I`` - Use the cycle for an Index Read.
    | * A number - Repeat the previous character the indicated number of times.
    | * ``*`` - Repeat the previous character until the end of the read or index (length per ``RunInfo.xml``).
    | Commas separate read masks. The format for dual indexing is the following syntax or specified variations:
    | ``--use-bases-mask Y*,I*,I*,Y*``
    | You can also specify `--use-bases-mask` multiple times for separate lanes. In the following example, ``1:`` indicates that the setting applies to lane 1. The second ``--use-bases-mask`` parameter applies to all other lanes.
    | ``--use-bases-mask 1:y*,i*,i*,y* --use-bases-mask y*,n*,n*,y*``
    | If this option is not specified, ``RunInfo.xml`` determines the mask. If it cannot determine the mask, specify the `--use-bases-mask` option. When specified, the number of index cycles and the index length in the sample sheet must match.


* ``--tiles``

    | Selects a subset of available tiles for processing. To make multiple selections, separate the regular expressions with commas. For example:
    | To select all tiles ending with 5 in all lanes:
    | ``--tiles [0–9][0–9][0–9]5``
    | To select tile 2 in lane 1 and all the tiles in the other lanes:
    | ``--tiles s_1_0002,s_[2-8]``

Running
-------

**Case 1. MiSeq, 2x300 bp reads, dual indexing**

.. code-block:: console

    $ bcl2fastq \
      --output-dir $OUTPUT_DIR \
      --sample-sheet $SAMPLE_SHEET \
      --runfolder-dir $RUNFOLDER_DIR \
      --interop-dir $OUTPUT_DIR/Interop \
      --stats-dir $OUTPUT_DIR/Stats \
      --reports-dir $OUTPUT_DIR/Reports \
      --no-lane-splitting \
      --use-bases-mask Y301,I8,I8,Y301 \
      --barcode-mismatches 0 \
      --processing-threads 10


**Case 2. NextSeq, 2x150 bp reads, single indexing**

.. code-block:: console

    $ bcl2fastq \
      --output-dir $OUTPUT_DIR\
      --sample-sheet $SAMPLE_SHEET \
      --runfolder-dir $RUNFOLDER_DIR \
      --interop-dir $OUTPUT_DIR/Interop \
      --stats-dir $OUTPUT_DIR/Stats \
      --reports-dir $OUTPUT_DIR/Reports \
      --no-lane-splitting \
      --tiles s_1,s_2,s_3,s_4 \
      --use-bases-mask Y151,I8,Y151 \
      --barcode-mismatches 0 \
      --processing-threads 20

Cutadapt
========

Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads.

Trim Galore!
============

Trim Galore! is a wrapper tool around Cutadapt and FastQC to consistently apply quality and adapter trimming to FastQ files, with some extra functionality for MspI-digested RRBS-type (Reduced Representation Bisufite-Seq) libraries.

FastQC
======

FastQC is a quality control tool for high throughput sequence data.

SAMtools
========

Frequently used commands for SAMtools
-------------------------------------

* To extract the sequence reads of a BAM file:

    .. code-block:: console

        $ samtools view sample.bam

* To extract the header of a BAM file:

    .. code-block:: console

        $ samtools view -H sample.bam

* To index a BAM file:

    .. code-block:: console

        $ samtools index sample.bam

* To index a FASTA file:

    .. code-block:: console

        $ samtools faidx ref.fa -o ref.fa.fai

* To slice a BAM file:

    .. code-block:: console

        $ samtools view -b sample.bam "chr1:1000-2000" > sliced.bam

* To merge two BAM files:

    .. code-block:: console

        $ samtools merge merged.bam run1.bam run2.bam

* To extract the sample ID from a BAM file:

    .. code-block:: console

        $ samtools view -H sample.bam | grep "^@RG" | sed "s/.*SM:\([^\t]*\).*/\1/g" | uniq

* To estimate the read length of a BAM file:

    .. code-block:: console

        $ samtools view sample.bam | head -n 1000000 | cut -f 10 | perl -ne 'chomp;print length($_) . "\n"' | sort | uniq -c

BCFtools
========

Variant calling pipeline
------------------------

1. Calculate genotype likelihoods at each genomic position with coverage. Note that the reference FASTA file and the input BAM file(s) must have the same chromosome string style.

    .. code-block:: console

        $ bcftools mpileup -Ou -q 1 -a AD --max-depth 1000 -f ref.fa -r chr1:1000-2000 -o sample.bcf sample.bam

2. Make variant calls.

    .. code-block:: console

        $ bcftools call -Oz -mv -o sample.vcf.gz sample.bcf

3. Index the VCF file.

    .. code-block:: console

        $ bcftools index sample.vcf.gz

4. Left-align and normalize indels.

    .. code-block:: console

        $ bcftools norm -Ob -f ref.fa -o sample.normalized.bcf sample.vcf.gz

5. Filter variants.

    .. code-block:: console

        $ bcftools filter -Ov --IndelGap 5 -o sample.normalized.filtered.vcf sample.normalized.bcf

Illumina
========

Adapter sequences
-----------------

Here's the `link <https://www.eurofinsgenomics.eu/media/1610545/illumina-adapter-sequences.pdf>`__ to Illumina's adapter sequences.
