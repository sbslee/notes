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


Illumina
========

Adapter sequences
-----------------

Here's the `link <https://www.eurofinsgenomics.eu/media/1610545/illumina-adapter-sequences.pdf>`__ to Illumina's adapter sequences.
