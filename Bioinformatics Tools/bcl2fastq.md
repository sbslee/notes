# bcl2fastq

## Table of contents

* [Introduction](#Introduction)
* [Platforms](#Platforms)
* [Commonly used options](#Commonly-used-options)
* [Running](#Running)

## Introduction <a name="Introduction"></a>

bcl2fastq is a software tool developed by Illumina Inc. for demultiplexing sequence read data. The official documentation is available [here](https://sapac.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/bcl2fastq/bcl2fastq2-v2-20-software-guide-15051736-03.pdf).

## Platforms <a name="Platforms"></a>

| Platform                 | Number of lanes | Fluidically distinct |
| -----------------------  | --------------- | -------------------- |
| MiSeq                    | 1               | N/A                  |
| HiSeq - Rapid Mode       | 2               | Yes                  |
| HiSeq - High Output Mode | 8               | Yes                  |
| HiSeq X Ten              | 8               | Yes                  |
| NextSeq                  | 4               | No                   |
| NovaSeq S1, S2           | 2               | Yes with XP protocol |
| NovaSeq S3, S4           | 4               | Yes with XP protocol |

## Commonly used options <a name="Commonly-used-options"></a>

* `--no-lane-splitting`

    Do not split FASTQ files by lane.

* `--barcode-mismatches`

    Specifies how to process each cycle:
    
    * n - Ignore the cycle.
    * Y (or y) - Use the cycle.
    * I - Use the cycle for an Index Read.
    * A number - Repeat the previous character the indicated number of times.
    * \* - Repeat the previous character until the end of the read or index (length per `RunInfo.xml`).
    
    Commas separate read masks. The format for dual indexing is the following syntax or specified variations:
    `--use-bases-mask Y*,I*,I*,Y*`

* `--tiles`

    Selects a subset of available tiles for processing. To make multiple selections, separate the regular expressions with commas.

    For example:
    
    To select all tiles ending with 5 in all lanes:
    
    `--tiles [0–9][0–9][0–9]5`
    
    To select tile 2 in lane 1 and all the tiles in the other lanes:
    
    `--tiles s_1_0002,s_[2-8]`

## Running <a name="Running"></a>

**Case 1. MiSeq, 2x300 bp reads, dual indexing**

```
bcl2fastq \
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
```

**Case 2. NextSeq, 2x150 bp reads, single indexing**

```
bcl2fastq \
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
```
