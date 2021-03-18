# bcl2fastq

bcl2fastq is a software tool developed by Illumina Inc. for demultiplexing sequence read data.


| Platform                 | Number of lanes | Fluidically distinct |
| -----------------------  | --------------- | -------------------- |
| MiSeq                    | 1               | N/A                  |
| HiSeq - Rapid Mode       | 2               | Yes                  |
| HiSeq - High Output Mode | 8               | Yes                  |
| HiSeq X Ten              | 8               | Yes                  |
| NextSeq                  | 4               | No                   |
| NovaSeq S1, S2           | 2               | Yes with XP protocol |
| NovaSeq S3, S4           | 4               | Yes with XP protocol |

## Case 1. MiSeq, 2x300 bp reads, single indexing

```
bcl2fastq \
--output-dir $OUTPUT_DIR \
--sample-sheet $SAMPLE_SHEET \
--runfolder-dir $RUNFOLDER_DIR \
--interop-dir $OUTPUT_DIR/Interop \
--stats-dir $OUTPUT_DIR/Stats \
--reports-dir $OUTPUT_DIR/Reports \
--no-lane-splitting \
--tiles s_1 \
--use-bases-mask Y301,I8,I8,Y301 \
--barcode-mismatches 0 \
--processing-threads 10
```

## Case 2. NextSeq, 2x150 bp reads, dual indexing

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