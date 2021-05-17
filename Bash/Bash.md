






# -- Tricks ------------------------------------------------------------------

for SAMPLE in `awk '{print $2}' trio.ped.lims.txt`
do
  echo /net/grc/vol5/data/processed/samples/$SAMPLE/WHOLE_GENOME/$SAMPLE.merged.sorted.markeddups.realigned.recal.bam
done > weiss_bam.list

for FRACTION in `cat fractions.txt`; do
  for SAMPLE in `cat samples.txt`; do
    echo /nfs/home/sbslee/projects/cyp2d6/20170405_sean_bam_downsampling/$FRACTION/$SAMPLE/$SAMPLE.merged.sorted.markeddups.realigned.recal.$FRACTION.sorted.bam;
  done > $FRACTION.list;
done

# Find a file.
for DIR in `cat weiss_sample_dir.list`;
do
  ls -lt $DIR | grep ".bai" | grep -v ".md5"
done

# Delete Everything After a Certain Pattern in a String
sed 's/WHOLE_GENOME.*/WHOLE_GENOME/' weiss_bam.list  	

# How to count unique lines in file
sort ips.txt | uniq -c | sort -bgr

# -- Awk ---------------------------------------------------------------------

awk 'NR==FNR{a[$0];next} !($0 in a)' weiss_bams_truncated.list weiss_bams_all.list > weiss_bams_good.list

```
