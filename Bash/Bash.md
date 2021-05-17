# Bash

## Table of contents

* [Frequently used commands](#Frequently-used-commands)
    * [Checksum](#Checksum)
    * [List things](#List-things)
    * [Zipped files](#Zipped-files)
    * [Count things](#Count-things)
    * [Estimate size](#Estimate-size)
    * [Comparison](#Comparison)
    * [Text manipulation](#Text-manipulation)
    * [Check things](#Check-things)
* [Neat tricks](#Neat-tricks)
* [guration](#Bash-configuration)
* [System permission](#System-permission)
* [File transfer](#File-transfer)



## Neat tricks <a name="Neat-tricks"></a>

To access hard drives:

```
$ cd /
$ cd Volumes
$ cd ls
```

To copy a directory:

```
rsync -avzP source destination
```

The `-avzP` means:

* `a` - use archive mode
* `v` - be verbose
* `z` - compress file data during the transfer
* `P` - display progress and preserve partial files

To move the cursor forward by one word:

`esc` key + `F` key

To move the cursor backward by one word:

`esc` key + `B` key

To extract lines repeated at least three times:

```
awk '++a[$0] == 3 { print $0 }' example.txt
```

To print every fifth line:

```
awk 'NR % 5 == 0' example.txt
```

To skip the first two lines when priting a file:

```
tail -n +3 example.txt
```

To concatenate a string to each line of ls command output:

```
ls | xargs -i echo "Hello World {}"
```

To combine arrays as columns:

```
a=(A B C)
b=(1 2 3)
paste <(printf "%s\n" "${a[@]}") <(printf "%s\n" "${b[@]}")
```

To echo tab characters:

```
echo Hello$'\t'World
```

To read file names in the current directory into an array:

```
a=(*)
```


















## FUCs for Bash <a name="FUCs-for-Bash"></a>

```
# -- System ------------------------------------------------------------------

# Change group ownership.
chgrp -R shendure-pipeline *

ln -s original-fn new-fn

# -- Output files ------------------------------------------------------------

# Redirect stdout and stderr.
$ command > <out_file> 2><error_file>

# -- Vi and Vim --------------------------------------------------------------

# Find each occurrence of 'foo' (in all lines), and replace it with 'bar'.
:%s/foo/bar/g

# -- module ------------------------------------------------------------------

# List all available modules.
$ module list

# Load the latest version of <tool_name>
$ module load <tool_name>/latest

# Load all the available modules from the GS server.
$ module load modules modules-init modules-gs

$ module avail

# -- File transfer -----------------------------------------------------------

# Move only the files, not directories, from one directory to another.

$ find . -maxdepth 1 -type f -exec mv {} <directory_name> \;

rsync --progress -av /net/grc/vol5/nobackup/nocleanup/downsamplesForSteven ~/

scp sbslee@nexus.gs.washington.edu:/nfs/home/sbslee/liverbank_cyp2a7_rpkm.txt /Users/seungbeenlee/Desktop

wget -r --no-parent http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/b37.vcf/

[sbslee@grchead ~]$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR195/003/ERR1955323/ERR1955323_1.fastq.gz
--2018-10-25 10:31:50--  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR195/003/ERR1955323/ERR1955323_1.fastq.gz
           => “ERR1955323_1.fastq.gz”
Resolving ftp.sra.ebi.ac.uk... 193.62.192.7
Connecting to ftp.sra.ebi.ac.uk|193.62.192.7|:21... connected.
Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
==> TYPE I ... done.  ==> CWD (1) /vol1/fastq/ERR195/003/ERR1955323 ... done.
==> SIZE ERR1955323_1.fastq.gz ... 18111752809
==> PASV ... done.    ==> RETR ERR1955323_1.fastq.gz ... done.
Length: 18111752809 (17G) (unauthoritative)

 0% [                                                                                                             ] 32,134,016  4.33M/s  eta 92m 33s ^C
C

[sbslee@grchead ~]$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR195/003/ERR1955323/ERR1955323_2.fastq.gz
--2018-10-25 10:33:54--  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR195/003/ERR1955323/ERR1955323_2.fastq.gz
           => “ERR1955323_2.fastq.gz”
Resolving ftp.sra.ebi.ac.uk... 193.62.192.7
Connecting to ftp.sra.ebi.ac.uk|193.62.192.7|:21... connected.
Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
==> TYPE I ... done.  ==> CWD (1) /vol1/fastq/ERR195/003/ERR1955323 ... done.
==> SIZE ERR1955323_2.fastq.gz ... 24694164973
==> PASV ... done.    ==> RETR ERR1955323_2.fastq.gz ... done.
Length: 24694164973 (23G) (unauthoritative)

 0% [                                                                                                             ] 7,006,872   1.65M/s  eta 5h 0m   ^C
C



# -- Globus ------------------------------------------------------------------

[sbslee@grchead ~]$ wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR195/003/ERR1955323/ERR1955323.fastq.gz
--2018-10-25 10:31:09--  ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR195/003/ERR1955323/ERR1955323.fastq.gz
           => “ERR1955323.fastq.gz”
Resolving ftp.sra.ebi.ac.uk... 193.62.192.7
Connecting to ftp.sra.ebi.ac.uk|193.62.192.7|:21... connected.
Logging in as anonymous ... Logged in!
==> SYST ... done.    ==> PWD ... done.
==> TYPE I ... done.  ==> CWD (1) /vol1/fastq/ERR195/003/ERR1955323 ... done.
==> SIZE ERR1955323.fastq.gz ... 167227
==> PASV ... done.    ==> RETR ERR1955323.fastq.gz ... done.
Length: 167227 (163K) (unauthoritative)

100%[============================================================================================================>] 167,227      196K/s   in 0.8s    

2018-10-25 10:31:11 (196 KB/s) - “ERR1955323.fastq.gz” saved [167227]

# -- Text manipulation -------------------------------------------------------

awk '{print $1":"$2"-"$3}' pgrn_seq_nimblegen_968kb_v2.bed | sed 's/chr//g' > pgrnseq_v2_probes_all.list

# Turn ls output to space-separated.
ls /mnt/garnet/Users/sbslee/SNUH/30x/bam/*.bam | tr '\n' ' '

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

Related posts:

* [Search and replace](https://vim.fandom.com/wiki/Search_and_replace)




## File transfer <a name="File-transfer"></a>

```
lftp sftp://user_id:user_pw.@host_name:port_number
mirror -c target_dir destination_dir
```
