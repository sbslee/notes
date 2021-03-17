# Bash

## Frequently used commands

### Checksum

To determine SHA-512 checksum:

```
shasum -a 256 test.txt
```

To determine MD5 checksum:

```
md5sum test.txt
```

To determine MD5 checksum (macOS):

```
md5 test.txt
```

### Zipped files

To create a .tar.gz file:

```
tar -czvf dir.tar.gz dir
```

To unzip a .tar.gz file:

```
tar -xf dir.tar.gz
```

### Count things or estimate size

To count unique lines in a file:

```
sort test.txt | uniq -c | sort -bgr
```

To count files in a directory:

```
find . | wc -l
```

To estimate storage size:

```
df -h
```

To estimate directory size:

```
$ du -sh dir
```

### Comparison

To find difference between two directories:

```
diff -qr dir1 dir2
```
