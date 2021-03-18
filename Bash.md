# Bash

## Table of Contents

* [Frequently used commands](#Frequently-used-commands)
    * [Checksum](#Checksum)
    * [Zipped files](#Zipped-files)
    * [Count things or estimate size](#Count-things-or-estimate-size)
    * [Comparison](#Comparison)

## Frequently used commands <a name="Frequently-used-commands"></a>

### Checksum <a name="Checksum"></a>

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

### Zipped files <a name="Zipped-files"></a>

To create a .tar.gz file:

```
tar -czvf dir.tar.gz dir
```

To unzip a .tar.gz file:

```
tar -xf dir.tar.gz
```

### Count things or estimate size <a name="Count-things-or-estimate-size"></a>

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
du -sh dir
```

### Comparison <a name="Comparison"></a>

To find difference between two directories:

```
diff -qr dir1 dir2
```
