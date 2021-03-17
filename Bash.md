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

### Comparison

To find difference between two directories:

```
diff -qr dir1 dir2
```
