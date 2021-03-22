# Bash

## Table of Contents

* [Frequently used commands](#Frequently-used-commands)
    * [Checksum](#Checksum)
    * [List things](#List-things)
    * [Zipped files](#Zipped-files)
    * [Count things](#Count-things)
    * [Estimate size](#Estimate-size)
    * [Comparison](#Comparison)
* [Neat tricks](#Neat-tricks)
* [Bash configuration](#Bash-configuration)

## Frequently used commands <a name="Frequently-used-commands"></a>

### Checksum <a name="Checksum"></a>

To determine SHA-512 checksum:

```
shasum -a 256 example.txt
```

To determine MD5 checksum:

```
md5sum example.txt
```

To determine MD5 checksum (macOS):

```
md5 example.txt
```

### List things <a name="List-things"></a>

To list one file per line:

```
ls -1 dir
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

### Count things <a name="Count-things"></a>

To count unique lines in a file:

```
sort example.txt | uniq -c | sort -bgr
```

To count files in a directory:

```
find . | wc -l
```

### Estimate size <a name="Estimate-size"></a>

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

## Neat tricks <a name="Neat-tricks"></a>

To move the cursor forward by one word:

`esc` key + `F` key

To move the cursor backward by one word:

`esc` key + `B` key

To extract lines repeated at least three times:

```
awk '++a[$0] == 3 { print $0 }' example.txt
```

## Bash configuration <a name="Bash-configuration"></a>

The `.bashrc` file is used to provide a place where you can set up variables, functions and aliases, define your (PS1) prompt and define other settings that you want to use every time you open a new terminal window. The following command will activate the configuration:

```
source .bashrc
```

There is also the `.bash_profile` file, which is executed for login shells, while `.bashrc` is executed for interactive non-login shells. When an installed program cannot be called from the command line, add `export PATH=~/.local/bin:$PATH` to the `.bash_profile` file.
