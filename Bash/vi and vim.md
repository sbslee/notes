# vi and vim

To search a pattern:

* Press `/`.
* Type the search pattern.
* Press `Enter` to perform the search.
* Press `n` to find the next occurrence or `N` to find the previous occurrence.

To search and replace in the entire file:

```
:%s/foo/bar/g
```

To search and replace a pattern involving the `/` character:

```
:%s#/foo#/bar#g
```

To move the cursor to end of the file:

Press the `Esc` key and then press the `Shift` and `G` keys together
