# phyloseq

## Table of Contents

* [Frequently used commands](#Frequently-used-commands)

## Frequently used commands <a name="Frequently-used-commands"></a>

To create a phyloseq object from a .qza file:

```
ps <- qza_to_phyloseq(
  features=features_file,
  tree=tree_file,
  taxonomy=taxonomy_file,
  metadata=metadata_file
)
```

To subset some samples that meet certain criteria:

```
ps <- subset_samples(ps, subset_exp)
```

To filter out taxa whose combined read count is less than 10:

```
ps <- filter_taxa(ps, sum(x) >= 10, TRUE)
```
