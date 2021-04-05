# phyloseq

* [Frequently used commands](#Frequently-used-commands)

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
ps <- subset_samples(ps, subsetting_expression)
```
