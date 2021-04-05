# phyloseq

To create a phyloseq object from a .qza file:

```
ps <- qza_to_phyloseq(
  features=features_file,
  tree=tree_file,
  taxonomy=taxonomy_file,
  metadata=metadata_file
)
```
