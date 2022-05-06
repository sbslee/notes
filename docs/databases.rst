Databases
*********

Catalogue Of Somatic Mutations In Cancer (COSMIC)
=================================================

`COSMIC <https://cancer.sanger.ac.uk/cosmic>`__, the Catalogue Of Somatic Mutations In Cancer, is the world's largest and most comprehensive resource for exploring the impact of somatic mutations in human cancer.

Single Base Substitution (SBS) Signatures
-----------------------------------------

https://cancer.sanger.ac.uk/signatures/sbs/

Single base substitutions (SBS), also known as single nucleotide variants, are defined as a replacement of a certain nucleotide base. Considering the pyrimidines of the Watson-Crick base pairs, there are only six different possible substitutions: C>A, C>G, C>T, T>A, T>C, and T>G. These SBS classes can be further expanded considering the nucleotide context.

Current SBS signatures have been identified using 96 different contexts, considering not only the mutated base, but also the bases immediately 5’ and 3’.

cBioPortal
==========

https://www.cbioportal.org/

The cBioPortal for Cancer Genomics was originally developed at Memorial Sloan Kettering Cancer Center (MSK). The public cBioPortal site is hosted by the Center for Molecular Oncology at MSK. The cBioPortal software is now available under an open source license via GitHub. The software is now developed and maintained by a multi-institutional team, consisting of MSK, the Dana Farber Cancer Institute, Princess Margaret Cancer Centre in Toronto, Children's Hospital of Philadelphia, The Hyve in the Netherlands, and Bilkent University in Ankara, Turkey.

European Nucleotide Archive (ENA)
=================================

- `ENA: Guidelines and Tutorials <https://ena-docs.readthedocs.io/en/latest/>`__

Integrated Microbial Genomes (IMG)
==================================

The mission of the Integrated Microbial Genomes & Microbiomes(IMG/M) system is to support the annotation, analysis and distribution of microbial genome and microbiome datasets sequenced at DOE's Joint Genome Institute (JGI).

- `Official website <https://img.jgi.doe.gov/>`__

Ensembl
=======

This `page <http://asia.ensembl.org/info/website/archives/index.html>`__ says: "Ensembl aims to maintain stable identifiers for genes (ENSG), transcripts (ENST), proteins (ENSP) and exons (ENSE) as long as possible. Changes within the genome sequence assembly or an updated genome annotation may dramatically change a gene model. In these cases, the old set of stable IDs is retired and a new one assigned. Gene and transcript pages both have an ID History view which maps changes in the ID from the earliest version in Ensembl."

Variant Effect Predictor (VEP)
------------------------------

Order of annotations
^^^^^^^^^^^^^^^^^^^^

The ordering of the results per line simply uses the ENST IDs. For example:

- ENST00000572062
- ENST00000572573
- ENST00000572608
- ENST00000575820

Within a result, the consequences are ordered by severity. For example:

intron_variant&non_coding_transcript_variant


References:

  - `Order of annotation <https://github.com/Ensembl/ensembl-vep/issues/193>`__
  - `Ensembl Variation - Calculated variant consequences <https://m.ensembl.org/info/genome/variation/prediction/predicted_data.html>`__
  - `Cool stuff the Ensembl VEP can do: take your pick <https://www.ensembl.info/2019/03/22/cool-stuff-the-ensembl-vep-can-do-take-your-pick/>`__

Data Slicer
-----------

The `Data Slicer <http://grch37.ensembl.org/Homo_sapiens/Tools/DataSlicer?db=core;expand_form=true;tl=p4LmwgtfOgvfuAbL-7339566>`__ provides an interface which allows users to get subsections of either VCF (VCFtools) or BAM (SAMtools) files based on genomic coordinates.

References:

  - `Data Slicer <https://www.internationalgenome.org/data-slicer>`__
  - `How do I get a sub-section of a VCF file? <https://www.internationalgenome.org/faq/how-do-i-get-sub-section-vcf-file/>`__

MutSig
======

https://software.broadinstitute.org/cancer/cga/mutsig

MutSig stands for "Mutation Significance".  MutSig analyzes lists of mutations discovered in DNA sequencing, to identify genes that were mutated more often than expected by chance given background mutation processes.

UniProt
=======

https://www.uniprot.org/

The mission of UniProt is to provide the scientific community with a comprehensive, high-quality and freely accessible resource of protein sequence and functional information.
